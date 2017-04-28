#!/usr/bin/env python

import os
import sys
from argparse import ArgumentParser
from glob import glob
from re import findall
from os.path import splitext
from netCDF4 import Dataset
from averager import MeanAverager
from filespecs import AggregationFile
from aggmaskloader import AggMaskLoader
from numpy.ma import masked_array, masked_where
from numpy import zeros, ones, mod, where, isnan, logical_not, logical_and, arange
from warnings import filterwarnings


def findfile(filelist, scenario, variable):
    for ncfile in filelist:
        if '_%s_%s_' % (scenario, variable) in ncfile:
            return ncfile
    return None


def shiftdata(data, pl_date, harv_date, year=1):
    latidx, lonidx = where(logical_and(pl_date >= harv_date, harv_date >= year))
    datac = data.copy()
    shiftd = datac[:, latidx, lonidx]
    shiftd[1:] = shiftd[: -1]  # shift data down one
    shiftd[0].mask = True
    datac[:, latidx, lonidx] = shiftd
    return datac

# parse inputs
parser = ArgumentParser("GGCMI Aggregator")
parser.add_argument("--adaptation", required=True,
                    help="Adaptation level (0=standard simulation, 1=full growing season adaptation)")
parser.add_argument("--agg", required=True, help="Aggregation mask file:var")
parser.add_argument("--calc_area", action="store_true", default=False,
                    help="Flag to indicate weights are fractions")
parser.add_argument("--co2", required=True, help="co2 level (C360, C510, C660, C810)")
parser.add_argument("--crop", required=True, help="Crop")
parser.add_argument("--gsfile", required=True, help="Growing season file")
parser.add_argument("--indir", required=True, help="Input directory")
parser.add_argument("--lufile", required=True, help="Landuse weight file")
parser.add_argument("--outfile", required=True, help="Output file")
parser.add_argument("--nitrogen", required=True, help="Nitrogen level (N10, N60, N200)")
parser.add_argument("--precipitation", required=True,
                    help="Precipitation level (W-50, W-30, W-20, W-10, W0, W10, W20, W30, Winf)")
parser.add_argument("--temperature", required=True, help="Temperature level (T-1, T0, T1, T2, T3, T4, T6)")
args = parser.parse_args()

adaptation = args.adaptation
agg = args.agg
calc_area = args.calc_area
co2 = args.co2
crop = args.crop
gsfile = args.gsfile
indir = args.indir
lufile = args.lufile
nitrogen = args.nitrogen
outfile = args.outfile
precipitation = args.precipitation
temperature = args.temperature
crop_names = {'mai': 'maize', 'soy': 'soy', 'swh': 'spring_wheat', 'wwh': 'winter_wheat', 'ric': 'rice'}
crop_long = crop_names[crop]

filterwarnings('ignore')

# Winf uses irrigated, everything else uses rainfed
if precipitation == "Winf":
    irrigation = "firr"
else:
    irrigation = "noirr"

# some constants
yieldthr1 = 0.1
yieldthr2 = 30
yearthr = 1

files = glob(os.path.join(indir, crop_long, adaptation, "*/*%s_%s_%s_%s*[!old]" %
                          (co2, temperature, precipitation, nitrogen)))
if not files:
    sys.exit("No files found")

# extract time and units from first filename
y1, y2 = [int(y) for y in findall(r'\d+', splitext(files[0])[0])[-2:]]
years = arange(y1, y2 + 1)
tunits = 'growing seasons since %d-01-01 00:00:00' % y1

# check time against time in file
with Dataset(files[0]) as f:
    tvar = f.variables['time']
    if len(tvar) != len(years):
        tunits = tvar.units
        years = tvar[:] + int(findall(r'\d+', tunits)[0]) - 1
        y1 = years.min()
        y2 = years.max()

# load landuse mask
with Dataset(lufile) as f:
    lats = f.variables['lat'][:]
    lons = f.variables['lon'][:]
    wir = f.variables['irrigated'][:]
    wrf = f.variables['rainfed'][:]

    weights = masked_array(zeros((2,) + wir.shape), mask=ones((2,) + wir.shape))
    weights[0] = wir
    weights[1] = wrf

    if 'time' in f.variables:
        ltime = f.variables['time'][:]
        ltunits = f.variables['time'].units
        ltime += int(findall(r'\d+', ltunits)[0])
        weights = weights[:, logical_and(ltime >= y1, ltime <= y2)]  # restrict range
    else:
        ltime = years.copy()

# restrict to overlapping years
yrsinrange = logical_and(years >= min(ltime), years <= max(ltime))
years = years[yrsinrange]
t0 = min(years) - y1 + 1
time = arange(t0, t0 + len(years))

# load aggregation mask
afile, avar = [a.strip() for a in agg.split(':')]
aggloader = AggMaskLoader(afile, avar)
adata = aggloader.data()[0]
audata = aggloader.udata()[0]
aname = aggloader.names()[0]
aunits = aggloader.units()[0]
alongname = aggloader.longnames()[0]

# load growing season file
with Dataset(gsfile) as f:
    pdate = f.variables['planting_day'][:]
    hdate = f.variables['growing_season_length'][:]

# get variables and scenarios
variables = []
scens = []
scens_full = []

for i in range(len(files)):
    file_split = files[i].split('_')
    variables.append(file_split[3])
    scens.append(file_split[2])
    scens_full.append(file_split[2] + '_' + irrigation)

variables = list(set(variables))
scens = list(set(scens))
scens_full = list(set(scens_full))

variables.sort()
scens.sort()
scens_full.sort()

# number of variables, times, scenarios, aggregation levels
nvariables = len(variables)
ntimes = len(time)
nscens = len(scens)
naudata = len(audata)

# preallocate final averages and areas
averages = masked_array(zeros((nvariables, naudata, ntimes, nscens, 3)),
                        mask=ones((nvariables, naudata, ntimes, nscens, 3)))
areas = masked_array(zeros((naudata, ntimes, nscens, 3)), mask=ones((naudata, ntimes, nscens, 3)))

# aggregator object
avobj = MeanAverager()

vunits = [''] * nvariables
for i in range(len(scens)):
    scen = scens[i]

    # find scenario and irr
    sidx = scens.index(scen)
    iidx = int(irrigation != 'firr')

    # Load planting file
    plantingfile = findfile(files, scen, 'plant-day')
    if plantingfile:
        with Dataset(plantingfile) as f:
            planting_date = f.variables['plant-day_' + crop][0]
    else:
        planting_date = pdate[iidx]

    # load harvest file
    harvestfile = findfile(files, scen, 'maty-day')
    if harvestfile:
        with Dataset(harvestfile) as f:
            harvest_date = f.variables['maty-day_' + crop][0]
    else:
        harvest_date = hdate[iidx]

    # convert to Julian day
    harvest_date = mod(planting_date + harvest_date, 366)
    harvest_date[harvest_date == 0] = 1

    # load yield file
    yieldfile = findfile(files, scen, 'yield')

    if not yieldfile:
        yvar = masked_array(zeros((ntimes,) + planting_date.shape), mask=zeros((ntimes,) + planting_date.shape))
    else:
        with Dataset(yieldfile) as f:
            yvar = f.variables['yield_' + crop][:]
            yvar = shiftdata(yvar, planting_date, harvest_date, yearthr)     # shift data
            yvar = masked_where(yvar < yieldthr1, yvar)  # mask yields below threshold
            yvar = masked_where(yvar > yieldthr2, yvar)  # mask yields above threshold
            yvar = yvar[yrsinrange]                      # restrict to overlapping years
    ymsk = logical_not(yvar.mask)

    # compute areas
    print ymsk.sum(), ymsk.shape, adata.shape, weights[iidx].shape
    areas[:, :, sidx, iidx] = avobj.areas(yvar, adata, lats, weights[iidx], calc_area)

    for j in range(nvariables):
        varfile = findfile(files, scen, variables[j])
        print "file %s" % varfile
        if not varfile:
            continue
        with Dataset(varfile) as f:
            var = f.variables[variables[j] + '_' + crop]
            if not i:
                vunits[j] = var.units if 'units' in var.ncattrs() else ''
            var = var[:]
        var[isnan(var)] = 0.                   # change NaNs to zero
        var = shiftdata(var, planting_date, harvest_date, yearthr)  # shift data
        var = var[yrsinrange]                  # restrict to overlapping years

        # aggregate
        print "sum is %s" % avobj.sum(var, adata, lats, weights[iidx], calc_area, ymsk)
        averages[j, :, :, sidx, iidx] = avobj.sum(var, adata, lats, weights[iidx], calc_area, ymsk)
        averages[j, :, :, sidx, iidx] /= areas[:, :, sidx, iidx]

# sum areas
area1, area2 = areas[:, :, :, 0], areas[:, :, :, 1]
area1[logical_and(area1.mask, ~area2.mask)] = 0.
area2[logical_and(area2.mask, ~area1.mask)] = 0.
areas[:, :, :, 2] = area1 + area2

# create output file
fout = AggregationFile(outfile, time, tunits, scens, audata, aname, aunits, alongname)

# add area
fout.append('area_' + aname, areas, (aname, 'time', 'scen', 'irr'), 'hectares', aname + ' harvested area')

# add variables
for i in range(nvariables):
    # average
    av1 = areas[:, :, :, 0] * averages[i, :, :, :, 0]
    av2 = areas[:, :, :, 1] * averages[i, :, :, :, 1]
    av1[logical_and(av1.mask, ~av2.mask)] = 0.
    av2[logical_and(av2.mask, ~av1.mask)] = 0.
    averages[i, :, :, :, 2] = (av1 + av2) / areas[:, :, :, 2]
    fout.append(variables[i] + '_' + aname, averages[i], (aname, 'time', 'scen', 'irr'), vunits[i],
                'average ' + aname + ' ' + variables[i])
