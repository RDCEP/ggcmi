#!/usr/bin/env python

# add paths
import os, sys
for p in os.environ['PATH'].split(':'): sys.path.append(p)

# import modules
from re import findall
from os import listdir
from os.path import splitext
from optparse import OptionParser
from netCDF4 import Dataset as nc
from averager import MeanAverager
from filespecs import AggregationFile
from aggmaskloader import AggMaskLoader
from os.path import sep, isfile, basename
from numpy.ma import masked_array, masked_where
from numpy import zeros, ones, mod, where, isnan, logical_not, logical_and, append, newaxis, arange

def filterfiles(listing):
    files = []
    for l in listing:
        if isfile(l) and not l.endswith('.old'): # skip old files
            files.append(l)
    return files

def findfile(files, scen_irr, var):
    for f in files:
        if '_%s_%s_' % (scen_irr, var) in f:
            return f
    return []

def shiftdata(data, pdate, hdate, yearthr = 1):
    latidx, lonidx = where(logical_and(pdate >= hdate, hdate >= yearthr))
    datac = data.copy()
    shiftd = datac[:, latidx, lonidx]
    shiftd[1 :] = shiftd[: -1] # shift data down one
    shiftd[0].mask = True
    datac[:, latidx, lonidx] = shiftd
    return datac

# parse inputs
parser = OptionParser()
parser.add_option("-i", "--indir", dest = "indir", default = "pDSSAT/AgCFSR/maize", type = "string",
                  help = "Input directory")
parser.add_option("-c", "--crop", dest = "crop", default = "mai", type = "string",
                  help = "Crop")
parser.add_option("-l", "--lufile", dest = "lufile", default = "maize.nc4", type = "string",
                  help = "Landuse weight file", metavar = "FILE")
parser.add_option("-a", "--agg", dest = "agg", default = "", type = "string",
                  help = "Aggregation mask file:var")
parser.add_option("-g", "--gsfile", dest = "gsfile", default = "maize_growing_season_dates.nc4", type = "string",
                  help = "Growing season file", metavar = "FILE")
parser.add_option("--calc_area", action = "store_true", dest = "calcarea", default = False,
                  help = "Flag to indicate weights are fractions (optional)")
parser.add_option("-o", "--outfile", dest = "outfile", default = "", type = "string",
                  help = "Output file", metavar = "FILE")
options, args = parser.parse_args()

indir    = options.indir
crop     = options.crop
lufile   = options.lufile
agg      = options.agg
gsfile   = options.gsfile
calcarea = options.calcarea
outfile  = options.outfile

# some constants
yieldthr1, yieldthr2, yearthr = 0.1, 30, 1

# files in directory
listing = [indir + sep + l for l in listdir(indir)]
files   = filterfiles(listing)
files   = [basename(l) for l in files]
files   = [f for f in files if '_' + crop + '_' in f]

if not len(files):
    print 'No files found. Skipping directory . . .'
    exit(0)

# extract time and units from first filename
y1, y2 = [int(y) for y in findall(r'\d+', splitext(files[0])[0])[-2 :]]
time   = arange(1, y2 - y1 + 2)
years  = arange(y1, y2 + 1)
tunits = 'growing seasons since %d-01-01 00:00:00' % y1

# load landuse mask
with nc(lufile) as f:
    lats, lons = f.variables['lat'][:], f.variables['lon'][:]

    weights = f.variables['irrigated'][:][newaxis]
    weights = append(weights, f.variables['rainfed'][:][newaxis], axis = 0)

    if 'time' in f.variables:
        ltime   = f.variables['time'][:]
        ltunits = f.variables['time'].units
        ltime  += int(findall(r'\d+', ltunits)[0])
        weights = weights[:, logical_and(ltime >= y1, ltime <= y2)] # restrict range
    else:
        ltime = years.copy()

# restrict to overlapping years
yrsinrange = logical_and(years >= min(ltime), years <= max(ltime))
years      = years[yrsinrange]
t0         = min(years) - y1 + 1
time       = arange(t0, t0 + len(years))

# load aggregation mask
afile, avar = [a.strip() for a in agg.split(':')]
aggloader   = AggMaskLoader(afile, avar)
adata       = aggloader.data()[0]
audata      = aggloader.udata()[0]
aname       = aggloader.names()[0]
aunits      = aggloader.units()[0]
alongname   = aggloader.longnames()[0]

# load growing season file
with nc(gsfile) as f:
    pdate = f.variables['planting_day'][:]
    hdate = f.variables['growing_season_length'][:]

# get variables and scenarios
variables = []; scens = []; scens_full = []
for i in range(len(files)):
    fs = files[i].split('_')

    if 'noirr' in fs:
        fs.remove('noirr')
        ir = 'noirr'
    else:
        fs.remove('firr')
        ir = 'firr'

    if len(fs) == 9:
        variables.append(fs[4])
        scens.append(fs[3])
        scens_full.append(fs[3] + '_' + ir)
    elif len(fs) == 10: # pt files, etc.
        variables.append(fs[5])
        scens.append(fs[3] + '_' + fs[4])
        scens_full.append('_'.join([fs[3], ir, fs[4]]))
    else:
	continue # irregular files

variables  = list(set(variables));  variables.sort()
scens      = list(set(scens));      scens.sort()
scens_full = list(set(scens_full)); scens_full.sort()

# number of variables, times, scenarios, aggregation levels
nv, nt, ns, na = len(variables), len(time), len(scens), len(audata)

# preallocate final averages and areas
averages = masked_array(zeros((nv, na, nt, ns, 3)), mask = ones((nv, na, nt, ns, 3)))
areas    = masked_array(zeros((na, nt, ns, 3)),     mask = ones((na, nt, ns, 3)))

# aggregator object
avobj = MeanAverager()

vunits = [''] * nv
for i in range(len(scens_full)):
    scen_irr = scens_full[i]

    # find scenario and irr
    scen_irr_split = scen_irr.split('_')
    if len(scen_irr_split) == 2:
        sidx = scens.index(scen_irr_split[0])
    else:
        sidx = scens.index('_'.join(scen_irr_split[0 :: 2]))
    iidx = int(scen_irr_split[1] != 'firr')

    # load planting file
    plantingfile = findfile(files, scen_irr, 'plant-day')
    if plantingfile:
        with nc(indir + sep + plantingfile) as f:
            pd = f.variables['plant-day_' + crop][0]
    else:
        pd = pdate[iidx]

    # load harvest file
    harvestfile = findfile(files, scen_irr, 'maty-day')
    if harvestfile:
        with nc(indir + sep + harvestfile) as f:
            hd = f.variables['maty-day_' + crop][0]
    else:
        hd = hdate[iidx]

    # convert to Julian day
    hd = mod(pd + hd, 366)
    hd[hd == 0] = 1

    # load yield file
    yieldfile = findfile(files, scen_irr, 'yield')
    if yieldfile == []:
        yvar = masked_array(zeros((nt,) + pd.shape), mask = zeros((nt,) + pd.shape))
    else:
        with nc(indir + sep + yieldfile) as f:
            yvar = f.variables['yield_' + crop][:]
            yvar = shiftdata(yvar, pd, hd, yearthr)     # shift data
            yvar = masked_where(yvar < yieldthr1, yvar) # mask yields below threshold
            yvar = masked_where(yvar > yieldthr2, yvar) # mask yields above threshold
            yvar = yvar[yrsinrange]                     # restrict to overlapping years
    ymsk = logical_not(yvar.mask)

    # compute areas
    areas[:, :, sidx, iidx] = avobj.areas(yvar, adata, lats, weights[iidx], calcarea)

    for j in range(nv):
        # load variable
        varfile = findfile(files, scen_irr, variables[j])
        if not len(varfile): continue
        with nc(indir + sep + varfile) as f:
            var = f.variables[variables[j] + '_' + crop]
            if not i: vunits[j] = var.units if 'units' in var.ncattrs() else ''
            var = var[:]
        var[isnan(var)] = 0.                  # change NaNs to zero
        var = shiftdata(var, pd, hd, yearthr) # shift data
        var = var[yrsinrange]                 # restrict to overlapping years

        # aggregate
        averages[j, :, :, sidx, iidx]  = avobj.sum(var, adata, lats, weights[iidx], calcarea, ymsk)
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
for i in range(nv):
    # average
    av1 = areas[:, :, :, 0] * averages[i, :, :, :, 0]
    av2 = areas[:, :, :, 1] * averages[i, :, :, :, 1]
    av1[logical_and(av1.mask, ~av2.mask)] = 0.
    av2[logical_and(av2.mask, ~av1.mask)] = 0.
    averages[i, :, :, :, 2] = (av1 + av2) / areas[:, :, :, 2]

    fout.append(variables[i] + '_' + aname, averages[i], (aname, 'time', 'scen', 'irr'), vunits[i], 'average ' + aname + ' ' + variables[i])
