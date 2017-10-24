#!/usr/bin/env python

import sys
from argparse import ArgumentParser
from netCDF4 import Dataset
from os import listdir
from re import findall
from os.path import sep, isfile, basename
from numpy.ma import masked_array, masked_where
from numpy import zeros, ones, mod, where, isnan, logical_not, logical_and, arange
from aggmaskloader import AggMaskLoader
from averager import MeanAverager
from filespecs import AggregationFile
import warnings
warnings.filterwarnings('ignore')


def filterfiles(listing):
    files = []
    for l in listing:
        # skip old files
        if isfile(l) and not l.endswith('.old'):
            files.append(l)
    return files


def findfile(files, scen_irr, var):
    for f in files:
        if '_%s_%s_' % (scen_irr, var) in f:
            return f
    return []


def shiftdata(data, pdate, hdate, yearthr=1):
    latidx, lonidx = where(logical_and(pdate >= hdate, hdate >= yearthr))
    datac = data.copy()
    shiftd = datac[:, latidx, lonidx]
    # shift data down one
    shiftd[1:] = shiftd[: -1]
    shiftd[0].mask = True
    datac[:, latidx, lonidx] = shiftd
    return datac


def main():
    parser = ArgumentParser()
    parser.add_argument("-i", "--indir", default="pDSSAT/AgCFSR/maize", help="Input directory")
    parser.add_argument("-c", "--crop", default="mai", help="Crop")
    parser.add_argument("-l", "--lufile", default="maize.nc4", help="Landuse weight file")
    parser.add_argument("-a", "--agg", default="", help="Aggregation mask file:var")
    parser.add_argument("-g", "--gsfile", default="maize_growing_season_dates.nc4", help="Growing season file")
    parser.add_argument("--calc_area", action="store_true", dest="calcarea", default=False,
                        help="Flag to indicate weights are fractions (optional)")
    parser.add_argument("-o", "--outfile", default="", help="Output file")
    parser.add_argument("-y", "--year_start", type=int, required=True, help="Start year")
    args = parser.parse_args()

    indir = args.indir
    crop = args.crop
    lufile = args.lufile
    agg = args.agg
    gsfile = args.gsfile
    calcarea = args.calcarea
    outfile = args.outfile
    year_start = args.year_start

    yieldthr1 = 0.1
    yieldthr2 = 30
    yearthr = 1

    # files in directory
    listing = [indir + sep + l for l in listdir(indir)]
    files = filterfiles(listing)
    files = [basename(l) for l in files]
    files = [f for f in files if '_' + crop + '_' in f]

    if not len(files):
        sys.exit('No files found. Skipping directory . . .')

    # check time against time in file
    with Dataset(indir + sep + files[0]) as f:
        tvar = f.variables['time']
        tunits = tvar.units
        years = tvar[:] + int(findall(r'\d+', tunits)[0]) - 1
        y1 = years.min()
        y2 = years.max()
        # If the first year in the data is not the first year of the climate data, don't process
        if y1 != year_start:
            print "Data/Climate time mismatch in %s" % files[0]
            sys.exit(0)

    # load landuse mask
    with Dataset(lufile) as f:
        lats, lons = f.variables['lat'][:], f.variables['lon'][:]

        wir = f.variables['irrigated'][:]
        wrf = f.variables['rainfed'][:]

        weights = masked_array(zeros((2,) + wir.shape), mask=ones((2,) + wir.shape))
        weights[0] = wir
        weights[1] = wrf

        if 'time' in f.variables:
            ltime = f.variables['time'][:]
            ltunits = f.variables['time'].units
            ltime += int(findall(r'\d+', ltunits)[0])
            # restrict range
            weights = weights[:, logical_and(ltime >= y1, ltime <= y2)]
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
        elif len(fs) == 10:  # pt files, etc.
            variables.append(fs[5])
            scens.append(fs[3] + '_' + fs[4])
            scens_full.append('_'.join([fs[3], ir, fs[4]]))
        else:
            continue

    variables = list(set(variables))
    scens = list(set(scens))
    scens_full = list(set(scens_full))

    variables.sort()
    scens.sort()
    scens_full.sort()

    # number of variables, times, scenarios, aggregation levels
    nv = len(variables)
    nt = len(time)
    ns = len(scens)
    na = len(audata)

    # preallocate final averages and areas
    averages = masked_array(zeros((nv, na, nt, ns, 3)), mask=ones((nv, na, nt, ns, 3)))
    areas = masked_array(zeros((na, nt, ns, 3)), mask=ones((na, nt, ns, 3)))

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
            sidx = scens.index('_'.join(scen_irr_split[0::2]))
        iidx = int(scen_irr_split[1] != 'firr')

        # load planting file
        plantingfile = findfile(files, scen_irr, 'plant-day')
        if plantingfile:
            with Dataset(indir + sep + plantingfile) as f:
                pd = f.variables['plant-day_' + crop][0]
        else:
            pd = pdate[iidx]

        # load harvest file
        harvestfile = findfile(files, scen_irr, 'maty-day')
        if harvestfile:
            with Dataset(indir + sep + harvestfile) as f:
                hd = f.variables['maty-day_' + crop][0]
        else:
            hd = hdate[iidx]

        # convert to Julian day
        hd = mod(pd + hd, 366)
        hd[hd == 0] = 1

        # load yield file
        yieldfile = findfile(files, scen_irr, 'yield')
        if not yieldfile:
            yvar = masked_array(zeros((nt,) + pd.shape), mask=zeros((nt,) + pd.shape))
        else:
            with Dataset(indir + sep + yieldfile) as f:
                yvar = f.variables['yield_' + crop][:]
                yvar = shiftdata(yvar, pd, hd, yearthr)      # shift data
                yvar = masked_where(yvar < yieldthr1, yvar)  # mask yields below threshold
                yvar = masked_where(yvar > yieldthr2, yvar)  # mask yields above threshold
                yvar = yvar[yrsinrange]                      # restrict to overlapping years
        ymsk = logical_not(yvar.mask)

        # compute areas
        print ymsk.sum(), ymsk.shape, adata.shape, weights[iidx].shape
        areas[:, :, sidx, iidx] = avobj.areas(yvar, adata, lats, weights[iidx], calcarea)

        for j in range(nv):
            # load variable
            varfile = findfile(files, scen_irr, variables[j])
            if not len(varfile):
                continue
            with Dataset(indir + sep + varfile) as f:
                var = f.variables[variables[j] + '_' + crop]
                if not i:
                    vunits[j] = var.units if 'units' in var.ncattrs() else ''
                var = var[:]
            var[isnan(var)] = 0.                   # change NaNs to zero
            var = shiftdata(var, pd, hd, yearthr)  # shift data
            var = var[yrsinrange]                  # restrict to overlapping years

            # aggregate
            averages[j, :, :, sidx, iidx] = avobj.sum(var, adata, lats, weights[iidx], calcarea, ymsk)
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

        fout.append(variables[i] + '_' + aname, averages[i], (aname, 'time', 'scen', 'irr'), vunits[i], 'average ' +
                    aname + ' ' + variables[i])


if __name__ == "__main__":
    main()

