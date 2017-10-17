#!/usr/bin/env python

from argparse import ArgumentParser
from netCDF4 import Dataset
from numpy import setdiff1d, unique, resize, where, zeros, diff, ma
from shutil import copyfile

parser = ArgumentParser()
parser.add_argument("-i", "--inputfile", default="", help="Input file", required=True)
parser.add_argument("-o", "--outputfile", default="", help="Output file", required=True)
args = parser.parse_args()

inputfile = args.inputfile
outputfile = args.outputfile

copyfile(inputfile, outputfile)

with Dataset(outputfile, 'a') as f:
    lats = f.variables['lat'][:]
    lons = f.variables['lon'][:]
    dlat = abs(diff(lats)[0])
    dlon = abs(diff(lons)[0])
    latd = resize(lats, (len(lons), len(lats))).T
    lond = resize(lons, (len(lats), len(lons)))
    levels = [str(s) for s in setdiff1d(f.variables.keys(), ['lat', 'lon'])]

    f.createDimension('bnds', 2)
    bndvar = f.createVariable('bnds', 'i4', 'bnds')
    bndvar[:] = [1, 2]
    bndvar.units = 'mapping'
    bndvar.long_name = 'lower, upper'

    for lev in levels:
        levg = f.variables[lev][:]

        levvals = unique(levg)
        if ma.isMaskedArray(levvals):
            levvals = levvals[~levvals.mask]

        nlevs = len(levvals)
        latbnds = zeros((nlevs, 2))
        lonbnds = zeros((nlevs, 2))
        for i in range(nlevs):
            latidx, lonidx = where(levg == levvals[i])
            latlev = latd[latidx, lonidx]
            lonlev = lond[latidx, lonidx]
            latbnds[i] = [latlev.min() - dlat / 2., latlev.max() + dlat / 2.]
            lonbnds[i] = [lonlev.min() - dlon / 2., lonlev.max() + dlon / 2.]

        f.createDimension('lev_' + lev, nlevs)
        levvar = f.createVariable('lev_' + lev, 'f8', ('lev_' + lev), zlib=True, shuffle=False, complevel=9,
                                  fill_value=1e20)
        levvar[:] = levvals

        latbndvar = f.createVariable('latbnds_' + lev, 'f4', ('lev_' + lev, 'bnds'), zlib=True, shuffle=False,
                                     complevel=9, fill_value=1e20)
        latbndvar[:] = latbnds
        latbndvar.units = 'degrees_north'
        latbndvar.long_name = 'latitude bounds for ' + lev

        lonbndvar = f.createVariable('lonbnds_' + lev, 'f4', ('lev_' + lev, 'bnds'), zlib=True, shuffle=False,
                                     complevel=9, fill_value=1e20)
        lonbndvar[:] = lonbnds
        lonbndvar.units = 'degrees_east'
        lonbndvar.long_name = 'longitude bounds for ' + lev
