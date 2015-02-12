#!/usr/bin/env python

# add paths
import os, sys
for p in os.environ['PATH'].split(':'): sys.path.append(p)

# import modules
import fnmatch
from re import findall
from os import listdir, sep
from os.path import basename
from numpy.ma import masked_array
from netCDF4 import Dataset as nc
from optparse import OptionParser
from filespecs import RescaledFile
from numpy import where, unique, array, ones, zeros, logical_and

parser = OptionParser()
parser.add_option("-i", "--infile", dest = "infile", default = "tscorr_agmerra_hist_mai_annual_1980_2010.ensemble.nc4", type = "string",
                  help = "Input model ensemble file", metavar = "FILE")
parser.add_option("-d", "--indir", dest = "indir", default = "", type = "string",
                  help = "Directory of model rescaled files")
parser.add_option("-m", "--mkfile", dest = "mkfile", default = "", type = "string",
                  help = "Mask file", metavar = "FILE")
parser.add_option("-c", "--crmthd", dest = "crmthd", default = "2,0,2", type = "string",
                  help = "Comma-separated list of dt,mp,cr indices")
parser.add_option("-a", "--agglvl", dest = "agglvl", default = "gadm0", type = "string",
                  help = "Aggregation level (e.g., gadm0, fpu, kg)")
parser.add_option("-o", "--outfile", dest = "outfile", default = "", type = "string",
                  help = "Output ensemble rescaled spatial file")
options, args = parser.parse_args()

infile  = options.infile
indir   = options.indir
mkfile  = options.mkfile
crmthd  = options.crmthd
agglvl  = options.agglvl
outfile = options.outfile

dtidx, mpidx, cridx = [int(m) for m in crmthd.split(',')]

climate, crop = [basename(infile).split('_')[idx] for idx in [1, 3]]

vname = 'yield_' + crop

with nc(mkfile) as f:
    aggmap = f.variables[agglvl][:]

with nc(infile) as f:
    aggs = f.variables[agglvl][:]

    time  = f.variables['time'][:]
    time += int(findall(r'\d+', f.variables['time'].units)[0])

    modelnames = f.variables['model_order'].long_name.split(', ')
    modelorder = f.variables['model_order'][:, dtidx, mpidx, cridx, 0] # top model

    scennames = f.variables['top_scens'].long_name.split(', ')
    topscens  = f.variables['top_scens'][:, dtidx, mpidx, cridx, 0] # top scenario

mask = modelorder.mask # remove masked values
aggs = aggs[~mask]
modelorder = modelorder[~mask]
topscens = topscens[~mask]

naggs = len(aggs)

modelscen = [0] * naggs
for i in range(naggs):
    modelscen[i] = modelnames[int(modelorder[i] - 1)] + ',' + scennames[int(topscens[i] - 1)]

umodelscen = unique(modelscen) # unique model-scenarios

nuniq = len(umodelscen)

uaggs = [0] * nuniq
for i in range(nuniq):
    idx = where(array(modelscen) == umodelscen[i])[0]
    uaggs[i] = aggs[idx]

dirlist = listdir(indir)

files = [0] * nuniq
for i in range(nuniq):
    model, scen = umodelscen[i].split(',')

    f = fnmatch.filter(dirlist, '%s_%s*%s_yield*%s*' % (model, climate, scen, crop))
    if len(f) != 1:
        raise Exception('Could not find single file match for %s, %s, %s, %s' % (model, climate, scen, crop))

    files[i] = indir + sep + f[0]

with nc(files[0]) as f: # read first file for dimensions, etc.
    lats, lons = f.variables['lat'][:], f.variables['lon'][:]

    irr = f.variables['irr'].long_name.split(', ')

    vunits = f.variables[vname].units
    vlname = f.variables[vname].long_name

sh = (len(time), len(lats), len(lons), len(irr))
yld = masked_array(zeros(sh), mask = ones(sh))

for i in range(nuniq):
    with nc(files[i]) as f:
        timef  = f.variables['time'][:]
        timef += int(findall(r'\d+', f.variables['time'].units)[0])

        yldf = f.variables[vname][:] # yield(time, lat, lon, irr)

    tmin = max(time.min(), timef.min())
    tmax = min(time.max(), timef.max())

    tidx  = logical_and(time  >= tmin, time  <= tmax)
    tidxf = logical_and(timef >= tmin, timef <= tmax)

    for j in range(len(uaggs[i])):
        latidx, lonidx = where(aggmap == uaggs[i][j])

        ytmp       = yld[:,  latidx, lonidx, :].copy()
        ytmp[tidx] = yldf[:, latidx, lonidx, :][tidxf]

        yld[:, latidx, lonidx, :] = ytmp

fout = RescaledFile(outfile, time, lats, lons, irr)
fout.append(vname, yld, ('time', 'lat', 'lon', 'irr'), vunits, vlname)