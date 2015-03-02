#!/usr/bin/env python

# add paths
import os, sys
for p in os.environ['PATH'].split(':'): sys.path.append(p)

# import modules
from re import findall
from netCDF4 import Dataset as nc
from optparse import OptionParser
from filespecs import RescaledFile
from os.path import basename, splitext, isfile
from numpy.ma import masked_array, masked_where
from numpy import ones, zeros, where, isnan, cos, pi, resize, logical_and

parser = OptionParser()
parser.add_option("-i", "--irfile", dest = "irfile", default = "", type = "string",
                  help = "Input irrigated spatial file", metavar = "FILE")
parser.add_option("-r", "--rffile", dest = "rffile", default = "", type = "string",
                  help = "Input rainfed spatial file", metavar = "FILE")
parser.add_option("-b", "--bcfile", dest = "bcfile", default = "", type = "string",
                  help = "Bias corrected file", metavar = "FILE")
parser.add_option("-m", "--mkfile", dest = "mkfile", default = "", type = "string",
                  help = "Mask file", metavar = "FILE")
parser.add_option("-w", "--wtfile", dest = "wtfile", default = "", type = "string",
                  help = "Weights file", metavar = "FILE")
parser.add_option("-c", "--crmthd", dest = "crmthd", default = "2,0,2", type = "string",
                  help = "Comma-separated list of dt,mp,cr indices")
parser.add_option("-a", "--agglvl", dest = "agglvl", default = "gadm0", type = "string",
                  help = "Aggregation level (e.g., gadm0, fpu, kg)")
parser.add_option("-v", "--vname", dest = "vname", default = "", type = "string",
                  help = "Variable name")
parser.add_option("-o", "--outfile", dest = "outfile", default = "", type = "string",
                  help = "Output rescaled spatial file")
options, args = parser.parse_args()

irfile  = options.irfile
rffile  = options.rffile
bcfile  = options.bcfile
mkfile  = options.mkfile
wtfile  = options.wtfile
crmthd  = options.crmthd
agglvl  = options.agglvl
vname   = options.vname
outfile = options.outfile

dtidx, mpidx, cridx = [int(m) for m in crmthd.split(',')]

fsplit = basename(irfile).split('_')
scen   = fsplit[3] if len(fsplit) == 10 else '_'.join([fsplit[i] for i in [3, 5]])

with nc(wtfile) as f: # load weights
    lat, lon = f.variables['lat'][:], f.variables['lon'][:]
    wir, wrf = f.variables['irrigated'][:], f.variables['rainfed'][:]

with nc(mkfile) as f: # load mask
    aggmap = f.variables[agglvl][:]

with nc(bcfile) as f: # load bias-corrected file
    aggs   = f.variables[agglvl][:]
    time   = f.variables['time'][:]
    tunits = f.variables['time'].units

    snames = f.variables['scen'].long_name.split(', ')
    if not scen in snames:
        print 'Cannot find scenario %s. Exiting . . .' % scen
        sys.exit()
    else:
        sidx = snames.index(scen)

    ydtr = f.variables['yield_detrend'][:, :, sidx, 0,     0,     0]
    yrtr = f.variables['yield_retrend'][:, :, sidx, dtidx, mpidx, cridx]

time += int(findall(r'\d+', tunits)[0]) # get time in years

ydtr = masked_where(isnan(ydtr), ydtr) # convert NaNs to masked
yrtr = masked_where(isnan(yrtr), yrtr)

(nlats, nlons), nt, nirr = aggmap.shape, len(time), 3 # get dimensions

area    = 100 * (111.2 / 2) ** 2 * cos(pi * lat / 180) # get areas
area    = resize(area, (nlons, nlats)).T
areair  = resize(area * wir, (nt, nlats, nlons))
arearf  = resize(area * wrf, (nt, nlats, nlons))
areatot = areair + arearf
areatot = masked_where(areatot == 0, areatot)

varr = masked_array(zeros((nt, nlats, nlons, nirr)), mask = ones((nt, nlats, nlons, nirr)))

y1, y2 = [int(y) for y in findall(r'\d+', splitext(basename(irfile))[0])[-2 :]]
ftime  = range(y1, y2 + 1)

with nc(irfile) as f:
    var = f.variables[vname]

    units = var.units     if 'units'     in var.ncattrs() else ''
    lname = var.long_name if 'long_name' in var.ncattrs() else ''

    varr[:, :, :, 0] = var[logical_and(ftime >= time[0], ftime <= time[-1])]

if isfile(rffile):
    with nc(rffile) as f:
        varr[:, :, :, 1] = f.variables[vname][logical_and(ftime >= time[0], ftime <= time[-1])]

varr[isnan(varr)] = 0.

for i in range(len(aggs)):
    latidx, lonidx = where(aggmap == aggs[i])
    for j in range(len(latidx)):
        varr[:, latidx[j], lonidx[j], 0] *= yrtr[i, :] / ydtr[i, :]
        varr[:, latidx[j], lonidx[j], 1] *= yrtr[i, :] / ydtr[i, :]

varr[:, :, :, 2] = (areair * varr[:, :, :, 0] + arearf * varr[:, :, :, 1]) / areatot

fout = RescaledFile(outfile, time, lat, lon, ['ir', 'rf', 'sum'])
fout.append(vname, varr, ('time', 'lat', 'lon', 'irr'), units, lname)
