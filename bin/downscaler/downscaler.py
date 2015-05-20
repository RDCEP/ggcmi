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
from numpy import ones, zeros, where, isnan, cos, pi, resize, logical_and, intersect1d

parser = OptionParser()
parser.add_option("-i", "--irfile", dest = "irfile", default = "", type = "string",
                  help = "Input irrigated spatial file", metavar = "FILE")
parser.add_option("-r", "--rffile", dest = "rffile", default = "", type = "string",
                  help = "Input rainfed spatial file", metavar = "FILE")
parser.add_option("-a", "--agfile", dest = "agfile", default = "", type = "string",
                  help = "Aggregated file", metavar = "FILE")
parser.add_option("-f", "--reffile", dest = "reffile", default = "", type = "string",
                  help = "Reference file", metavar = "FILE")
parser.add_option("-m", "--mkfile", dest = "mkfile", default = "", type = "string",
                  help = "Mask file", metavar = "FILE")
parser.add_option("-w", "--wtfile", dest = "wtfile", default = "", type = "string",
                  help = "Weights file", metavar = "FILE")
parser.add_option("--agglvl", dest = "agglvl", default = "gadm0", type = "string",
                  help = "Aggregation level (e.g., gadm0, fpu, kg)")
parser.add_option("-v", "--vname", dest = "vname", default = "", type = "string",
                  help = "Variable name")
parser.add_option("-o", "--outfile", dest = "outfile", default = "", type = "string",
                  help = "Output downscaled spatial file")
options, args = parser.parse_args()

irfile  = options.irfile
rffile  = options.rffile
agfile  = options.agfile
reffile = options.reffile
mkfile  = options.mkfile
wtfile  = options.wtfile
agglvl  = options.agglvl
vname   = options.vname
outfile = options.outfile

fsplit = basename(irfile).split('_')
scen   = fsplit[3] if len(fsplit) == 10 else '_'.join([fsplit[i] for i in [3, 5]])
crop   = fsplit[6] if len(fsplit) == 10 else fsplit[7]

with nc(wtfile) as f: # load weights
    lat, lon = f.variables['lat'][:], f.variables['lon'][:]
    wir, wrf = f.variables['irrigated'][:], f.variables['rainfed'][:]

with nc(mkfile) as f: # load mask
    aggmap = f.variables[agglvl][:]

with nc(agfile) as f: # load aggregated file
    ain       = f.variables[agglvl][:]
    tin       = f.variables['time'][:]
    tin_units = f.variables['time'].units

    snames = f.variables['scen'].long_name.split(', ')
    if not scen in snames:
        print 'Cannot find scenario %s. Exiting . . .' % scen
        sys.exit()
    else:
        sidx = snames.index(scen)

    sum_idx  = f.variables['irr'].long_name.split(', ').index('sum')
    yield_in = f.variables['yield_' + agglvl][:, :, sidx, sum_idx]

with nc(reffile) as f: # load reference file
    aref       = f.variables[agglvl][:]
    tref       = f.variables['time'][:]
    tref_units = f.variables['time'].units
    dtidx      = f.variables['dt'].long_name.split(', ').index('none')
    mpidx      = f.variables['mp'].long_name.split(', ').index('true')

    var = 'yield_' + crop
    if var in f.variables:
        yield_ref = f.variables[var][:, :, dtidx, mpidx]
    else:
        print 'Crop %s unavailable in reference file %s. Exiting . . .' % (crop, reffile)
        sys.exit()

tref += int(findall(r'\d+', tref_units)[0])    # get reference time
tin  += int(findall(r'\d+', tin_units)[0]) - 1 # get simulation time

yield_ref = masked_where(isnan(yield_ref), yield_ref) # convert NaNs to masked
yield_in  = masked_where(isnan(yield_in),  yield_in)

time = intersect1d(tref, tin) # common times
aggs = intersect1d(aref, ain) # common aggregations

yield_ref = yield_ref[:, logical_and(tref >= time[0], tref <= time[-1])]
yield_in  = yield_in[:,  logical_and(tin  >= time[0], tin  <= time[-1])]

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

    if len(var) != len(ftime): # happens with epic-test
        ftime = f.variables['time'][:] + int(findall(r'\d+', f.variables['time'].units)[0]) - 1

    varr[:, :, :, 0] = var[logical_and(ftime >= time[0], ftime <= time[-1])]

if isfile(rffile):
    with nc(rffile) as f:
        varr[:, :, :, 1] = f.variables[vname][logical_and(ftime >= time[0], ftime <= time[-1])]

varr[isnan(varr)] = 0.

for i in range(len(aggs)):
    latidx, lonidx = where(aggmap == aggs[i])
    aidx1 = where(ain  == aggs[i])[0][0]
    aidx2 = where(aref == aggs[i])[0][0]
    for j in range(len(latidx)):
        varr[:, latidx[j], lonidx[j], 0] *= yield_ref[aidx2, :] / yield_in[aidx1, :]
        varr[:, latidx[j], lonidx[j], 1] *= yield_ref[aidx2, :] / yield_in[aidx1, :]

varr[:, :, :, 2] = (areair * varr[:, :, :, 0] + arearf * varr[:, :, :, 1]) / areatot

fout = RescaledFile(outfile, time, lat, lon, ['ir', 'rf', 'sum'])
fout.append(vname, varr, ('time', 'lat', 'lon', 'irr'), units, lname)
