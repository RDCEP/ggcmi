#!/usr/bin/env python

# add paths
import os, sys
for p in os.environ['PATH'].split(':'): sys.path.append(p)

# import modules
from numpy import zeros, ones
from itertools import product
from optparse import OptionParser
from netCDF4 import Dataset as nc
from detrender import DetrenderWrapper
from numpy.ma import masked_array, isMaskedArray

# parse inputs
parser = OptionParser()
parser.add_option("-i", "--inputfile", dest = "inputfile", default = "mai_weight_ray_1961-2008.nc4", type = "string",
                  help = "Input file", metavar = "FILE")
parser.add_option("-v", "--variables", dest = "variables", default = "yield_mai", type = "string",
                  help = "Comma-separated list of variables")
parser.add_option("-o", "--outputfile", dest = "outputfile", default = "out.nc4", type = "string",
                  help = "Output file", metavar = "FILE")
options, args = parser.parse_args()

inputfile  = options.inputfile
variables  = options.variables
outputfile = options.outputfile

dt = ['none', 'lin', 'quad', 'ma', 'ffd', 'ffdtr']
mp = ['true', 'false']
ndt, nmp = len(dt), len(mp)

variables = variables.split(',')

fi = nc(inputfile,  'r')
fo = nc(outputfile, 'w')

time  = fi.variables['time'][:]
ntime = len(time)

if 'lat' in fi.variables and 'lon' in fi.variables:
    latname = 'lat'
    lonname = 'lon'
elif 'latitude' in fi.variables and 'longitude' in fi.variables:
    latname = 'latitude'
    lonname = 'longitude'
else:
    raise Exception('Uncannot find latitude, longitude variables')

lats, lons   = fi.variables[latname][:], fi.variables[lonname][:]
nlats, nlons = len(lats), len(lons) 

sh = (ntime, nlats, nlons, ndt, nmp)

fo.createDimension(latname, nlats)
latvar = fo.createVariable(latname, 'f8', latname, zlib = True, shuffle = False, complevel = 9)
latvar[:] = lats
latvar.units = 'degrees_north'
latvar.long_name = 'latitude'

fo.createDimension(lonname, nlons)
lonvar = fo.createVariable(lonname, 'f8', lonname, zlib = True, shuffle = False, complevel = 9)
lonvar[:] = lons
lonvar.units = 'degrees_east'
lonvar.long_name = 'longitude'

fo.createDimension('time', ntime)
timevar = fo.createVariable('time', 'i4', 'time', zlib = True, shuffle = False, complevel = 9)
timevar[:] = time
timevar.units = fi.variables['time'].units if 'units' in fi.variables['time'].ncattrs() else ''
timevar.long_name = 'time'

fo.createDimension('dt', ndt)
dtvar = fo.createVariable('dt', 'i4', 'dt', zlib = True, shuffle = False, complevel = 9)
dtvar[:] = range(1, ndt + 1)
dtvar.units = 'mapping'
dtvar.long_name = ', '.join(dt)

fo.createDimension('mp', nmp)
mpvar = fo.createVariable('mp', 'i4', 'mp', zlib = True, shuffle = False, complevel = 9)
mpvar[:] = range(1, nmp + 1)
mpvar.units = 'mapping'
mpvar.long_name = ', '.join(mp)

for v in variables:
    var = fi.variables[v]

    vardt = masked_array(zeros(sh), mask = ones(sh))
    for latidx, lonidx in product(range(nlats), range(nlons)):
        varmat = var[:, latidx, lonidx]
        if not isMaskedArray(varmat) or not varmat.mask[0]:
            for d, c in product(range(ndt), range(nmp)):
                detrender = DetrenderWrapper(dt[d], mp[c])
                vardt[:, latidx, lonidx, d, c] = detrender.detrend(varmat)[0]

    vvar = fo.createVariable(v, 'f4', ('time', latname, lonname, 'dt', 'mp'), zlib = True, shuffle = False, complevel = 9, fill_value = 1e20)
    vvar[:] = vardt
    vvar.units = var.units if 'units' in var.ncattrs() else ''
    vvar.long_name = var.long_name if 'long_name' in var.ncattrs() else ''

fi.close()
fo.close()