#!/usr/bin/env python

# import modules
from optparse import OptionParser
from netCDF4 import Dataset as nc
from numpy import zeros, isnan, cos, pi
from numpy.ma import resize, masked_array, where, logical_or, logical_and, logical_not, masked_where

# parse inputs
parser = OptionParser()
parser.add_option("-i", "--inputfile", dest = "inputfile", default = "mai_weight_ray_1961-2008.nc4", type = "string",
                  help = "Input file", metavar = "FILE")
parser.add_option("-w", "--weightfile", dest = "weightfile", default = "maize.nc4", type = "string",
                  help = "Weight mask file", metavar = "FILE")
parser.add_option("-v", "--variable", dest = "variable", default = "area_mai", type = "string",
                  help = "Area variable")
parser.add_option("-y", "--year", dest = "year", default = 1961, type = "int",
                  help = "Start year")
parser.add_option("-o", "--outputfile", dest = "outputfile", default = "out.nc4", type = "string",
                  help = "Output file", metavar = "FILE")
options, args = parser.parse_args()

inputfile  = options.inputfile
weightfile = options.weightfile
variable   = options.variable
year       = options.year
outputfile = options.outputfile

with nc(inputfile) as f:
    lats, lons = f.variables['lat'][:], f.variables['lon'][:]
    areasum1 = f.variables[variable][:]
    aunits = f.variables[variable].units

areasum1 = masked_where(isnan(areasum1), areasum1)

mk = areasum1.mask
sh = areasum1.shape

if aunits == '%': # convert percent to area
    area = 100 * (111.2 / 2) ** 2 * cos(pi * lats / 180)
    area = resize(area, (len(lons), len(lats))).T
    area = resize(area, sh) # resize
    area = masked_where(mk, area) # mask
    areasum1 = areasum1 * area / 100.

with nc(weightfile) as f:
    areair2 = f.variables['irrigated'][:]
    arearf2 = f.variables['rainfed'][:]

areair2  = resize(areair2, sh)
arearf2  = resize(arearf2, sh)
areasum2 = areair2 + arearf2

areair1 = masked_array(zeros(sh), mask = mk)
arearf1 = masked_array(zeros(sh), mask = mk)

# no data -> all rainfed
nodata = logical_and(~mk, logical_or(areasum2.mask, areasum2 == 0))
idx1, idx2, idx3 = where(nodata)
areair1[idx1, idx2, idx3] = 0.
arearf1[idx1, idx2, idx3] = areasum1[idx1, idx2, idx3]

# has data
hasdata = logical_and(~mk, logical_not(nodata))
idx1, idx2, idx3 = where(hasdata)
areair1[idx1, idx2, idx3] = areasum1[idx1, idx2, idx3] * areair2[idx1, idx2, idx3] / areasum2[idx1, idx2, idx3]
arearf1[idx1, idx2, idx3] = areasum1[idx1, idx2, idx3] * arearf2[idx1, idx2, idx3] / areasum2[idx1, idx2, idx3]

with nc(outputfile, 'w') as f:
    f.createDimension('time', sh[0])
    tvar = f.createVariable('time', 'i4', 'time')
    tvar[:] = range(0, sh[0])
    tvar.units = 'years since %d' % year
    tvar.long_name = 'time'

    f.createDimension('lat', sh[1])
    latvar = f.createVariable('lat', 'f8', 'lat')
    latvar[:] = lats
    latvar.units = 'degrees_north'
    latvar.long_name = 'latitude'

    f.createDimension('lon', sh[2])
    lonvar = f.createVariable('lon', 'f8', 'lon')
    lonvar[:] = lons
    lonvar.units = 'degrees_east'
    lonvar.long_name = 'longitude'

    ivar = f.createVariable('irrigated', 'f8', ('time', 'lat', 'lon'), zlib = True, shuffle = False, complevel = 9, fill_value = 1e20)
    ivar[:] = areair1
    ivar.units = 'hectares'
    ivar.long_name = 'irrigated harvested area'

    rvar = f.createVariable('rainfed', 'f8', ('time', 'lat', 'lon'), zlib = True, shuffle = False, complevel = 9, fill_value = 1e20)
    rvar[:] = arearf1
    rvar.units = 'hectares'
    rvar.long_name = 'rainfed harvested area'

    svar = f.createVariable('sum', 'f8', ('time', 'lat', 'lon'), zlib = True, shuffle = False, complevel = 9, fill_value = 1e20)
    svar[:] = areasum1
    svar.units = 'hectares'
    svar.long_name = 'total harvested area'