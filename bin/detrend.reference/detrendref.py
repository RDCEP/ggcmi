#!/usr/bin/env python

# add paths
import os, sys
for p in os.environ['PATH'].split(':'): sys.path.append(p)

# import modules
from numpy import zeros, ones
from itertools import product
from numpy.ma import masked_array
from optparse import OptionParser
from netCDF4 import Dataset as nc
from detrender import DetrenderWrapper

# parse inputs
parser = OptionParser()
parser.add_option("-i", "--inputfile", dest = "inputfile", default = "ray.1961-2008.fpu.nc4", type = "string",
                  help = "Input file", metavar = "FILE")
parser.add_option("-a", "--agglvl", dest = "agglvl", default = "fpu", type = "string",
                  help = "Aggregation level")
parser.add_option("-c", "--crops", dest = "crops", default = "mai,ric,soy,whe", type = "string",
                  help = "Comma-separated list of crops")
parser.add_option("-o", "--outputfile", dest = "outputfile", default = "out.nc4", type = "string",
                  help = "Output file", metavar = "FILE")
options, args = parser.parse_args()

inputfile  = options.inputfile
agglvl     = options.agglvl
crops      = options.crops.split(',')
outputfile = options.outputfile 

dt = ['none', 'lin', 'quad', 'ma', 'ffd', 'ffdtr']
mp = ['true', 'false']
ndt, nmp = len(dt), len(mp)

fi = nc(inputfile,  'r')
fo = nc(outputfile, 'w')

aggs  = fi.variables[agglvl][:]
naggs = len(aggs)

time  = fi.variables['time'][:]
ntime = len(time)

sh = (naggs, ntime, ndt, nmp)

fo.createDimension(agglvl, naggs)
aggvar = fo.createVariable(agglvl, 'i4', agglvl, zlib = True, shuffle = False, complevel = 9)
aggvar[:] = aggs
aggvar.units = fi.variables[agglvl].units
aggvar.long_name = fi.variables[agglvl].long_name

fo.createDimension('time', ntime)
timevar = fo.createVariable('time', 'i4', 'time', zlib = True, shuffle = False, complevel = 9)
timevar[:] = time
timevar.units = fi.variables['time'].units
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

for cr in crops:
    yld = fi.variables['yield_' + cr][:]

    ylddt = masked_array(zeros(sh), mask = ones(sh))
    for d, c in product(range(ndt), range(nmp)):
        detrender = DetrenderWrapper(dt[d], mp[c])
        for i in range(naggs):
            ylddt[i, :, d, c] = detrender.detrend(yld[i, :])[0]

    yvar = fo.createVariable('yield_' + cr, 'f4', (agglvl, 'time', 'dt', 'mp'), zlib = True, shuffle = False, complevel = 9, fill_value = 1e20)
    yvar[:] = ylddt
    yvar.units = fi.variables['yield_' + cr].units
    yvar.long_name = fi.variables['yield_' + cr].long_name

fi.close()
fo.close()