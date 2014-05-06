#!/usr/bin/env python

# import modules
import numpy.ma as ma
from numpy import array
from shutil import copyfile
from optparse import OptionParser
from netCDF4 import Dataset as nc

parser = OptionParser()
parser.add_option("-i", "--inputfile", dest = "inputfile", default = 64, type = "string",
                  help = "Input file", metavar = "FILE")
parser.add_option("-g", "--gfile", dest = "growingfile", default = 64, type = "string",
                  help = "Growing season dates file", metavar = "FILE")
parser.add_option("-o", "--outputfile", dest = "outputfile", default = 64, type = "string",
                  help = "Output file", metavar = "FILE")
(options, args) = parser.parse_args()

inputfile   = options.inputfile # get options
growingfile = options.growingfile
outputfile  = options.outputfile

with nc(inputfile) as f: # load data
    var = array(f.variables.keys())
    idx = array([v != 'lon' and v != 'lat' and v != 'time' for v in var])
    varname = var[idx][0]
    inputd = f.variables[varname][:]
with nc(growingfile) as f:
    lats = f.variables['lat'][:]
    pdates = f.variables['planting day'][:]
    mdates = f.variables['growing season length'][:]
si = [i[0] for i in sorted(enumerate(lats), reverse = True, key = lambda x : x[1])]
pdates = pdates[:, si, :] # flip latitudes
mdates = mdates[:, si, :]

pdates = ma.masked_where(pdates < 0, pdates) # mask < 0
mdates = ma.masked_where(mdates < 0, mdates)

latidx, lonidx = ma.where(365 - pdates < mdates)
shiftd = inputd[:, latidx, lonidx]
shiftd[: -1] = shiftd[1 :] # shift data
shiftd[-1].mask = True
inputd[:, latidx, lonidx] = shiftd

copyfile(inputfile, outputfile)
with nc(outputfile, 'r+') as f: # modify existing variable
    v = f.variables[varname]
    v[:] = inputd
