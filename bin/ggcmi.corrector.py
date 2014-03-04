#!/usr/bin/env python

# import modules
import numpy.ma as ma
from numpy import array
from shutil import copyfile
from optparse import OptionParser
from netCDF4 import Dataset as nc

def loadData(file):
    with nc(file) as f:
        var = array(f.variables.keys())
        idx = array([v != 'lon' and v != 'lat' and v != 'time' for v in var])
        varname = var[idx][0]
        data = f.variables[varname][:]
    return data, varname

parser = OptionParser()
parser.add_option("-i", "--inputfile", dest = "inputfile", default = 64, type = "string",
                  help = "Input file", metavar = "FILE")
parser.add_option("-p", "--pfile", dest = "plantingfile", default = 64, type = "string",
                  help = "Planting date file", metavar = "FILE")       
parser.add_option("-m", "--mfile", dest = "maturityfile", default = 64, type = "string",
                  help = "Maturity date file", metavar = "FILE")
parser.add_option("-o", "--outputfile", dest = "outputfile", default = 64, type = "string",
                  help = "Output file", metavar = "FILE")             
(options, args) = parser.parse_args()

inputfile    = options.inputfile # get options
plantingfile = options.plantingfile
maturityfile = options.maturityfile
outputfile   = options.outputfile

inputd, varname = loadData(inputfile) # load data
pdates = loadData(plantingfile)[0]
mdates = loadData(maturityfile)[0]

pdates = pdates[0] # select first year
mdates = mdates[0]

latidx, lonidx = ma.where(365 - pdates < mdates)
shiftd = inputd[:, latidx, lonidx]
shiftd[: -1] = shiftd[1 :] # shift data
shiftd[-1].mask = True
inputd[:, latidx, lonidx] = shiftd

copyfile(inputfile, outputfile)
with nc(outputfile, 'r+') as f: # modify existing variable
    v = f.variables[varname]
    v[:] = inputd