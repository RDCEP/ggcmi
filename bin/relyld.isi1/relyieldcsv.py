#!/usr/bin/env python

# import modules
from csv import writer
from re import findall
from os import listdir
from os.path import isdir
from fnmatch import filter
from netCDF4 import Dataset as nc
from optparse import OptionParser
from numpy import zeros, where, append

# parse inputs
parser = OptionParser()
parser.add_option("-i", "--inputdir", dest = "inputdir", default = ".", type = "string",
                  help = "Input directory")
parser.add_option("-w", "--weather", dest = "weather", default = "GFDL-ESM2M,HadGEM2-ES,IPSL-CM5A-LR,MIROC-ESM-CHEM,NorESM1-M", type = "string",
                  help = "Comma-separated list of GCMs")
parser.add_option("-g", "--cropmodels", dest = "cropmodels", default = "EPIC,GEPIC,IMAGE_LEITAP,LPJ-GUESS,LPJmL,pDSSAT,PEGASUS", type = "string",
                  help = "Comma-separated list of GGCMs")
parser.add_option("-c", "--crop", dest = "crop", default = "maize", type = "string",
                  help = "Crop (e.g., maize, rice, soy, wheat)")
parser.add_option("-r", "--rcps", dest = "rcps", default = "rcp2p6,rcp4p5,rcp6p0,rcp8p5", type = "string",
                  help = "Comma-separated list of RCPs")
parser.add_option("--co2", dest = "co2", default = "co2", type = "string",
                  help = "co2 setting (e.g., co2, noco2)")
parser.add_option("-d", "--decade", dest = "decade", default = "2010", type = "int",
                  help = "Start of decade")
parser.add_option("-o", "--outputdir", dest = "outputdir", default = ".", type = "string",
                  help = "Output directory", metavar = "FILE")
options, args = parser.parse_args()

inputdir   = options.inputdir
weather    = options.weather
cropmodels = options.cropmodels
crop       = options.crop
rcps       = options.rcps
co2        = options.co2
decade     = options.decade
outputdir  = options.outputdir

weather    = weather.split(',')
cropmodels = cropmodels.split(',')
rcps       = rcps.split(',')

files   = []
headers = []

cf2cs = {'maize': 'mai', 'rice': 'ric', 'soy': 'soy', 'wheat': 'whe'}

for w in range(len(weather)):
    for c in range(len(cropmodels)):
        for r in range(len(rcps)):
            drct = '%s/%s/%s/%s/%s/%s' % (inputdir, cropmodels[c], weather[w], crop, rcps[r], co2)
            if isdir(drct):
                fs = filter(listdir(drct), '*yield*')
                if len(fs) == 1:
                    files.append('%s/%s' % (drct, fs[0]))
                    headers.append('%s_%s_%s' % (weather[w], cropmodels[c], rcps[r]))

if not len(files):
    print 'No file matches. Exiting . . .'
    exit

with nc(files[0]) as f:
    gadm0 = f.variables['gadm0'][:]

data = zeros((len(gadm0), len(files)))

for i in range(len(files)):
    with nc(files[i]) as f:
        time  = f.variables['time'][:]
        time += int(findall(r'\d+', f.variables['time'].units)[0])

        try:
            htidx0, htidx1 = where(time == 1980)[0][0], where(time == 2009)[0][0] + 1
            ftidx0, ftidx1 = where(time == decade)[0][0], where(time == decade + 9)[0][0] + 1

            yld = f.variables['yield_%s_gadm0' % cf2cs[crop]][:]

            data[:, i] = yld[ftidx0 : ftidx1, :, 2].mean(axis = 0) / yld[htidx0 : htidx1, :, 2].mean(axis = 0) # sum
        except:
            print 'Warning with %s' % files[i]

outputfile = '%s/%s_%d_%s.csv' % (outputdir, crop, decade, co2)
with open(outputfile, 'w') as f:
    w = writer(f, delimiter = ',')
    w.writerow(['gadm0'] + headers)
    for i in range(len(data)):
        w.writerow(append(gadm0[i], data[i]))