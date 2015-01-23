#!/usr/bin/env python

# import modules
import matplotlib
from re import findall
from numpy import zeros, ones
from itertools import product
import matplotlib.pyplot as plt
from optparse import OptionParser
from netCDF4 import Dataset as nc
from matplotlib.patches import Polygon
from numpy.ma import masked_array, reshape, resize, arange, masked_where

# parse inputs
parser = OptionParser()
parser.add_option("-i", "--infile", dest = "infile", default = "metrics.nc4", type = "string",
                  help = "Metrics file", metavar = "FILE")
parser.add_option("-c", "--crop", dest = "crop", default = "maize", type = "string",
                  help = "Crop (or all)")
parser.add_option("-b", "--boxfile", dest = "boxfile", default = "box.png", type = "string",
                  help = "Output boxplot file", metavar = "FILE")
parser.add_option("-n", "--ncfile", dest = "ncfile", default = "box.nc4", type = "string",
                  help = "Output netcdf data file", metavar = "FILE")
options, args = parser.parse_args()

infile  = options.infile
crop    = options.crop
boxfile = options.boxfile
ncfile  = options.ncfile

cals = {'maize': 3.60, 'wheat': 3.34, 'soy': 3.35, 'rice': 3.59}

with nc(infile) as f:
    decade = f.variables['decade'][:]
    dunits = f.variables['decade'].units
    decade = decade * 10 + int(findall(r'\d+', dunits)[0])

models = ['epic', 'gepic', 'image_leitap', 'lpj-guess', 'lpjml', 'pdssat', 'pegasus']
gcms   = ['gfdl-esm2m', 'hadgem2-es', 'ipsl-cm5a-lr', 'miroc-esm-chem', 'noresm1-m']
crops  = ['maize', 'wheat', 'soy', 'rice'] if crop == 'all' else [crop]
co2s   = ['co2', 'noco2']

nm, ng, ncr, nd, nco2 = len(models), len(gcms), len(crops), len(decade), len(co2s)

# benefit
sh = (nm, ng, ncr, nd, nco2)
benefit = masked_array(zeros(sh), mask = ones(sh))
with nc(infile) as f:
    for m, g, c, co in product(range(nm), range(ng), range(ncr), range(nco2)):
        var = 'benefit_global_%s_%s_%s_%s' % (models[m], gcms[g], crops[c], co2s[co])
        if var in f.variables:
            benefit[m, g, c, :, co] = f.variables[var][:, 0, 2] # all decades, sum

# weights
weights = masked_array(zeros(sh), mask = ones(sh))
for i in range(ncr):
    weights[:, :, i] = resize(cals[crops[i]], (nm, ng, nd, 2))
weights = masked_where(benefit.mask, weights) # mask

# average
wbenefit  = (benefit * weights).sum(axis = 2)
wbenefit /= weights.sum(axis = 2)

hadgemidx = gcms.index('hadgem2-es')

bps = [0] * 3

plt.figure(figsize = (10, 5))
ax = plt.subplot(111)

# hadgem noco2
p1 = wbenefit[[0, 1, 3, 4, 5, 6], hadgemidx, :, 1]
bp1 = ax.boxplot(p1, positions = range(1, 4 * nd, 4))
bps[0] = bp1

# hadgem co2
p2 = wbenefit[[0, 1, 3, 4, 5, 6], hadgemidx, :, 0]
bp2 = ax.boxplot(p2, positions = range(2, 4 * nd, 4))
bps[1] = bp2

# all co2
p3 = reshape(wbenefit[:, :, :, 0], (nm * ng, nd))
bp3 = ax.boxplot(p3, positions = range(3, 4 * nd, 4))
bps[2] = bp3

colors = ['r', 'y', 'b']

# change colors
for i in range(len(bps)):
    for j in range(len(bps[i]['boxes'])):
        box = bps[i]['boxes'][j]
        boxx, boxy = box.get_xdata(), box.get_ydata()
        boxcoords = zip(boxx, boxy)
        boxpolygon = Polygon(boxcoords, facecolor = colors[i])
        ax.add_patch(boxpolygon)

        med = bps[i]['medians'][j]
        medianx, mediany = med.get_xdata(), med.get_ydata()
        plt.plot(medianx, mediany, 'k')

plt.xlim([0, 4 * nd])
plt.xticks(arange(2, 4 * nd, 4), ['%d-%d' % (d, d + 9) for d in decade], rotation = 45)
plt.ylim([-0.5, 0.5])
plt.yticks(arange(-0.5, 0.6, 0.25))
plt.grid(which = 'major', axis = 'y')
plt.tight_layout()

l1, = plt.plot([1, 1], 'r')
l2, = plt.plot([1, 1], 'y')
l3, = plt.plot([1, 1], 'b')
plt.legend((l1, l2, l3), ['no co2 (6)', 'co2 HadGEM2-ES (6)', 'co2 (35)'], loc = 'lower left')
l1.set_visible(False)
l2.set_visible(False)
l3.set_visible(False)

# save
plt.savefig(boxfile)
plt.close()

# write data file
with nc(ncfile, 'w') as f:
    f.createDimension('models', nm)
    mvar = f.createVariable('models', 'i4', 'models')
    mvar[:] = range(1, 1 + nm)
    mvar.units = 'mapping'
    mvar.long_name = ', '.join(models)

    f.createDimension('gcms', ng)
    gvar = f.createVariable('gcms', 'i4', 'gcms')
    gvar[:] = range(1, 1 + ng)
    gvar.units = 'mapping'
    gvar.long_name = ', '.join(gcms)

    f.createDimension('decade', nd)
    dvar = f.createVariable('decade', 'i4', 'decade')
    dvar[:] = range(nd)
    dvar.units = 'decades since from 1980'
    dvar.long_name = 'decade'

    f.createDimension('co2', nco2)
    cvar = f.createVariable('co2', 'i4', 'co2')
    cvar[:] = range(1, 1 + nco2)
    cvar.units = 'mapping'
    cvar.long_name = ', '.join(co2s)

    bvar = f.createVariable('benefit', 'f8', ('models', 'gcms', 'decade', 'co2'), zlib = True, shuffle = False, complevel = 9, fill_value = 1e20)
    bvar[:] = wbenefit
    bvar.long_name = 'benefit of RCP 2.6 over RCP 8.5 at global level'