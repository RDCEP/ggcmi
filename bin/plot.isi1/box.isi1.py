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
from numpy.ma import masked_array, reshape, arange, masked_where, median

# parse inputs
parser = OptionParser()
parser.add_option("-i", "--infile", dest = "infile", default = "metrics.nc4", type = "string",
                  help = "Metrics file", metavar = "FILE")
parser.add_option("-c", "--crop", dest = "crop", default = "maize", type = "string",
                  help = "Crop (or all)")
parser.add_option("-r", "--hareafile", dest = "hareafile", default = "all.global.nc4", type = "string",
                  help = "Harvested area file", metavar = "FILE")
parser.add_option("-v", "--variable", dest = "variable", default = "beta", type = "string",
                  help = "Variable to plot")
parser.add_option("-b", "--boxfile", dest = "boxfile", default = "box.png", type = "string",
                  help = "Output boxplot file", metavar = "FILE")
parser.add_option("-n", "--ncfile", dest = "ncfile", default = "box.nc4", type = "string",
                  help = "Output netcdf data file", metavar = "FILE")
options, args = parser.parse_args()

infile    = options.infile
crop      = options.crop
hareafile = options.hareafile
variable  = options.variable
boxfile   = options.boxfile
ncfile    = options.ncfile

cals = {'maize': 3.60, 'wheat': 3.34, 'soy': 3.35, 'rice': 2.80}

with nc(infile) as f:
    decade = f.variables['decade'][:]
    dunits = f.variables['decade'].units
    decade = decade * 10 + int(findall(r'\d+', dunits)[0])

careas = {}
with nc(hareafile) as f:
    for c in ['maize', 'wheat', 'soy', 'rice']:
        careas[c] = f.variables['area_' + c][:]

models = ['epic', 'gepic', 'lpj-guess', 'lpjml', 'pdssat', 'pegasus'] # exclude image
gcms   = ['gfdl-esm2m', 'hadgem2-es', 'ipsl-cm5a-lr', 'miroc-esm-chem', 'noresm1-m']
crops  = ['maize', 'wheat', 'soy', 'rice'] if crop == 'all' else [crop]
co2s   = ['co2', 'noco2']

hadgemidx  = gcms.index('hadgem2-es')

nm, ng, ncr, nd, nco2 = len(models), len(gcms), len(crops), len(decade), len(co2s)

# variable
sh = (nm, ng, ncr, nd, nco2)
varr = masked_array(zeros(sh), mask = ones(sh))
with nc(infile) as f:
    for m, g, c, co in product(range(nm), range(ng), range(ncr), range(nco2)):
        var = '%s_global_%s_%s_%s_%s' % (variable, models[m], gcms[g], crops[c], co2s[co])
        if var in f.variables:
            varr[m, g, c, :, co] = f.variables[var][0, :] # all decades

if crop == 'all':
    # replace rice for pegasus with median over other crops
    pegasusidx = models.index('pegasus')
    riceidx    = crops.index('rice')
    notriceidx = range(len(crops))
    notriceidx.pop(riceidx);

    for g, co in product(range(ng), range(nco2)):
        varr[pegasusidx, g, riceidx, :, co] = median(varr[pegasusidx, g, notriceidx, :, co], axis = 0)

# weights
weights = masked_array(zeros(sh), mask = ones(sh))
areas   = masked_array(zeros(sh), mask = ones(sh))
for i in range(ncr):
    weights[:, :, i] = cals[crops[i]]
    areas[:, :, i]   = careas[crops[i]][0]
weights = masked_where(varr.mask, weights) # mask
areas   = masked_where(varr.mask, areas)

# convert to Pcal
wvarr = (varr * weights * areas).sum(axis = 2) / 1e6 # Gcal -> Pcal

bps       = [0] * 3
medianarr = zeros((3, nd))
maxarr    = zeros((3, nd))
minarr    = zeros((3, nd))

plt.figure(figsize = (10, 5))
ax = plt.subplot(111)

# hadgem noco2
p1           = wvarr[:, hadgemidx, :, 1]
medianarr[0] = median(p1, axis = 0)
maxarr[0]    = p1.max(axis = 0)
minarr[0]    = p1.min(axis = 0)
bp1 = ax.boxplot(p1, positions = range(1, 4 * nd, 4))
bps[0] = bp1

# hadgem co2
p2           = wvarr[:, hadgemidx, :, 0]
medianarr[1] = median(p2, axis = 0)
maxarr[1]    = p2.max(axis = 0)
minarr[1]    = p2.min(axis = 0)
bp2 = ax.boxplot(p2, positions = range(2, 4 * nd, 4))
bps[1] = bp2

# all co2
p3           = reshape(wvarr[:, :, :, 0], (nm * ng, nd))
medianarr[2] = median(p3, axis = 0)
maxarr[2]    = p3.max(axis = 0)
minarr[2]    = p3.min(axis = 0)
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
plt.ylim([2000, 10000])
plt.yticks(arange(2000, 11000, 1000))
plt.ylabel('Pcal')
plt.grid(which = 'major', axis = 'y')
plt.tight_layout()

l1, = plt.plot([1, 1], 'r')
l2, = plt.plot([1, 1], 'y')
l3, = plt.plot([1, 1], 'b')
plt.legend((l1, l2, l3), ['no co2 (6)', 'co2 HadGEM2-ES (6)', 'co2 (30)'], loc = 'upper left')
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

    f.createDimension('scen', len(medianarr))
    svar = f.createVariable('scen', 'i4', 'scen')
    svar[:] = range(1, 1 + len(medianarr))
    svar.units = 'mapping'
    svar.long_name = 'HadGEM-noCO2, HadGEM-CO2, CO2'

    bvar = f.createVariable(variable, 'f8', ('models', 'gcms', 'decade', 'co2'), zlib = True, shuffle = False, complevel = 9, fill_value = 1e20)
    bvar[:] = wvarr
    bvar.long_name = '%s at global level' % variable

    medianvar = f.createVariable('median_%s' % variable, 'f8', ('scen', 'decade'), zlib = True, shuffle = False, complevel = 9, fill_value = 1e20)
    medianvar[:] = medianarr
    medianvar.long_name = 'ensemble median of %s at global level' % variable

    maxvar = f.createVariable('max_%s' % variable, 'f8', ('scen', 'decade'), zlib = True, shuffle = False, complevel = 9, fill_value = 1e20)
    maxvar[:] = maxarr
    maxvar.long_name = 'ensemble maximum of %s at global level' % variable

    minvar = f.createVariable('min_%s' % variable, 'f8', ('scen', 'decade'), zlib = True, shuffle = False, complevel = 9, fill_value = 1e20)
    minvar[:] = minarr
    minvar.long_name = 'ensemble minimum of %s at global level' % variable