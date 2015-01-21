#!/usr/bin/env python

# import modules
import matplotlib
from itertools import product
import matplotlib.pyplot as plt
from optparse import OptionParser
from netCDF4 import Dataset as nc
from numpy.ma import masked_array
from mpl_toolkits.basemap import Basemap
from numpy import zeros, ones, resize, meshgrid, arange

# parse inputs
parser = OptionParser()
parser.add_option("-i", "--infile", dest = "infile", default = "metrics.nc4", type = "string",
                  help = "Metrics file", metavar = "FILE")
parser.add_option("-c", "--crop", dest = "crop", default = "maize", type = "string",
                  help = "Crop (or all)")
parser.add_option("-a", "--aggfile", dest = "aggfile", default = "fpu.mask.nc4", type = "string",
                  help = "Aggregation file", metavar = "FILE")
parser.add_option("--noco2", action = "store_true", dest = "noco2", default = False,
                  help = "Flag to indicate to plot noco2 scenario")
parser.add_option("-o", "--outfile", dest = "outfile", default = "", type = "string",
                  help = "Output plot file", metavar = "FILE")
options, args = parser.parse_args()

infile  = options.infile
crop    = options.crop
aggfile = options.aggfile
noco2   = options.noco2
outfile = options.outfile

cals = {'maize': 1, 'wheat': 1, 'soy': 1, 'rice': 1}

with nc(infile) as f:
    fpu = f.variables['fpu'][:]

with nc(aggfile) as f:
    lats, lons = f.variables['lat'][:], f.variables['lon'][:]
    fpumap     = f.variables['fpu'][:]

models = ['epic', 'gepic', 'image_leitap', 'lpj-guess', 'lpjml', 'pdssat', 'pegasus']
gcms   = ['gfdl-esm2m', 'hadgem2-es', 'ipsl-cm5a-lr', 'miroc-esm-chem', 'noresm1-m']
crops  = ['maize', 'wheat', 'soy', 'rice'] if crop == 'all' else [crop]
co2    = 'noco2' if noco2 else 'co2'

nm, ng, ncr, nfpu = len(models), len(gcms), len(crops), len(fpu)

# benefit
sh = (nm, ng, ncr, 3, nfpu)
benefit = masked_array(zeros(sh), mask = ones(sh))
with nc(infile) as f:
    for m, g, c in product(range(nm), range(ng), range(ncr)):
        var = 'benefit_fpu_%s_%s_%s_%s' % (models[m], gcms[g], crops[c], co2)
        if var in f.variables:
            benefit[m, g, c] = f.variables[var][-3 :, :, 2] # last three decades, sum

# weights
weights = masked_array(zeros(sh), mask = ones(sh))
for i in range(ncr):
    weights[:, :, i] = resize(cals[crops[i]], (nm, ng, 3, nfpu))

# average
wbenefit  = (benefit * weights).sum(axis = 3).sum(axis = 2).sum(axis = 1).sum(axis = 0)
wbenefit /= weights.sum(axis = 3).sum(axis = 2).sum(axis = 1).sum(axis = 0)

# rasterize
benefitmap = masked_array(zeros((len(lats), len(lons))), mask = ones((len(lats), len(lons))))
for i in range(nfpu):
    benefitmap[fpumap == fpu[i]] = wbenefit[i]

# plot
m = Basemap(llcrnrlon = -180, llcrnrlat = -90, urcrnrlon = 180, urcrnrlat = 90, \
            resolution = 'c', projection = 'cyl')
glon, glat = meshgrid(lons, lats)
x, y = m(glon, glat)
cs = m.pcolor(x, y, benefitmap, vmin = -1, vmax = 0.6, cmap = matplotlib.cm.jet)
cbar = m.colorbar(cs, location = 'right')
cbar.set_ticks(arange(-1, 0.7, 0.2))
m.drawcoastlines()
m.drawstates(linewidth = 0.2)
m.drawmapboundary()
m.drawcountries(linewidth = 0.2)
m.drawparallels(arange(90, -110, -30), labels = [1, 0, 0, 0])
m.drawmeridians(arange(-180, 180, 60), labels = [0, 0, 0, 1])

# save
plt.savefig(outfile)
plt.close()

#with nc('test.nc4', 'w') as f:
#    f.createDimension('lat', len(lats))
#    latvar = f.createVariable('lat', 'f8', 'lat')
#    latvar[:] = lats
#    latvar.units = 'degrees_north'
#    latvar.long_name = 'latitude'
#
#    f.createDimension('lon', len(lons))
#    lonvar = f.createVariable('lon', 'f8', 'lon')
#    lonvar[:] = lons
#    lonvar.units = 'degrees_east'
#    lonvar.long_name = 'longitude'
#
#    bfvar = f.createVariable('benefit', 'f8', ('lat', 'lon'), zlib = True, shuffle = False, complevel = 9, fill_value = 1e20)
#    bfvar[:] = benefitmap