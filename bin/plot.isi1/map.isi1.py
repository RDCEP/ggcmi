#!/usr/bin/env python

# import modules
import matplotlib
from shapefile import Reader
from itertools import product
import matplotlib.pyplot as plt
from optparse import OptionParser
from netCDF4 import Dataset as nc
from numpy.ma import masked_array
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import LineCollection
from numpy import zeros, ones, resize, meshgrid, arange, array

# parse inputs
parser = OptionParser()
parser.add_option("-i", "--infile", dest = "infile", default = "metrics.nc4", type = "string",
                  help = "Metrics file", metavar = "FILE")
parser.add_option("-c", "--crop", dest = "crop", default = "maize", type = "string",
                  help = "Crop (or all)")
parser.add_option("-a", "--aggfile", dest = "aggfile", default = "fpu.mask.nc4", type = "string",
                  help = "Aggregation file", metavar = "FILE")
parser.add_option("-s", "--shapefile", dest = "shapefile", default = "fpu", type = "string",
                  help = "Shape file", metavar = "FILE")
parser.add_option("--co2", dest = "co2", default = "co2", type = "string",
                  help = "co2 setting (e.g., co2 or noco2)")
parser.add_option("-m", "--mapfile", dest = "mapfile", default = "map.png", type = "string",
                  help = "Output map file", metavar = "FILE")
parser.add_option("-n", "--ncfile", dest = "ncfile", default = "map.nc4", type = "string",
                  help = "Output netcdf data file", metavar = "FILE")
options, args = parser.parse_args()

infile    = options.infile
crop      = options.crop
aggfile   = options.aggfile
shapefile = options.shapefile
co2       = options.co2
mapfile   = options.mapfile
ncfile    = options.ncfile

cals = {'maize': 1, 'wheat': 1, 'soy': 1, 'rice': 1}

with nc(infile) as f:
    fpu = f.variables['fpu'][:]

with nc(aggfile) as f:
    lats, lons = f.variables['lat'][:], f.variables['lon'][:]
    fpumap     = f.variables['fpu'][:]

models = ['epic', 'gepic', 'image_leitap', 'lpj-guess', 'lpjml', 'pdssat', 'pegasus']
gcms   = ['gfdl-esm2m', 'hadgem2-es', 'ipsl-cm5a-lr', 'miroc-esm-chem', 'noresm1-m']
crops  = ['maize', 'wheat', 'soy', 'rice'] if crop == 'all' else [crop]

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

# load shape file
r = Reader(shapefile)
shapes  = r.shapes()
records = r.records()

# plot map and fpu boundaries
plt.figure()
ax = plt.subplot(111)
m = Basemap(llcrnrlon = -180, llcrnrlat = -60, urcrnrlon = 180, urcrnrlat = 90, \
            resolution = 'c', projection = 'cyl')
for record, shape in zip(records, shapes):
    slons, slats = zip(*shape.points)
    data = array(m(slons, slats)).T

    if len(shape.parts) == 1:
        segs = [data,]
    else:
        segs = []
        for i in range(1, len(shape.parts)):
            index  = shape.parts[i - 1]
            index2 = shape.parts[i]
            segs.append(data[index : index2])
        segs.append(data[index2 :])

    lines = LineCollection(segs, antialiaseds = (1,))
    lines.set_edgecolors('k')
    lines.set_linewidth(0.1)
    ax.add_collection(lines)

# plot benefit map
glon, glat = meshgrid(lons, lats)
x, y = m(glon, glat)
cs = m.pcolor(x, y, benefitmap, vmin = -0.5, vmax = 0.5, cmap = matplotlib.cm.RdYlGn)
cbar = m.colorbar(cs, location = 'right')
cbar.set_ticks(arange(-0.5, 0.6, 0.25))
m.drawcoastlines()
m.drawmapboundary()
m.drawparallels(arange(90, -90, -30), labels = [1, 0, 0, 0])
m.drawmeridians(arange(-180, 180, 60), labels = [0, 0, 0, 1])

# save
plt.savefig(mapfile)
plt.close()

with nc(ncfile, 'w') as f:
    f.createDimension('lat', len(lats))
    latvar = f.createVariable('lat', 'f8', 'lat')
    latvar[:] = lats
    latvar.units = 'degrees_north'
    latvar.long_name = 'latitude'

    f.createDimension('lon', len(lons))
    lonvar = f.createVariable('lon', 'f8', 'lon')
    lonvar[:] = lons
    lonvar.units = 'degrees_east'
    lonvar.long_name = 'longitude'

    f.createDimension('fpu', nfpu)
    fpuvar = f.createVariable('fpu', 'i4', 'fpu')
    fpuvar[:] = fpu
    fpuvar.units = 'FPU index'
    fpuvar.long_name = '309 Food Producing Units'

    bmvar = f.createVariable('benefit', 'f8', ('lat', 'lon'), zlib = True, shuffle = False, complevel = 9, fill_value = 1e20)
    bmvar[:] = benefitmap
    bmvar.long_name = 'benefit of RCP 2.6 over RCP 8.5'

    bfvar = f.createVariable('benefit_fpu', 'f8', 'fpu', zlib = True, shuffle = False, complevel = 9, fill_value = 1e20)
    bfvar[:] = wbenefit
    bfvar.long_name = 'benefit of RCP 2.6 over RCP 8.5 at fpu level'