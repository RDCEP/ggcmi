#!/usr/bin/env python

# import modules
import matplotlib
from os.path import splitext
from shapefile import Reader
from itertools import product
import matplotlib.pyplot as plt
from optparse import OptionParser
from netCDF4 import Dataset as nc
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import LineCollection
from numpy.ma import masked_array, masked_where, median, reshape
from numpy import zeros, ones, resize, meshgrid, arange, array, cos, pi, where

# parse inputs
parser = OptionParser()
parser.add_option("-i", "--infile", dest = "infile", default = "metrics.nc4", type = "string",
                  help = "Metrics file", metavar = "FILE")
parser.add_option("-c", "--crop", dest = "crop", default = "maize", type = "string",
                  help = "Crop (or all)")
parser.add_option("-a", "--aggfile", dest = "aggfile", default = "fpu.mask.nc4", type = "string",
                  help = "Aggregation file", metavar = "FILE")
parser.add_option("-r", "--hareafile", dest = "hareafile", default = "all.fpu.nc4", type = "string",
                  help = "Harvested area file", metavar = "FILE")
parser.add_option("-s", "--shapefile", dest = "shapefile", default = "fpu", type = "string",
                  help = "Shape file", metavar = "FILE")
parser.add_option("-w", "--weightfile", dest = "weightfile", default = "maize.nc4", type = "string",
                  help = "Weight file", metavar = "FILE")
parser.add_option("-p", "--percent", dest = "percent", default = "0.1", type = "float",
                  help = "Percent threshold")
parser.add_option("-m", "--mapfile", dest = "mapfile", default = "map.png", type = "string",
                  help = "Output map file", metavar = "FILE")
parser.add_option("-n", "--ncfile", dest = "ncfile", default = "map.nc4", type = "string",
                  help = "Output netcdf data file", metavar = "FILE")
options, args = parser.parse_args()

infile     = options.infile
crop       = options.crop
aggfile    = options.aggfile
hareafile  = options.hareafile
shapefile  = options.shapefile
weightfile = options.weightfile
percent    = options.percent
mapfile    = options.mapfile
ncfile     = options.ncfile

cals = {'maize': 3.60, 'wheat': 3.34, 'soy': 3.35, 'rice': 2.80}

with nc(infile) as f:
    fpu = f.variables['fpu'][:]
nfpu = len(fpu)

# load aggregation file
with nc(aggfile) as f:
    lats, lons = f.variables['lat'][:], f.variables['lon'][:]
    fpumap     = f.variables['fpu'][:]
nlats, nlons = len(lats), len(lons)

# load weight file
with nc(weightfile) as f:
    harea = f.variables['sum'][:]

# load area file
careas = {}
with nc(hareafile) as f:
    cfpu = f.variables['fpu'][:]
    for c in ['maize', 'wheat', 'soy', 'rice']:
        careas[c] = f.variables['area_' + c][:]

# find valid fpus
tarea = 100 * (111.2 / 2) ** 2 * cos(pi * lats / 180)
tarea = resize(tarea, (nlons, nlats)).T
validfpus = []
for i in range(nfpu):
    hareafpu = harea[fpumap == fpu[i]].sum()
    tareafpu = tarea[fpumap == fpu[i]].sum()
    if hareafpu / tareafpu > percent / 100.:
        validfpus.append(fpu[i])

# load shape file
r = Reader(shapefile)
shapes  = r.shapes()
records = r.records()

models = ['epic', 'gepic', 'lpj-guess', 'lpjml', 'pdssat', 'pegasus'] # exclude image
gcms   = ['gfdl-esm2m', 'hadgem2-es', 'ipsl-cm5a-lr', 'miroc-esm-chem', 'noresm1-m']
crops  = ['maize', 'wheat', 'soy', 'rice'] if crop == 'all' else [crop]
co2s   = ['co2', 'noco2']

nm, ng, ncr, nco2 = len(models), len(gcms), len(crops), len(co2s)

# variables
sh = (nm, ng, ncr, 3, nfpu, nco2)
dy26arr = masked_array(zeros(sh), mask = ones(sh))
dy85arr = masked_array(zeros(sh), mask = ones(sh))
with nc(infile) as f:
    for m, g, c, co in product(range(nm), range(ng), range(ncr), range(nco2)):
        var = 'delta_yield_26_fpu_%s_%s_%s_%s' % (models[m], gcms[g], crops[c], co2s[co])
        if var in f.variables:
            dy26arr[m, g, c, :, :, co] = f.variables[var][:, -3 :].T # last three decades
        var = 'delta_yield_85_fpu_%s_%s_%s_%s' % (models[m], gcms[g], crops[c], co2s[co])
        if var in f.variables:
            dy85arr[m, g, c, :, :, co] = f.variables[var][:, -3 :].T

# weights
weights = masked_array(zeros(sh), mask = ones(sh))
areas   = masked_array(zeros(sh), mask = ones(sh))
for i in range(ncr):
    weights[:, :, i] = cals[crops[i]]
    for f in range(nfpu):
        fpuidx = where(cfpu == fpu[f])[0][0]
        areas[:, :, i, :, f] = careas[crops[i]][fpuidx]
weights = masked_where(dy26arr.mask, weights) # mask
areas   = masked_where(dy26arr.mask, areas)

# average over crops and decades
dy26arr = (dy26arr * weights * areas).sum(axis = 3).sum(axis = 2) / areas.sum(axis = 3).sum(axis = 2)
dy85arr = (dy85arr * weights * areas).sum(axis = 3).sum(axis = 2) / areas.sum(axis = 3).sum(axis = 2)

hadgemidx = gcms.index('hadgem2-es')

barr = masked_array(zeros((3, nfpu)), mask = ones((3, nfpu)))
larr = masked_array(zeros((3, nfpu)), mask = ones((3, nfpu)))

# hadgem noco2
dy26m = median(dy26arr[:, hadgemidx, :, 1], axis = 0)
dy85m = median(dy85arr[:, hadgemidx, :, 1], axis = 0)
negy = dy85m < -0.01
barr[0, negy] = 100 * (1 - dy26m[negy] / dy85m[negy])
posy = dy85m > 0.01
larr[0, posy] = 100 * (dy26m[posy] / dy85m[posy] - 1)

# hadgem co2
dy26m = median(dy26arr[:, hadgemidx, :, 0], axis = 0)
dy85m = median(dy85arr[:, hadgemidx, :, 0], axis = 0)
negy = dy85m < -0.01
barr[1, negy] = 100 * (1 - dy26m[negy] / dy85m[negy])
posy = dy85m > 0.01
larr[1, posy] = 100 * (dy26m[posy] / dy85m[posy] - 1)

# all co2
dy26m = median(reshape(dy26arr[:, :, :, 0], (nm * ng, nfpu)), axis = 0)
dy85m = median(reshape(dy85arr[:, :, :, 0], (nm * ng, nfpu)), axis = 0)
negy = dy85m < -0.01
barr[2, negy] = 100 * (1 - dy26m[negy] / dy85m[negy])
posy = dy85m > 0.01
larr[2, posy] = 100 * (dy26m[posy] / dy85m[posy] - 1)

filename, ext = splitext(mapfile)
mapfiles = [filename + '.noco2' + ext, filename + '.co2.hadgem' + ext, filename + '.co2' + ext]
filename, ext = splitext(ncfile)
ncfiles = [filename + '.noco2' + ext, filename + '.co2.hadgem' + ext, filename + '.co2' + ext]

for i in range(len(barr)):
    # rasterize
    bmap = masked_array(zeros((nlats, nlons)), mask = ones((nlats, nlons)))
    lmap = masked_array(zeros((nlats, nlons)), mask = ones((nlats, nlons)))
    for j in range(len(validfpus)):
        fpuidx = where(fpu == validfpus[j])[0][0]
        bmap[fpumap == validfpus[j]] = barr[i, fpuidx]
        lmap[fpumap == validfpus[j]] = larr[i, fpuidx]

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
            for j in range(1, len(shape.parts)):
                index  = shape.parts[j - 1]
                index2 = shape.parts[j]
                segs.append(data[index : index2])
            segs.append(data[index2 :])

        lines = LineCollection(segs, antialiaseds = (1,))
        lines.set_edgecolors('k')
        lines.set_linewidth(0.1)
        ax.add_collection(lines)

    # plot variable map
    glon, glat = meshgrid(lons, lats)
    x, y = m(glon, glat)
    cs1 = m.pcolor(x, y, bmap, vmin = -100, vmax = 100, cmap = matplotlib.cm.seismic)
    cs2 = m.pcolor(x, y, lmap, vmin = -100, vmax = 0, cmap = matplotlib.cm.Greys_r)
    m.drawcoastlines()
    m.drawmapboundary()
    m.drawparallels(arange(90, -90, -30),  labels = [1, 0, 0, 0])
    m.drawmeridians(arange(-180, 180, 60), labels = [0, 0, 0, 1])

    # save
    plt.savefig(mapfiles[i])
    plt.close()

    # write data file
    with nc(ncfiles[i], 'w') as f:
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

        bmvar = f.createVariable('beta', 'f8', ('lat', 'lon'), zlib = True, shuffle = False, complevel = 9, fill_value = 1e20)
        bmvar[:] = bmap
        bmvar.long_name = '100 * (1 - delta yield RCP 2.6 / delta yield RCP 8.5)'

        lmvar = f.createVariable('lambda', 'f8', ('lat', 'lon'), zlib = True, shuffle = False, complevel = 9, fill_value = 1e20)
        lmvar[:] = lmap
        lmvar.long_name = '100 * (delta yield RCP 2.6 / delta yield RCP 8.5 - 1)'