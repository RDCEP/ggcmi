#!/usr/bin/env python

# add paths
import os, sys
for p in os.environ['PATH'].split(':'): sys.path.append(p)

# import modules
from optparse import OptionParser
from netCDF4 import Dataset as nc
from averager import MeanAverager
from numpy.ma import masked_where
from numpy import logical_not, where
from aggmaskloader import AggMaskLoader

def createAggFile(filename, time, tunits, adata, anames, aunits, alongnames, scens, irr, leaddim, hasscen):
    if leaddim == 'scen':
        nscens = None
        ntime  = len(time)
    else:
        nscens = len(scens) if hasscen else 0
        ntime  = None

    with nc(filename, 'w', format = 'NETCDF4_CLASSIC') as f:
        f.createDimension('time', ntime)
        timevar = f.createVariable('time', 'i4', 'time')
        timevar[:] = time
        timevar.units = tunits
        timevar.long_name = 'time'

        if hasscen:
            f.createDimension('scen', nscens)
            scenvar = f.createVariable('scen', 'i4', 'scen')
            scenvar[:] = scens
            scenvar.units = 'no'
            scenvar.long_name = 'scenarios'

        f.createDimension('irr', len(irr) + 1)
        irrvar = f.createVariable('irr', 'i4', 'irr')
        irrvar[:] = range(1, len(irr) + 2)
        irrvar.units = 'mapping'
        irrvar.long_name = ', '.join(irr + ['sum'])

        for i in range(len(anames)):
            rname = str(anames[i])
            f.createDimension(rname, len(adata[i]))
            rvar = f.createVariable(rname, 'i4', rname)
            rvar[:] = adata[i]
            rvar.units = aunits[i]
            rvar.long_name = alongnames[i]

def getyieldmask(yieldvar, yieldthr1 = 0.1, yieldthr2 = 30):
    yieldm = masked_where(yieldvar < yieldthr1, yieldvar)
    yieldm = masked_where(yieldm > yieldthr2, yieldm)
    return logical_not(yieldm.mask)

# parse inputs
parser = OptionParser()
parser.add_option("-i", "--input", dest = "input", default = "", type = "string",
                  help = "File to aggregate:var 1,var 2,...,var M")
parser.add_option("-w", "--weights", dest = "weights", default = None, type = "string",
                  help = "Weights file (optional)", metavar = "FILE")
parser.add_option("-a", "--agg", dest = "agg", default = "", type = "string",
                  help = "Aggregation file")
parser.add_option("-n", "--numchunks", dest = "numchunks", default = 1, type = int,
                  help = "Number of chunks to split the data into (default = 1)")
parser.add_option("-s", "--scen", dest = "scen", default = None, type = "string",
                  help = "Comma-separated list of scenarios to aggregate (default = all)")
parser.add_option("-l", "--leaddim", dest = "leaddim", default = "scen", type = "string",
                  help = "Lead dimension of output file (default = scen)")
parser.add_option("-o", "--output", dest = "output", default = "", type = "string",
                  help = "Output file", metavar = "FILE")
parser.add_option("-y", "--yieldvar", dest = "yieldvar", default = None, type = "string",
                  help = "Name of yield variable (optional)")
parser.add_option("--calc_area", action = "store_true", dest = "calcarea", default = False,
                  help = "Flag to indicate weights are fractions (optional)")
options, args = parser.parse_args()

inputf    = options.input
weightsf  = options.weights
aggf      = options.agg
numchunks = options.numchunks
scen      = options.scen
leaddim   = options.leaddim
outputf   = options.output
yieldvar  = options.yieldvar
calcarea  = options.calcarea

if not leaddim in ['scen', 'time']: raise Exception('Unknown lead dimension')

yieldthr1, yieldthr2 = 0.1, 33 # yield thresholds (t/ha)

ifile, ivars = [i.strip() for i in inputf.split(':')]
ivars = [v.strip() for v in ivars.split(',')]
nvars = len(ivars)
var = [0] * nvars; dims = [0] * nvars; vunits = [0] * nvars
with nc(ifile) as f:
    lats, lons = f.variables['lat'][:], f.variables['lon'][:]

    t = f.variables['time']
    time = t[:]
    tunits = t.units if 'units' in t.ncattrs() else ''

    irr = f.variables['irr'].long_name.split(', ')

    hasscen = 'scen' in f.variables
    if hasscen:
        scenall = f.variables['scen'][:]
        if scen is None:
            scensel = scenall.copy()
            scenidx = range(len(scenall))
        else:
            scensel = [int(s) for s in scen.split(',')]
            scenidx = [0] * len(scensel)
            for i in range(len(scensel)):
                scenidx[i] = where(scenall == scensel[i])[0][0]
    else:
        scensel = [0] # single scenario
        scenidx = None

    for i in range(nvars):
        v = f.variables[ivars[i]]
        dims[i] = v.dimensions

        slicer = [slice(0, n) for n in v.shape]
        if hasscen: slicer[dims[i].index('scen')] = scenidx # select scenarios
        var[i] = v[slicer]

        vunits[i] = v.units if 'units' in v.ncattrs() else ''

    if yieldvar: # pull yield
        yieldv = f.variables[yieldvar]
        slicer = [slice(0, n) for n in yieldv.shape]
        if hasscen: slicer[yieldv.dimensions.index('scen')] = scenidx
        yieldv = yieldv[slicer]

iridx, rfidx = irr.index('ir'), irr.index('rf')

if weightsf:
    with nc(weightsf) as f:
        rfweights = f.variables['rainfed'][:]
        irweights = f.variables['irrigated'][:]
else:
    rfweights = irweights = None

aggloader  = AggMaskLoader(aggf, lats = lats, lons = lons)
adata      = aggloader.data()
audata     = aggloader.udata()
anames     = aggloader.names()
aunits     = aggloader.units()
alongnames = aggloader.longnames()

createAggFile(outputf, time, tunits, audata, anames, aunits, alongnames, scensel, irr, leaddim, hasscen)
f = nc(outputf, 'a')

avobj = MeanAverager()
for i in range(len(audata)):
    if leaddim == 'scen':
        dimsv = ('scen', 'time', anames[i], 'irr')
    elif hasscen:
        dimsv = ('time', 'scen', anames[i], 'irr')
    else:
        dimsv = ('time', anames[i], 'irr')

    for j in range(nvars):
        avev = f.createVariable(ivars[j] + '_' + anames[i], 'f4', dimsv, fill_value = 1e20, zlib = True, complevel = 9)
        avev.units = vunits[j]
        avev.long_name = 'average ' + anames[i] + ' ' + ivars[j]

        slicer = [slice(0, n) for n in var[j].shape]
        dimslist = list(dims[j])
        if hasscen: dimslist.remove('scen')
        dimslist.remove('irr')
        timeidx, latidx, lonidx = dimslist.index('time'), dimslist.index('lat'), dimslist.index('lon')

        for k in range(len(scensel)):
            if hasscen: slicer[dims[j].index('scen')] = k

            slicer[dims[j].index('irr')] = iridx
            irscen = var[j][slicer].transpose((timeidx, latidx, lonidx))
            slicer[dims[j].index('irr')] = rfidx
            rfscen = var[j][slicer].transpose((timeidx, latidx, lonidx))

            if yieldvar:
                slicer[dims[j].index('irr')] = iridx
                yir = yieldv[slicer].transpose((timeidx, latidx, lonidx))
                slicer[dims[j].index('irr')] = rfidx
                yrf = yieldv[slicer].transpose((timeidx, latidx, lonidx))
                ymaskir = getyieldmask(yir, yieldthr1, yieldthr2)
                ymaskrf = getyieldmask(yrf, yieldthr1, yieldthr2)
            else:
                ymaskir = ymaskrf = None

            if iridx > rfidx:
                data = avobj.combine(rfscen, irscen, adata[i], lats, rfweights, irweights, calcarea = calcarea, mask1 = ymaskrf, mask2 = ymaskir, numchunks = numchunks)
            else:
                data = avobj.combine(irscen, rfscen, adata[i], lats, irweights, rfweights, calcarea = calcarea, mask1 = ymaskir, mask2 = ymaskrf, numchunks = numchunks)

            if leaddim == 'scen':
                avev[k, :, :, :] = data.transpose((1, 0, 2))
            elif hasscen:
                avev[:, k, :, :] = data.transpose((1, 0, 2))
            else:
                avev[:] = data.transpose((1, 0, 2))

f.close()