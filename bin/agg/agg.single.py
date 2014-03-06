#!/usr/bin/env python

# import modules
from optparse import OptionParser
from netCDF4 import Dataset as nc
from aggmaskloader import AggMaskLoader
from averager import SumAverager, MeanAverager

def createAggFile(filename, time, tunits, adata, anames, aunits, alongnames):
    with nc(filename, 'w', format = 'NETCDF4_CLASSIC') as f:
        f.createDimension('time', len(time))
        timevar = f.createVariable('time', 'i4', ('time',))
        timevar[:] = time
        timevar.units = tunits
        timevar.long_name = 'time'
        for i in range(len(anames)):
            rname = anames[i] + '_index'
            f.createDimension(rname, len(adata[i]))
            rvar = f.createVariable(rname, 'i4', (rname,))
            rvar[:] = adata[i]
            rvar.units = aunits[i]
            rvar.long_name = alongnames[i]

# parse inputs
parser = OptionParser()
parser.add_option("-i", "--input", dest = "input", default = "", type = "string",
                  help = "File to aggregate: var1, ..., varN")
parser.add_option("-w", "--weights", dest = "weights", default = "none", type = "string",
                  help = "Weights file: var (none = all weights are one)")
parser.add_option("-a", "--agg", dest = "agg", default = "", type = "string",
                  help = "Aggregation file: var1, ..., varN")
parser.add_option("-t", "--type", dest = "type", default = "mean", type = "string",
                  help = "Aggregation type (mean or sum)")                  
parser.add_option("-o", "--output", dest = "output", default = "", type = "string",
                  help = "Output file", metavar = "FILE")
options, args = parser.parse_args()

if not options.type in ['mean', 'sum']: raise Exception('Invalid type')

weights = None # load weights
if not options.weights in ['none', 'None', '']:
    wfile, wvar = [w.strip() for w in options.weights.split(':')]
    with nc(wfile) as f: weights = f.variables[wvar][:]

afile, avars = [a.strip() for a in options.agg.split(':')] # load aggregation masks
avars = [v.strip() for v in avars.split(',')]
aggloader  = AggMaskLoader(afile, avars)
adata      = aggloader.data()
audata     = aggloader.udata()
anames     = aggloader.names()
aunits     = aggloader.units()
alongnames = aggloader.longnames()
nmasks     = len(audata)

ifile, ivars = [i.strip() for i in options.input.split(':')] # load input data
ivars = [v.strip() for v in ivars.split(',')]
nvars = len(ivars)
var = [0] * nvars; vunits = [0] * nvars
with nc(ifile) as f:
    lats = f.variables['lat'][:]
    t = f.variables['time']
    time = t[:]
    tunits = t.units if 'units' in t.ncattrs() else ''
    for i in range(nvars):
        v = f.variables[ivars[i]]
        var[i] = v[:]
        vunits[i] = v.units if 'units' in v.ncattrs() else ''

createAggFile(options.output, time, tunits, audata, anames, aunits, alongnames)
f = nc(options.output, 'a') # open file for appending

if options.type == 'mean': # averager instance
    avobj = MeanAverager()
else:
    avobj = SumAverager()

for i in range(nmasks):
    dims = (anames[i] + '_index', 'time')
    for j in range(nvars):
        avev = f.createVariable(ivars[j] + '_' + anames[i], 'f4', dims, fill_value = 1e20, zlib = True, complevel = 9)
        avev[:] = avobj.av(var[j], adata[i], lats, weights)
        avev.units = vunits[j]
        avev.long_name = ' '.join([options.type, anames[i], ivars[j]])
f.close()