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
        f.createDimension('irr', 3)
        irrvar = f.createVariable('irr', 'i4', ('irr',))
        irrvar[:] = range(1, 4)
        irrvar.units = 'mapping'
        irrvar.long_name = 'ir, rf, sum'
        for i in range(len(anames)):
            rname = anames[i] + '_index'
            f.createDimension(rname, len(adata[i]))
            rvar = f.createVariable(rname, 'i4', (rname,))
            rvar[:] = adata[i]
            rvar.units = aunits[i]
            rvar.long_name = alongnames[i]

# parse inputs
parser = OptionParser()
parser.add_option("--inputir", dest = "inputir", default = "", type = "string",
                  help = "Irrigated file to aggregate: var")
parser.add_option("--inputrf", dest = "inputrf", default = "", type = "string",
                  help = "Rainfed file to aggregate: var")
parser.add_option("--weightsir", dest = "weightsir", default = "none", type = "string",
                  help = "Irrigated weights file: var (none = all weights are one)")
parser.add_option("--weightsrf", dest = "weightsrf", default = "none", type = "string",
                  help = "Rainfed weights file: var (none = all weights are one)")
parser.add_option("-a", "--agg", dest = "agg", default = "", type = "string",
                  help = "Aggregation file: var1, ..., varN")
parser.add_option("-t", "--type", dest = "type", default = "mean", type = "string",
                  help = "Aggregation type (mean or sum)")                  
parser.add_option("-o", "--output", dest = "output", default = "", type = "string",
                  help = "Output file", metavar = "FILE")
options, args = parser.parse_args()

if not options.type in ['mean', 'sum']: raise Exception('Invalid type')

weightsir = None # load weights
if not options.weightsir in ['none', 'None', '']:
    wfile, wvar = [w.strip() for w in options.weightsir.split(':')]
    with nc(wfile) as f: weightsir = f.variables[wvar][:]
weightsrf = None
if not options.weightsrf in ['none', 'None', '']:
    wfile, wvar = [w.strip() for w in options.weightsrf.split(':')]
    with nc(wfile) as f: weightsrf = f.variables[wvar][:]

afile, avars = [a.strip() for a in options.agg.split(':')] # load aggregation masks
avars = [v.strip() for v in avars.split(',')]
aggloader  = AggMaskLoader(afile, avars)
adata      = aggloader.data()
audata     = aggloader.udata()
anames     = aggloader.names()
aunits     = aggloader.units()
alongnames = aggloader.longnames()
nmasks     = len(audata)

irfile, irvarname = [i.strip() for i in options.inputir.split(':')] # load input data
with nc(irfile) as f:
    lats = f.variables['lat'][:]
    t = f.variables['time']
    time = t[:]
    tunits = t.units if 'units' in t.ncattrs() else ''
    v = f.variables[irvarname]
    varir = v[:]
    vunits = v.units if 'units' in v.ncattrs() else ''
rffile, rfvarname = [i.strip() for i in options.inputrf.split(':')]
with nc(rffile) as f:
    varrf = f.variables[rfvarname][:]

varname = irvarname.split('_')[0]

createAggFile(options.output, time, tunits, audata, anames, aunits, alongnames)
f = nc(options.output, 'a') # open file for appending

if options.type == 'mean': # averager instance
    avobj = MeanAverager()
else:
    avobj = SumAverager()

for i in range(nmasks):
    dims = (anames[i] + '_index', 'time', 'irr')
    avev = f.createVariable(varname + '_' + anames[i], 'f4', dims, fill_value = 1e20, zlib = True, complevel = 9)
    avev[:] = avobj.combine(varir, varrf, adata[i], lats, weightsir, weightsrf)
    avev.units = vunits
    avev.long_name = ' '.join([options.type, anames[i], varname])
f.close()