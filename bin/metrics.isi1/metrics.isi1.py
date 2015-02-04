#!/usr/bin/env python

# import modules
from re import findall
from optparse import OptionParser
from netCDF4 import Dataset as nc
from numpy import where, zeros, ones
from numpy.ma import masked_array, reshape, resize

# parse inputs
parser = OptionParser()
parser.add_option("--rcp26file", dest = "rcp26file", default = "", type = "string",
                  help = "RCP 2.6 file", metavar = "FILE")
parser.add_option("--rcp85file", dest = "rcp85file", default = None, type = "string",
                  help = "RCP 8.5 file", metavar = "FILE")
parser.add_option("-v", "--variable", dest = "variable", default = "", type = "string",
                  help = "Variable to process")
parser.add_option("-o", "--outfile", dest = "outfile", default = "", type = "string",
                  help = "Output file", metavar = "FILE")
options, args = parser.parse_args()

rcp26file = options.rcp26file
rcp85file = options.rcp85file
variable  = options.variable
outfile   = options.outfile

with nc(rcp26file) as f:
    afpu, aglobal  = f.variables['fpu'][:], f.variables['global'][:]
    funits, flname = f.variables['fpu'].units, f.variables['fpu'].long_name
    gunits, glname = f.variables['global'].units, f.variables['global'].long_name

    time   = f.variables['time'][:]
    tunits = f.variables['time'].units
    time  += int(findall(r'\d+', tunits)[0])

    fpu26    = f.variables[variable + '_fpu'][:, :, 0, 0, 0, 0]    # fpu, time, scen, dt, mp, cr
    global26 = f.variables[variable + '_global'][:, :, 0, 0, 0, 0] # global, time, scen, dt, mp, cr

with nc(rcp85file) as f:
    fpu85    = f.variables[variable + '_fpu'][:, :, 0, 0, 0, 0]
    global85 = f.variables[variable + '_global'][:, :, 0, 0, 0, 0]

nt, nf, ng = len(time), len(afpu), len(aglobal)

tidx1, tidx2 = where(time == 1980)[0][0], where(time == 2009)[0][0] + 1

# number of decades
nd = nt / 10

# delta yield
dyfpu26 = reshape(fpu26, (nf, nd, 10)).mean(axis = 2) - resize(fpu26[:, tidx1 : tidx2].mean(axis = 1), (nd, nf)).T
dyfpu85 = reshape(fpu85, (nf, nd, 10)).mean(axis = 2) - resize(fpu85[:, tidx1 : tidx2].mean(axis = 1), (nd, nf)).T
dyglobal26 = reshape(global26, (ng, nd, 10)).mean(axis = 2) - resize(global26[:, tidx1 : tidx2].mean(axis = 1), (nd, ng)).T
dyglobal85 = reshape(global85, (ng, nd, 10)).mean(axis = 2) - resize(global85[:, tidx1 : tidx2].mean(axis = 1), (nd, ng)).T

# metrics
betafpu      = masked_array(zeros((nf, nd)), mask = ones((nf, nd)))
betaglobal   = masked_array(zeros((ng, nd)), mask = ones((ng, nd)))
lambdafpu    = masked_array(zeros((nf, nd)), mask = ones((nf, nd)))
lambdaglobal = masked_array(zeros((ng, nd)), mask = ones((ng, nd)))

# beta
negy             = dyfpu85 < -0.01
betafpu[negy]    = 100 * (1 - dyfpu26[negy] / dyfpu85[negy])
negy             = dyglobal85 < -0.01
betaglobal[negy] = 100 * (1 - dyglobal26[negy] / dyglobal85[negy])

# lambda
posy               = dyfpu85 > 0.01
lambdafpu[posy]    = 100 * (dyfpu26[posy] / dyfpu85[posy] - 1)
posy               = dyglobal85 > 0.01
lambdaglobal[posy] = 100 * (dyglobal26[posy] / dyglobal85[posy] - 1)

with nc(outfile, 'w') as f:
    f.createDimension('fpu', nf)
    fpuvar = f.createVariable('fpu', 'i4', 'fpu')
    fpuvar[:] = afpu
    fpuvar.units = funits
    fpuvar.long_name = flname

    f.createDimension('global', ng)
    globalvar = f.createVariable('global', 'i4', 'global')
    globalvar[:] = aglobal
    globalvar.units = gunits
    globalvar.long_name = glname

    f.createDimension('decade', nd)
    dvar = f.createVariable('decade', 'i4', 'decade')
    dvar[:] = range(nd)
    dvar.units = 'decades since from 1980'
    dvar.long_name = 'decade'

    bfvar = f.createVariable('beta_fpu', 'f8', ('fpu', 'decade'))
    bfvar[:] = betafpu
    bfvar.units = ''
    bfvar.long_name = '100 * (1 - delta yield rcp 2.6 / delta yield rcp 8.5) at fpu level'

    bgvar = f.createVariable('beta_global', 'f8', ('global', 'decade'))
    bgvar[:] = betaglobal
    bgvar.units = ''
    bgvar.long_name = '100 * (1 - delta yield rcp 2.6 / delta yield rcp 8.5) at global level'

    lfvar = f.createVariable('lambda_fpu', 'f8', ('fpu', 'decade'))
    lfvar[:] = lambdafpu
    lfvar.units = ''
    lfvar.long_name = '100 * (delta yield rcp 2.6 / delta yield rcp 8.5 - 1) at fpu level'

    lgvar = f.createVariable('lambda_global', 'f8', ('global', 'decade'))
    lgvar[:] = lambdaglobal
    lgvar.units = ''
    lgvar.long_name = '100 * (delta yield rcp 2.6 / delta yield rcp 8.5 - 1) at global level'