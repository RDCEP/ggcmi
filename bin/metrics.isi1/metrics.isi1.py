#!/usr/bin/env python

# import modules
from re import findall
from optparse import OptionParser
from netCDF4 import Dataset as nc
from numpy import where, zeros, ones
from numpy.ma import resize, reshape, masked_array

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

    irr      = f.variables['irr'][:]
    irrlname = f.variables['irr'].long_name

    fpu26    = f.variables[variable + '_fpu'][:]
    global26 = f.variables[variable + '_global'][:]

with nc(rcp85file) as f:
    fpu85    = f.variables[variable + '_fpu'][:]
    global85 = f.variables[variable + '_global'][:]

nt, nf, ng, nirr = len(time), len(afpu), len(aglobal), len(irr)

tidx1, tidx2 = where(time == 1980)[0][0], where(time == 2009)[0][0] + 1

# delta yield
dyfpu = masked_array(zeros((nt, nf, nirr, 2)), mask = ones((nt, nf, nirr, 2)))
dyfpu[:, :, :, 0] = fpu26 / resize(fpu26[tidx1 : tidx2].mean(axis = 0), fpu26.shape)
dyfpu[:, :, :, 1] = fpu85 / resize(fpu85[tidx1 : tidx2].mean(axis = 0), fpu85.shape)

dyglobal = masked_array(zeros((nt, ng, nirr, 2)), mask = ones((nt, ng, nirr, 2)))
dyglobal[:, :, :, 0] = global26 / resize(global26[tidx1 : tidx2].mean(axis = 0), global26.shape)
dyglobal[:, :, :, 1] = global85 / resize(global85[tidx1 : tidx2].mean(axis = 0), global85.shape)

# benefit
nd = nt / 10
bfpu  = reshape(dyfpu[:, :, :, 1], (nd, 10, nf, nirr)).mean(axis = 1)
bfpu /= reshape(dyfpu[:, :, :, 0], (nd, 10, nf, nirr)).mean(axis = 1)
bfpu  = 1. - bfpu

bglobal  = reshape(dyglobal[:, :, :, 1], (nd, 10, ng, nirr)).mean(axis = 1)
bglobal /= reshape(dyglobal[:, :, :, 0], (nd, 10, ng, nirr)).mean(axis = 1)
bglobal  = 1. - bglobal

with nc(outfile, 'w') as f:
    f.createDimension('time', nt)
    timevar = f.createVariable('time', 'i4', 'time')
    timevar[:] = time - time[0]
    timevar.units = tunits
    timevar.long_name = 'time'

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

    f.createDimension('irr', nirr)
    irrvar = f.createVariable('irr', 'i4', 'irr')
    irrvar[:] = irr
    irrvar.units = 'mapping'
    irrvar.long_name = irrlname

    f.createDimension('rcp', 2)
    rcpvar = f.createVariable('rcp', 'i4', 'rcp')
    rcpvar[:] = range(2)
    rcpvar.units = 'mapping'
    rcpvar.long_name = 'rcp26, rcp85'

    f.createDimension('decade', nd)
    dvar = f.createVariable('decade', 'i4', 'decade')
    dvar[:] = range(nd)
    dvar.units = 'decades since from 1980'
    dvar.long_name = 'decade'

    dyfvar = f.createVariable('delta_yield_fpu', 'f8', ('time', 'fpu', 'irr', 'rcp'))
    dyfvar[:] = dyfpu
    dyfvar.units = ''
    dyfvar.long_name = 'delta yield relative to 1980-2009 average at fpu level'

    dygvar = f.createVariable('delta_yield_global', 'f8', ('time', 'global', 'irr', 'rcp'))
    dygvar[:] = dyglobal
    dygvar.units = ''
    dygvar.long_name = 'delta yield relative to 1980-2009 average at global level'

    bfvar = f.createVariable('benefit_fpu', 'f8', ('decade', 'fpu', 'irr'))
    bfvar[:] = bfpu
    bfvar.units = ''
    bfvar.long_name = 'benefit of RCP 2.6 over RCP 8.5 at fpu level'

    bgvar = f.createVariable('benefit_global', 'f8', ('decade', 'global', 'irr'))
    bgvar[:] = bglobal
    bgvar.units = ''
    bgvar.long_name = 'benefit of RCP 2.6 over RCP 8.5 at global level'