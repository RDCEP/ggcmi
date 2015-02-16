#!/usr/bin/env python

# import modules
from re import findall
from numpy import where
from optparse import OptionParser
from netCDF4 import Dataset as nc
from numpy.ma import reshape, resize

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

# absolute global yield
global26 = reshape(global26, (ng, nd, 10)).mean(axis = 2)
global85 = reshape(global85, (ng, nd, 10)).mean(axis = 2)

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

    df26var = f.createVariable('delta_yield_26_fpu', 'f8', ('fpu', 'decade'))
    df26var[:] = dyfpu26
    df26var.units = 't ha-1 yr-1'
    df26var.long_name = 'delta between average decadal yield and average historical yield with RCP 2.6 at fpu level'

    df85var = f.createVariable('delta_yield_85_fpu', 'f8', ('fpu', 'decade'))
    df85var[:] = dyfpu85
    df85var.units = 't ha-1 yr-1'
    df85var.long_name = 'delta between average decadal yield and average historical yield with RCP 8.5 at fpu level'

    g26var = f.createVariable('yield_26_global', 'f8', ('global', 'decade'))
    g26var[:] = global26
    g26var.units = 't ha-1 yr-1'
    g26var.long_name = 'average decadal yield with RCP 2.6 at global level'

    g85var = f.createVariable('yield_85_global', 'f8', ('global', 'decade'))
    g85var[:] = global85
    g85var.units = 't ha-1 yr-1'
    g85var.long_name = 'average decadal yield with RCP 8.5 at global level'