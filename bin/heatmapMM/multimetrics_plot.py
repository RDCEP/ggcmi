#!/usr/bin/env python

# add paths
import os, sys
for p in os.environ['PATH'].split(':'): sys.path.append(p)

from os import sep, listdir
from itertools import product
from mmplotter import MMPlotter
from numpy.ma import masked_array
from netCDF4 import Dataset as nc
from optparse import OptionParser
from collections import OrderedDict as od
from numpy import zeros, ones, where, unique, union1d

def parseArg(var):
    vsplit = var.split('/')

    nv = len(vsplit)

    vnames, vvals = [0] * nv, [0] * nv
    for i in range(nv):
        vs = vsplit[i].split(':')
        vnames[i] = vs[0]
        vvals[i]  = '*' if len(vs) == 1 else vs[1].split(',')

    vnames = '/'.join(vnames)

    return {vnames: vvals}

parser = OptionParser()
parser.add_option("-d", "--dir", dest = "dir", default = ".", type = "string",
                  help = "Directory where multimetrics files are located")
parser.add_option("-m", "--metric", dest = "metric", default = "rmse", type = "string",
                  help = "Metric to plot")
parser.add_option("-x", "--xaxis", dest = "x", default = "clim", type = "string",
                  help = "Variable to plot along x-axis (followed optionally by comma-separated list of values to plot)")
parser.add_option("-y", "--yaxis", dest = "y", default = "model", type = "string",
                  help = "Variable to plot along y-axis (followed optionally by comma-separated list of values to plot)")
parser.add_option("-c", "--cross", dest = "cross", action = "append",
                  help = "Variable and value(s) to condition plot on")
parser.add_option("-a", "--ave", dest = "ave", action = "append",
                  help = "Variable to average out (followed optionally by values to average over)")
parser.add_option("-o", "--opt", dest = "opt", action = "append",
                  help = "Variable to optimize over (followed optionally by values to optimize over)")
parser.add_option("-w", "--wfile", dest = "wfile", default = None,
                  help = "Weight file (default = None)")
parser.add_option("--nm", dest = "nm", default = 1, type = int,
                  help = "Number of models to use in ensemble (default = 1)")
parser.add_option("--wt", dest = "wt", default = 1, type = int,
                  help = "Whether to use unweighted (1) or weighted (2) ensemble")
parser.add_option("-f", "--fmt", dest = "fmt", default = 'eps',
                  help = "Figure format (e.g., eps, png, jpg)")
parser.add_option("--outdir", dest = "outdir", default = ".", type = "string",
                  help = "Directory to save figures")
parser.add_option("--anon", action = "store_true", dest = "anon", default = False,
                  help = "Whether to anonymize x-axis labels")
options, args = parser.parse_args()

# x-axis
x = parseArg(options.x)

# y-axis
y = parseArg(options.y)

# conditioned dimensions
c = {}
if not options.cross is None:
    for ci in options.cross:
        c.update(parseArg(ci))

# averaged dimensions
a = {}
if not options.ave is None:
    for ai in options.ave:
        a.update(parseArg(ai))

# optimized dimensions
o = {}
if not options.opt is None:
    for oi in options.opt:
        o.update(parseArg(oi))

# find all models, climates, and crops
files    = listdir(options.dir)
models   = unique([f.split('_')[0] for f in files])
climates = unique([f.split('_')[1] for f in files])
crops    = unique([f.split('_')[3] for f in files])

nm, nw, ncp = len(models), len(climates), len(crops)

dims = od([])
dims['model']   = list(models)
dims['climate'] = list(climates)
dims['crop']    = list(crops)

with nc(options.dir + sep + files[0]) as f:
    dims['gadm0']      = list(f.variables['gadm0'][:])
    dims['scen']       = []
    dims['dt']         = f.variables['dt'].long_name.split(', ')
    dims['mp']         = f.variables['mp'].long_name.split(', ')
    dims['cr']         = f.variables['cr'].long_name.split(', ')
    dims['time_range'] = f.variables['time_range'].long_name.split(', ')

scens = []
for m, w, cp in product(range(nm), range(nw), range(ncp)):
    ftag = '%s_%s_hist_%s' % (models[m], climates[w], crops[cp])
    file = [f for f in files if f.startswith(ftag)]

    if len(file) == 1: # unique file exists
        with nc(options.dir + sep + file[0]) as f:
            if models[m] in ['rmse', 'tscorr']:
                scens = union1d(scens, ['optimal'])
            else:
                scens = union1d(scens, f.variables['scen'].long_name.split(', '))

dims['scen'] = list(scens)

ng, ns        = len(dims['gadm0']), len(dims['scen'])
ndt, nmp, ncr = len(dims['dt']), len(dims['mp']), len(dims['cr'])
nt            = len(dims['time_range'])

# load data
sh = (nm, nw, ncp, ng, ns, ndt, nmp, ncr, nt)
data = masked_array(zeros(sh), mask = ones(sh))
for m, w, cp in product(range(nm), range(nw), range(ncp)):
    ftag = '%s_%s_hist_%s' % (models[m], climates[w], crops[cp])
    file = [f for f in files if f.startswith(ftag)]

    if len(file) == 1:
        with nc(options.dir + sep + file[0]) as f:
            if models[m] in ['rmse', 'tscorr']:
                sidx  = where(scens == 'optimal')[0][0]
                nmidx = where(f.variables['nm'][:] == options.nm)[0][0]
                wtidx = where(f.variables['wt'][:] == options.wt)[0][0]
                data[m, w, cp, :, sidx] = f.variables[options.metric][:, :, :, :, :, nmidx, wtidx]
            else:
                scen = f.variables['scen'].long_name.split(', ')
                for s in range(len(scen)):
                    sidx = where(scens == scen[s])[0][0]
                    data[m, w, cp, :, sidx] = f.variables[options.metric][:, s]

# change ensemble model names
dims['model'][dims['model'].index('rmse')]   = 'rmse-weighted'
dims['model'][dims['model'].index('tscorr')] = 'tscorr-weighted'

# plot data
plotter = MMPlotter(options.metric, options.wfile, options.outdir, options.fmt, options.anon)
plotter.plot(data, dims, x, y, c, a, o)