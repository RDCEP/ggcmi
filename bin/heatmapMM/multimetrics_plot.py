#!/usr/bin/env python

from os import sep, listdir
from itertools import product
import matplotlib.pyplot as plt
from netCDF4 import Dataset as nc
from optparse import OptionParser
from collections import OrderedDict as od
from numpy.ma import masked_array, masked_where
from numpy import zeros, ones, double, where, resize, arange, unique, union1d

def str2num(arr):
    arr2 = arr[:]
    for i in range(len(arr2)):
        if arr2[i].isdigit(): arr2[i] = double(arr2[i])
    return arr2

def loadwts(weightfile, var, vals):
    if weightfile is None:
        return ones(len(vals))
    else:
        with nc(weightfile) as f:
            v = f.variables[var][:]
            w = f.variables['weights_' + var][:]
        nv = len(vals)
        wts = masked_array(zeros(nv), mask = ones(nv))
        for i in range(nv):
            idx = where(v == vals[i])[0]
            if idx.size: wts[i] = w[idx[0]]
        return wts

def condition(arr, dims, cross):
    """ arr is a numpy array, dims is an ordered dictionary of lists,
        cross is a tuple of (string, list (which can be * = all))
    """
    if not cross: return arr.copy(), dims.copy()

    var, vals = cross

    if not var in dims: raise Exception('Cannot slice along dimension')
    if vals == '*': return arr.copy(), dims.copy() # no change

    var_idx = dims.keys().index(var)
    dim_vals = dims[var]

    slice_in = [slice(0, size) for size in arr.shape]
    slice_out = slice_in[:] # copy

    sh = list(arr.shape)
    sh[var_idx] = len(vals)
    arr_out = masked_array(zeros(sh), mask = ones(sh))

    dims_out = dims.copy()
    dims_out[var] = vals

    for i in range(len(vals)):
        val = vals[i]
        if val in dim_vals:
            slice_in[var_idx] = [dim_vals.index(val)]
            slice_out[var_idx] = [i]
            arr_out[slice_out] = arr[slice_in]

    return arr_out, dims_out

def average(arr, dims, cross, wfile):
    """ arr is a numpy array, dims is an ordered dictionary of lists,
        cross is a tuple of (string, list (which can be * = all)), wfile is a filename string
    """
    arr_out, dims_out = condition(arr, dims, cross)

    if not cross: return arr.copy(), dims.copy()

    var, vals = cross

    if vals == '*': vals = dims_out[var]
    wts = loadwts(wfile, var, vals)

    var_idx = dims_out.keys().index(var)

    sh = arr_out.shape
    sh2 = list(sh)
    sh2[var_idx] = 1
    slice_arr = [slice(0, size) for size in arr_out.shape]

    wts_full = masked_array(zeros(sh), mask = ones(sh))
    for i in range(len(wts)):
        slice_arr[var_idx] = [i]
        wts_full[slice_arr] = resize(wts[i], sh2)

    wts_full = masked_where(arr_out.mask, wts_full) # mask

    arr_out = (wts_full * arr_out).mean(axis = var_idx) / wts_full.mean(axis = var_idx) # average over variable

    dims_out.pop(var); # remove variable from list

    return arr_out, dims_out

def optimize(arr, dims, cross, metric):
    """ arr is a numpy array, dims is an ordered dictionary of lists,
        cross is a tuple of (string, list (which can be * = all)), metric is a string
    """
    arr_out, dims_out = condition(arr, dims, cross)

    if not cross: return arr.copy(), dims.copy()

    var = cross[0]

    var_idx = dims_out.keys().index(var)

    if metric in ['tscorr', 'hitrate']:
        arr_out = arr_out.max(axis = var_idx) # maximize over variable
    else:
        arr_out = arr_out.min(axis = var_idx) # minimize over variable

    dims_out.pop(var);

    return arr_out, dims_out

def heatmap(arr, dims, x, y, metric, outdir, fmt, anon):
    vars = dims.keys()
    x_idx = vars.index(x); vars.remove(x)
    y_idx = vars.index(y); vars.remove(y)
    num_extra_dims = len(vars)
    if num_extra_dims:
        idx = [0] * num_extra_dims
        for i in range(num_extra_dims):
            idx[i] = range(0, len(dims[vars[i]]))
        for extra_dims in product(*idx):
            slice_arr = [slice(0, size) for size in arr.shape]
            title = '%s\n' % metric
            filename = '%s_x.%s_y.%s' % (metric, x, y)
            for i in range(len(extra_dims)):
                vidx = dims.keys().index(vars[i])
                slice_arr[vidx] = [extra_dims[i]]
                dim_val = str(dims[vars[i]][extra_dims[i]])
                title += '%s = %s, ' % (vars[i], dim_val)
                filename += '_%s' % dim_val
            arr_slice = arr[slice_arr].squeeze()
            if x_idx > y_idx: arr_slice = arr_slice.T # make x first
            ploth(dims[x], dims[y], arr_slice, x, y, metric, outdir + sep + filename + '.' + fmt, title = title[: -2], anon = anon)
    else:
        if x_idx > y_idx: arr = arr.T
        ploth(dims[x], dims[y], arr, x, y, metric, outdir + sep + '%s_x.%s_y.%s.%s' % (metric, x, y, fmt), title = metric, anon = anon)

def ploth(x, y, D, x_label, y_label, metric, filename, title = '', anon = False):
    if metric == 'tscorr':
        zmin, zmax = 0, 1
    elif metric in ['rmse', 'hitnumber']:
        zmin, zmax = 0, D.max()
    elif metric == 'varratio':
        zmin, zmax = 0, 2
    else:
        zmin, zmax = D.min(), D.max()
    fig, ax = plt.subplots()
    plt.pcolor(D.T, vmin = zmin, vmax = zmax)
    plt.axis([0, len(x), 0, len(y)])
    ax.set_xticks(arange(len(x)) + 0.5, minor = False)
    ax.set_yticks(arange(len(y)) + 0.5, minor = False)
    ax.set_xticks(arange(len(x)), minor = True)
    ax.set_yticks(arange(len(y)), minor = True)
    ax.invert_yaxis()
    if anon:
        ax.set_xticklabels(['%s %d' % (x_label, i) for i in range(1, len(x) + 1)], minor = False)
    else:
        ax.set_xticklabels(x, minor = False)
    ax.set_yticklabels(y, minor = False)
    plt.xticks(rotation = 90)
    plt.grid(which = 'minor', linestyle = '-')
    plt.tick_params(axis = 'both', which = 'major', bottom = 'off', top = 'off')
    plt.colorbar()
    plt.title(title, fontsize = 10)
    plt.savefig(filename, bbox_inches = 'tight')
    plt.close()

parser = OptionParser()
parser.add_option("-d", "--dir", dest = "dir", default = ".", type = "string",
                  help = "Directory where multimetrics files are located")
parser.add_option("-v", "--var", dest = "var", default = "rmse", type = "string",
                  help = "Variable to plot")
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
parser.add_option("-f", "--fmt", dest = "fmt", default = 'eps',
                  help = "Figure format (e.g., eps, png)")
parser.add_option("--outdir", dest = "outdir", default = ".", type = "string",
                  help = "Directory to save figures")
parser.add_option("--anon", action = "store_true", dest = "anon", default = False,
                  help = "Whether to anonymize x-axis labels")
options, args = parser.parse_args()

metric = options.var   # metric to plot
wfile  = options.wfile # weight file

# find all models, climates, and crops
files    = [f for f in listdir(options.dir) if not 'tscorr' in f and not 'rmse' in f]
models   = unique([f.split('_')[0] for f in files])
climates = unique([f.split('_')[1] for f in files])
crops    = unique([f.split('_')[3] for f in files])

nm, nw, ncp = len(models), len(climates), len(crops)

# x-axis
xsplit  = options.x.split(',')
x_var   = xsplit[0]
x_vals  = str2num(xsplit[1 :]) if len(xsplit) > 1 else '*'
x_tuple = (x_var, x_vals)

# y-axis
ysplit  = options.y.split(',')
y_var   = ysplit[0]
y_vals  = str2num(ysplit[1 :]) if len(ysplit) > 1 else '*'
y_tuple = (y_var, y_vals)

# conditioned dimensions
cross_tuples = []
if not options.cross is None:
    for c in options.cross:
        csplit = c.split(',')
        var, vals = csplit[0], str2num(csplit[1 :])
        cross_tuples.append((var, vals))

# averaged dimensions
ave_tuples = []
if not options.ave is None:
    for a in options.ave:
        asplit = a.split(',')
        var = asplit[0]
        vals = str2num(asplit[1 :]) if len(asplit) > 1 else '*'
        ave_tuples.append((var, vals))

# optimized dimensions
opt_tuples = []
if not options.opt is None:
    for o in options.opt:
        osplit = o.split(',')
        var = osplit[0]
        vals = str2num(osplit[1 :]) if len(osplit) > 1 else '*'
        opt_tuples.append((var, vals))

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
for m, w, c in product(range(nm), range(nw), range(ncp)):
    ftag = '%s_%s_hist_%s' % (models[m], climates[w], crops[c])
    file = [f for f in files if f.startswith(ftag)]

    if len(file) == 1: # unique file exists
        with nc(options.dir + sep + file[0]) as f:
            scens = union1d(scens, f.variables['scen'].long_name.split(', '))

dims['scen'] = scens

ng, ns        = len(dims['gadm0']), len(dims['scen'])
ndt, nmp, ncr = len(dims['dt']), len(dims['mp']), len(dims['cr'])
nt            = len(dims['time_range'])

# load data
sh = (nm, nw, ncp, ng, ns, ndt, nmp, ncr, nt)
data = masked_array(zeros(sh), mask = ones(sh))
for m, w, c in product(range(nm), range(nw), range(ncp)):
    ftag = '%s_%s_hist_%s' % (models[m], climates[w], crops[c])
    file = [f for f in files if f.startswith(ftag)]

    if len(file) == 1:
        with nc(options.dir + sep + file[0]) as f:
            scen = f.variables['scen'].long_name.split(', ')
            for s in range(len(scen)):
                sidx = where(scens == scen[s])[0][0]
                data[m, w, c, :, sidx] = f.variables[metric][:, s]

# condition on x-axis
data, dims = condition(data, dims, x_tuple)

# condition on y-axis
data, dims = condition(data, dims, y_tuple)

# condition on variables
for c in cross_tuples:
    data, dims = condition(data, dims, c)

# average over variables
for a in ave_tuples:
    data, dims = average(data, dims, a, wfile)

# optimize over variables
for o in opt_tuples:
    data, dims = optimize(data, dims, o, metric)

# plot data
heatmap(data, dims, x_var, y_var, metric, options.outdir, options.fmt, options.anon)