#!/usr/bin/env python

# This script only supports reference by value (not index)

# top maize gadm0_index: 240, 48, 32, 105, 145, 163, 106, 225, 11, 237

# import cPickle as pickle
import matplotlib.pyplot as plt
from os import sep, listdir
from itertools import product
from netCDF4 import Dataset as nc
from optparse import OptionParser
from collections import OrderedDict as od
from numpy.ma import masked_array, masked_where, concatenate
from numpy import zeros, ones, double, newaxis, where, resize, arange

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

    if metric in ['tscorr', 'hitnumber']:
        arr_out = arr_out.max(axis = var_idx) # maximize over variable
    else:
        arr_out = arr_out.min(axis = var_idx) # minimize over variable

    dims_out.pop(var);

    return arr_out, dims_out

def heatmap(arr, dims, x, y, metric):
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
            ploth(dims[x], dims[y], arr_slice, x, y, metric, filename + '.png', title = title[: -2])
    else:
        if x_idx > y_idx: arr = arr.T
        ploth(dims[x], dims[y], arr, x, y, metric, '%s_x.%s_y.%s.png' % (metric, x, y), title = metric)

def ploth(x, y, D, x_label, y_label, metric, filename, title = ''):
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
    ax.set_xticklabels(x, minor = False) # ['Model %d' % i for i in range(1, len(x) + 1)], minor = False)
    ax.set_yticklabels(y, minor = False)
    plt.xticks(rotation = 90)
    plt.grid(which = 'minor', linestyle = '-')
    plt.tick_params(axis = 'both', which = 'major', bottom = 'off', top = 'off')
    plt.colorbar()
    plt.title(title, fontsize = 10)
    plt.savefig(filename, bbox_inches='tight')
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
options, args = parser.parse_args()

metric = options.var # metric to plot
wfile = options.wfile # weight file

# find all models
mfiles = [f for f in listdir(options.dir) if not 'tscorr' in f and not 'rmse' in f]
models = [f.split('_')[0] for f in mfiles]

model_option = ''
model_vals = ''

# parse options
xsplit = options.x.split(',')
x_var = xsplit[0]
x_vals = str2num(xsplit[1 :]) if len(xsplit) > 1 else '*'
x_tuple = ()
if x_var == 'model':
    model_vals = x_vals
else:
    x_tuple = (x_var, x_vals)

ysplit = options.y.split(',')
y_var = ysplit[0]
y_vals = str2num(ysplit[1 :]) if len(ysplit) > 1 else '*'
y_tuple = ()
if y_var == 'model':
    model_vals = y_vals
else:
    y_tuple = (y_var, y_vals)

cross_tuples = []
if not options.cross is None:
    for c in options.cross:
        csplit = c.split(',')
        var, vals = csplit[0], str2num(csplit[1 :])
        if var == 'model':
            model_vals = vals
        else:
            cross_tuples.append((var, vals))

ave_tuples = []
if not options.ave is None:
    for a in options.ave:
        asplit = a.split(',')
        var = asplit[0]
        vals = str2num(asplit[1 :]) if len(asplit) > 1 else '*'
        if var == 'model':
            model_option, model_vals = 'ave', vals
        else:
            ave_tuples.append((var, vals))

opt_tuples = []
if not options.opt is None:
    for o in options.opt:
        osplit = o.split(',')
        var = osplit[0]
        vals = str2num(osplit[1 :]) if len(osplit) > 1 else '*'
        if var == 'model':
            model_option, model_vals = 'opt', vals
        else:
            opt_tuples.append((var, vals))

if model_vals == '*': model_vals = models # use all models

# load data
for m in range(len(model_vals)):
    dims = od([])

    with nc(options.dir + sep + model_vals[m] + '_metrics.nc4') as f:
        dims['gadm0_index'] = list(f.variables['gadm0_index'][:])
        dims['scen'] = f.variables['scen'].long_name.split(', ')
        dims['detr_methods'] = f.variables['detr_methods'].long_name.split(', ')
        dims['corr_methods'] = f.variables['corr_methods'].long_name.split(', ')
        dims['time_range'] = f.variables['time_range'].long_name.split(', ')
        dims['climate'] = f.variables['climate'].long_name.split(', ')
        dims['crop'] = f.variables['crop'].long_name.split(', ')
        data = f.variables[metric][:]

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

    # append to big data array
    if not m:
        all_data = data[..., newaxis].copy()
    else:
        data_mod = data[..., newaxis].copy()
        all_data = concatenate((all_data, data_mod), axis = len(all_data.shape) - 1)

# add model to dimensions
dims['model'] = model_vals

# average or optimize over models if necessary
if model_option == 'ave':
    all_data, dims = average(all_data, dims, ('model', model_vals), wfile)
elif model_option == 'opt':
    all_data, dims = optimize(all_data, dims, ('model', model_vals), metric)

# plot data
# pickle.dump([all_data, dims, x_var, y_var], open('debug.p', 'w'))
heatmap(all_data, dims, x_var, y_var, metric)