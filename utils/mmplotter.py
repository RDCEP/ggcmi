#!/usr/bin/env python

from os import sep
from itertools import product
import matplotlib.pyplot as plt
from netCDF4 import Dataset as nc
from numpy import zeros, ones, where, arange, double, array
from numpy.ma import masked_array, masked_where, resize, isMaskedArray

class MMPlotter(object):
    def __init__(self, metric, wfile, outdir, fmt, anon):
        self.metric = metric
        self.wfile  = wfile
        self.outdir = outdir
        self.fmt    = fmt
        self.anon   = anon

    def plot(self, data, dims, x, y, c, a, o):
        # condition on x-axis
        data, dims = self.__condition(data, dims, x)

        # condition on y-axis
        data, dims = self.__condition(data, dims, y)

        # condition on variables
        for key, value in c.iteritems():
            data, dims = self.__condition(data, dims, {key: value})

        # average over variables
        for key, value in a.iteritems():
            data, dims = self.__average(data, dims, {key: value})

        # optimize over variables
        for key, value in o.iteritems():
            data, dims = self.__optimize(data, dims, {key: value})

        self.__plotHeatmap(data, dims, x.keys()[0], y.keys()[0])

    def __plotHeatmap(self, data, dims, x, y):
        vars = dims.keys()

        x_idx = vars.index(x); vars.remove(x)
        y_idx = vars.index(y); vars.remove(y)

        num_extra_dims = len(vars)
        if num_extra_dims:
            idx = [0] * num_extra_dims
            for i in range(num_extra_dims):
                idx[i] = range(0, len(dims[vars[i]]))

            for extra_dims in product(*idx):
                slice_arr = [slice(0, size) for size in data.shape]
                title = '%s\n' % self.metric
                filename = '%s_x.%s_y.%s' % (self.metric, x.replace('/', ','), y.replace('/', ','))

                for i in range(len(extra_dims)):
                    vidx = dims.keys().index(vars[i])
                    slice_arr[vidx] = [extra_dims[i]]
                    dim_val = str(dims[vars[i]][extra_dims[i]])
                    title += '%s = %s, ' % (vars[i], dim_val)
                    filename += '_%s' % dim_val

                arr_slice = data[slice_arr].squeeze()
                if x_idx > y_idx: arr_slice = arr_slice.T # make x first
                filename = self.outdir + sep + filename + '.' + self.fmt
                self.__ploth(dims[x], dims[y], arr_slice, x, y, filename, title = title[: -2] + '\n')
        else:
            if x_idx > y_idx: arr = data.T
            filename = self.outdir + sep + '%s_x.%s_y.%s.%s' % (self.metric, x, y, self.fmt)
            self.__ploth(dims[x], dims[y], arr, x, y, filename, title = self.metric + '\n')

    def __ploth(self, x, y, D, x_label, y_label, filename, title = ''):
        x2 = x[:]
        y2 = y[:]

        if self.metric == 'tscorr':
            zmin, zmax = 0, 1
        elif self.metric in ['rmse', 'hitnumber']:
            zmin, zmax = 0, D.max()
        elif self.metric == 'varratio':
            zmin, zmax = 0, 2
        else:
            zmin, zmax = D.min(), D.max()

        # compute best rows
        D2 = masked_array(zeros((len(x) + 1, len(y) + 1)), mask = ones((len(x) + 1, len(y) + 1)))
        D2[: len(x), : len(y)] = D
        if self.metric in ['tscorr', 'hitrate']: # max
            D2[len(x), : len(y)] = D.max(axis = 0)
            D2[:, len(y)]        = D2.max(axis = 1)
        else: # min
            D2[len(x), : len(y)] = D.min(axis = 0)
            D2[:, len(y)]        = D2.min(axis = 1)
        D   = D2.copy()
        x2 += ['best']
        y2 += ['best']

        # remove entirely masked rows and columns
        yidx = where(~D.mask.all(axis = 0))[0]
        D    = D[:, yidx]
        xidx = where(~D.mask.all(axis = 1))[0]
        D    = D[xidx, :]
        x2   = [x2[i] for i in xidx]
        y2   = [y2[i] for i in yidx]

        fig, ax = plt.subplots()
        plt.pcolor(D.T, vmin = zmin, vmax = zmax)
        plt.axis([0, len(x2), 0, len(y2)])

        ax.set_xticks(arange(len(x2)) + 0.5, minor = False)
        ax.set_yticks(arange(len(y2)) + 0.5, minor = False)
        ax.set_xticks(arange(len(x2)), minor = True)
        ax.set_yticks(arange(len(y2)), minor = True)
        ax.invert_yaxis()
        if self.anon:
            x_labels = ['%s %d' % (x_label, i) for i in range(1, len(x2))] + ['best']
            ax.set_xticklabels(x_labels, minor = False)
        else:
            ax.set_xticklabels(x2, minor = False)
        ax.set_yticklabels(y2, minor = False)

        plt.xticks(rotation = 90)
        plt.grid(which = 'minor', linestyle = '-')
        plt.tick_params(axis = 'both', which = 'major', bottom = 'off', top = 'off')
        plt.colorbar()
        plt.title(title, fontsize = 10)
        plt.savefig(filename, bbox_inches = 'tight')
        plt.close()

    def __condition(self, arr, dims, cross):
        """ arr is a numpy array, dims is an ordered dictionary of lists,
            cross is a dictionary of (variable, values)
        """
        if not cross: return arr.copy(), dims.copy()

        var, vals = cross.items()[0]

        vars  = var.split('/')
        nvars = len(vars)

        # variable index
        varidx = [dims.keys().index(v) for v in vars]

        # expand '*'
        vals2 = vals[:]
        for i in range(nvars):
            if vals[i] == '*':
                vals2[i] = dims[vars[i]]

        # output dimensions
        dims_out = dims.copy()
        for i in range(nvars):
            dims_out.pop(vars[i]);
        dims_out[var] = ['/'.join(v) for v in list(product(*vals2))] # cross-product

        # dimension indices
        dimidx = [0] * nvars
        for i in range(nvars):
            v = vals2[i][:]
            dimidx[i] = [0] * len(v)
            for j in range(len(v)):
                if v[j].isdigit():
                    v[j] = double(v[j]) # convert to number
                dimidx[i][j] = dims[vars[i]].index(v[j])
        dimidx = list(product(*dimidx))

        # output array
        sh = [len(d) for d in dims_out.itervalues()]
        arr_out = masked_array(zeros(sh), mask = ones(sh))

        slice_in  = [slice(0, size) for size in arr.shape]
        slice_out = [slice(0, size) for size in arr_out.shape]

        for i in range(len(dimidx)):
            idx = dimidx[i]
            for j in range(nvars):
                slice_in[varidx[j]] = idx[j]
            slice_out[-1] = i # replace last element
            arr_out[slice_out] = arr[slice_in]

        return arr_out, dims_out

    def __average(self, arr, dims, cross):
        """ arr is a numpy array, dims is an ordered dictionary of lists,
            cross is a dictionary of (variable, values)
        """
        arr_out, dims_out = self.__condition(arr, dims, cross)

        if not cross: return arr.copy(), dims.copy()

        var    = cross.keys()[0]
        varidx = dims_out.keys().index(var)

        wts = self.__loadwts(var, dims_out[var])

        sh = arr_out.shape
        slice_arr = [slice(0, size) for size in arr_out.shape]

        wts_full = masked_array(zeros(sh), mask = ones(sh))
        for i in range(len(wts)):
            slice_arr[varidx] = i
            wts_full[slice_arr] = wts[i]
        wts_full = masked_where(arr_out.mask, wts_full) # mask

        arr_out = (wts_full * arr_out).sum(axis = varidx) / wts_full.sum(axis = varidx) # average over variable

        arr_out = resize(arr_out, arr_out.shape + (1,))
        dims_out[var] = ['ave']

        return arr_out, dims_out

    def __optimize(self, arr, dims, cross):
        """ arr is a numpy array, dims is an ordered dictionary of lists,
            cross is a dictionary of (variable, values)
        """
        arr_out, dims_out = self.__condition(arr, dims, cross)

        if not cross: return arr.copy(), dims.copy()

        var    = cross.keys()[0]
        varidx = dims_out.keys().index(var)

        if self.metric in ['tscorr', 'hitrate']: # TODO: make more generic
            arr_out = arr_out.max(axis = varidx) # maximize over variable
        else:
            arr_out = arr_out.min(axis = varidx) # minimize over variable

        arr_out = resize(arr_out, arr_out.shape + (1,))
        dims_out[var] = ['opt']

        return arr_out, dims_out

    def __loadwts(self, var, vals):
        if self.wfile is None:
            return ones(len(vals))
        else:
            vars  = var.split('/')
            nvars = len(vars)

            v = [0] * nvars
            w = [0] * nvars
            with nc(self.wfile) as f:
                for i in range(nvars):
                    if f.variables[vars[i]].units == 'mapping':
                        v[i] = array(f.variables[vars[i]].long_name.split(', '))
                    else:
                        v[i] = f.variables[vars[i]][:]
                    w[i] = f.variables['weights_' + vars[i]][:]

            nvals = len(vals)
            wts = masked_array(zeros(nvals), mask = ones(nvals))
            for i in range(nvals):
                svals = vals[i].split('/')
                for j in range(nvars):
                    if svals[j].isdigit():
                        svals[j] = double(svals[j]) # convert to number
                    idx = where(v[j] == svals[j])[0]
                    if idx.size:
                        if isMaskedArray(wts[i]):
                            wts[i] = w[j][idx[0]]
                        else:
                            wts[i] *= w[j][idx[0]]

            return wts