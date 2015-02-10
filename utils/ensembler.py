from re import findall
from os.path import basename
from itertools import product
from netCDF4 import Dataset as nc
from numpy.ma import masked_array, masked_where, argmin
from numpy import inf, zeros, ones, logical_and, arange, resize, where

def makeMasked(sh, dtype = float): return masked_array(zeros(sh), mask = ones(sh), dtype = dtype)

class Ensembler(object):
    def __init__(self, bcfiles, mmfiles, agglvl, metricname):
        self.nm = len(bcfiles)
        self.metricname = metricname

        with nc(bcfiles[0]) as f:
            self.aggs        = f.variables[agglvl][:]
            self.aggunits    = f.variables[agglvl].units
            self.agglongname = f.variables[agglvl].long_name
            self.dt          = f.variables['dt'].long_name.split(', ')
            self.mp          = f.variables['mp'].long_name.split(', ')
            self.cr          = f.variables['cr'].long_name.split(', ')
        naggs, ndt, nmp, ncr = len(self.aggs), len(self.dt), len(self.mp), len(self.cr)

        with nc(mmfiles[0]) as f:
            self.metricunits = f.variables[metricname].units

        self.models = [0] * self.nm
        self.tmin   = inf
        self.tmax   = -inf
        for i in range(self.nm):
            self.models[i] = basename(bcfiles[i]).split('_')[0]
            with nc(bcfiles[i]) as f:
                time        = f.variables['time'][:]
                time_units  = f.variables['time'].units
                time       += int(findall(r'\d+', time_units)[0])
                self.tmin = min(self.tmin, time[0])
                self.tmax = max(self.tmax, time[-1])
        self.time = arange(self.tmin, self.tmax + 1)
        ntime     = len(self.time)

        self.metrics    = makeMasked((self.nm, naggs, ndt, nmp, ncr))
        self.yield_detr = makeMasked((self.nm, naggs, ntime, ndt, nmp, ncr))
        self.yield_retr = makeMasked((self.nm, naggs, ntime, ndt, nmp, ncr))
        for i in range(self.nm):
            with nc(mmfiles[i]) as f:
                metric  = f.variables[metricname][:, :, :, :, :, 0]
                metric  = masked_where(metric < 0, metric) # mask negative values
                scenidx = argmin(self.__cost(metric), axis = 1)

            with nc(bcfiles[i]) as f:
                time        = f.variables['time'][:]
                time_units  = f.variables['time'].units
                time       += int(findall(r'\d+', time_units)[0])
                ydt         = f.variables['yield_detrend'][:]
                yrt         = f.variables['yield_retrend'][:]

                inrange = logical_and(self.time >= time[0], self.time <= time[-1])

                for a in range(naggs):
                    for d, m, c in product(range(ndt), range(nmp), range(ncr)):
                        si = scenidx[a, d, m, c]
                        self.metrics[i, a, d, m, c]             = metric[a, si, d, m, c]
                        self.yield_detr[i, a, inrange, d, m, c] = ydt[a, :, si, d, m, c]
                        self.yield_retr[i, a, inrange, d, m, c] = yrt[a, :, si, d, m, c]

    def average(self):
        nmodels, naggs, ntime, ndt, nmp, ncr = self.yield_detr.shape

        yield_detr_mean = makeMasked((naggs, ntime, ndt, nmp, ncr, nmodels, 2))
        yield_retr_mean = makeMasked((naggs, ntime, ndt, nmp, ncr, nmodels, 2))

        model_order   = makeMasked((naggs, ndt, nmp, ncr, nmodels))
        model_weights = makeMasked((naggs, ndt, nmp, ncr, nmodels))

        for a in range(naggs):
            for d, m, c in product(range(ndt), range(nmp), range(ncr)):
                met = self.metrics[:, a, d, m, c]
                ydt = self.yield_detr[:, a, :, d, m, c].T
                yrt = self.yield_retr[:, a, :, d, m, c].T

                # order models
                order, weights = self.__order_models(met)
                model_order[a, d, m, c]   = order + 1
                model_weights[a, d, m, c] = weights

                # compute ensemble
                yield_detr_mean[a, :, d, m, c] = self.__ensemble(met, ydt, order)
                yield_retr_mean[a, :, d, m, c] = self.__ensemble(met, yrt, order)

        return yield_detr_mean, yield_retr_mean, model_order, model_weights

    def __order_models(self, met):
        nmodels = len(met)

        order   = makeMasked((nmodels), dtype = int)
        weights = makeMasked((nmodels))

        if met.mask.all(): return order, weights

        midx = where(~met.mask)[0]
        metunmasked = met[~met.mask]

        sortidx = [i[0] for i in sorted(enumerate(self.__cost(metunmasked)), key = lambda x: x[1])]
        order[: len(sortidx)] = midx[sortidx]
        weights[: len(sortidx)] = metunmasked[sortidx]

        return order, weights

    def __ensemble(self, met, arr, order):
        """ met:   array of length nmodels of metrics,
            arr:   matrix of size (ntime, nmodels) of data,
            order: array of length nmodels indicating model order """
        ntime, nmodels = arr.shape

        sh = (ntime, nmodels, 2)
        mu = makeMasked(sh)

        if met.mask.all(): return mu # all metrics are masked

        order    = order[~order.mask]
        nmodels2 = len(order)

        arrsort = arr[:, order]
        metsort = resize(met[order], (ntime, nmodels2))

        metsort[arrsort.mask] = 0 # zero out masked
        arrsort[arrsort.mask] = 0

        wts = metsort.cumsum(axis = 1)
        wts = masked_where(wts == 0, wts) # mask zeros

        mu[:, : nmodels2, 0] = arrsort.cumsum(axis = 1) / ones((ntime, nmodels2)).cumsum(axis = 1)
        mu[:, : nmodels2, 1] = (arrsort * metsort).cumsum(axis = 1) / wts

        return mu

    def __cost(self, met): # metric "cost"
        if self.metricname in ['tscorr', 'hitrate']:
            return -met
        elif self.metricname == 'varratio':
            return abs(met - 1)
        else:
            return met