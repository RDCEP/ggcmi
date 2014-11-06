from re import findall
from itertools import product
from netCDF4 import Dataset as nc
from numpy.ma import masked_array, masked_where, argmin
from numpy import inf, zeros, ones, logical_and, arange, resize

def makeMasked(sh): return masked_array(zeros(sh), mask = ones(sh))

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

        self.tmin = inf
        self.tmax = -inf
        for i in range(self.nm):
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
        for a in range(naggs):
            for d, m, c in product(range(ndt), range(nmp), range(ncr)):
                met = self.metrics[:, a, d, m, c]
                ydt = self.yield_detr[:, a, :, d, m, c].T
                yrt = self.yield_retr[:, a, :, d, m, c].T

                yield_detr_mean[a, :, d, m, c] = self.__ensemble(met, ydt)
                yield_retr_mean[a, :, d, m, c] = self.__ensemble(met, yrt)

        return yield_detr_mean, yield_retr_mean

    def __ensemble(self, met, arr):
        """ met: array of length nmodels, arr: matrix of size (ntime, nmodels) """
        ntime, nmodels = arr.shape

        sh = (ntime, nmodels, 2)
        mu = makeMasked(sh)

        if met.mask.all(): return mu # all metrics are masked

        arr2 = arr[:, ~met.mask]
        met2 = met[~met.mask]

        nmodels2 = len(met2)

        sortidx = [i[0] for i in sorted(enumerate(self.__cost(met2)), key = lambda x: x[1])]
        metsort = met2[sortidx]
        arrsort = arr2[:, sortidx]
        metsort = resize(metsort, (ntime, nmodels2))

        mu[:, : nmodels2, 0] = arrsort.cumsum(axis = 1) / ones((ntime, nmodels2)).cumsum(axis = 1)
        mu[:, : nmodels2, 1] = (arrsort * metsort).cumsum(axis = 1) / metsort.cumsum(axis = 1)

        return mu

    def __cost(self, met): # metric "cost"
        if self.metricname in ['tscorr', 'hitrate']:
            return -met
        elif self.metricname == 'varratio':
            return abs(met - 1)
        else:
            return met