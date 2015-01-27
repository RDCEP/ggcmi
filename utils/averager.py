import abc, numpy.ma as ma
from numpy import zeros, ones, where, resize, cos, pi, logical_not, logical_and

class Averager(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def av(self, var, agg, lats, weights = None, calcarea = False, mask = None, numchunks = 1): return

    def combine(self, var1, var2, agg, lats, weights1 = None, weights2 = None, calcarea = True, mask1 = None, mask2 = None, numchunks = 1):
        nt, nlats, nlons = var1.shape

        av1 = self.av(var1, agg, lats, weights1, calcarea, mask1, numchunks)
        av2 = self.av(var2, agg, lats, weights2, calcarea, mask2, numchunks)

        area1 = self.areas(var1, agg, lats, weights1, calcarea, mask1)
        area2 = self.areas(var2, agg, lats, weights2, calcarea, mask2)

        av1[logical_and(av1.mask, ~av2.mask)] = 0. # zero out points masked in one but unmasked in other
        av2[logical_and(av2.mask, ~av1.mask)] = 0.
        area1[logical_and(area1.mask, ~area2.mask)] = 0.
        area2[logical_and(area2.mask, ~area1.mask)] = 0.

        totarea = area1 + area2
        totarea = ma.masked_where(totarea == 0, totarea)

        sz = len(av1)
        totav = ma.masked_array(zeros((sz, nt, 3)), mask = ones((sz, nt, 3)))
        totav[:, :, 0], totav[:, :, 1] = av1, av2
        totav[:, :, 2] = (area1 * av1 + area2 * av2) / totarea

        return totav

    def sum(self, var, agg, lats, weights = None, calcarea = False, mask = None, numchunks = 1):
        nt, nlats, nlons = var.shape

        if weights is None: # weights
            weights = ones((nt, nlats, nlons))
        elif len(weights.shape) == 2:
            weights = resize(weights, (nt, nlats, nlons))

        if calcarea: # area
            area = self.area(lats, nlats, nlons)
        else:
            area = ones((nlats, nlons))

        aggvals = self.__uniquevals(agg)
        sz = len(aggvals)

        if mask is None:
            varmask = ones((nt, nlats, nlons)) # no additional mask
        else:
            varmask = mask

        chunksize = sz / numchunks # chunk data to reduce memory usage

        sumv = ma.masked_array(zeros((sz, nt)), mask = ones((sz, nt)))

        maxchunksize = max(chunksize, chunksize + sz - chunksize * numchunks)

        aselect = ma.zeros((maxchunksize, nlats, nlons), dtype = bool) # preallocate
        vartmp  = ma.zeros((maxchunksize, nlats, nlons))

        cnt = 0
        for i in range(numchunks):
            startidx = cnt
            if i != numchunks - 1:
                endidx = cnt + chunksize
            else:
                endidx = sz

            aggvalsc = aggvals[startidx : endidx] # work on subset of aggregation values
            szc = len(aggvalsc)

            aselect[:] = 0 # clear
            for j in range(szc): aselect[j] = (agg == aggvalsc[j])
            ridx, latidx, lonidx = where(aselect)

            vartmp[:] = 0 # clear
            vartmp.mask = ones(vartmp.shape)
            for t in range(nt):
                vartmp[ridx, latidx, lonidx] = var[t, latidx, lonidx]        * \
                                               varmask[t, latidx, lonidx]    * \
                                               weights[t, latidx, lonidx]    * \
                                               area[latidx, lonidx]          * \
                                               aselect[ridx, latidx, lonidx]
                sumv[startidx : endidx, t] = vartmp.sum(axis = 2).sum(axis = 1)[: szc]

            cnt += chunksize

        return sumv

    def areas(self, var, agg, lats, weights = None, calcarea = False, mask = None):
        nt, nlats, nlons = var.shape

        if weights is None: # weights
            weights = ones((nt, nlats, nlons))
        elif len(weights.shape) == 2:
            weights = resize(weights, (nt, nlats, nlons))

        if calcarea: # area
            area = self.area(lats, nlats, nlons)
        else:
            area = ones((nlats, nlons))

        aggvals = self.__uniquevals(agg)
        sz = len(aggvals)

        varmask = logical_not(var.mask) if ma.isMaskedArray(var) else ones(var.shape) # use variable mask
        if not mask is None: varmask = logical_and(varmask, mask) # additional mask

        areas  = ma.masked_array(zeros((sz, nt)), mask = ones((sz, nt)))
        vartmp = zeros((nt, nlats, nlons))
        for i in range(len(aggvals)):
            warea = weights * area * (agg == aggvals[i])
            tidx, latidx, lonidx = ma.where(warea)

            vartmp[:] = 0
            vartmp[tidx, latidx, lonidx] = warea[tidx, latidx, lonidx] * \
                                           varmask[tidx, latidx, lonidx]
            areas[i] = vartmp.sum(axis = 2).sum(axis = 1)

        areas = ma.masked_where(areas == 0, areas)
        areas.mask = resize(areas.mask, areas.shape) # ensure mask is same size as data

        return areas

    def area(self, lats, nlats, nlons):
        A = 100 * (111.2 / 2) ** 2 * cos(pi * lats / 180)
        A = resize(A, (nlons, nlats)).T
        return A

    def __uniquevals(self, d):
        u = ma.unique(d)
        u = u[~u.mask]
        return u

class SumAverager(Averager):    
    def av(self, var, agg, lats, weights = None, calcarea = False, mask = None, numchunks = 1):
        return self.sum(var, agg, lats, weights, calcarea, mask, numchunks)

class MeanAverager(Averager):
    def av(self, var, agg, lats, weights = None, calcarea = False, mask = None, numchunks = 1):
        avv   = self.sum(var, agg, lats, weights, calcarea, mask, numchunks)
        areas = self.areas(var, agg, lats, weights, calcarea, mask)
        return avv / areas