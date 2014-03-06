import abc, numpy.ma as ma
from numpy import zeros, ones, where, resize, cos, pi

class Averager(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def av(self, var, agg, lats, weights = None): return

    def combine(self, var1, var2, agg, lats, weights1 = None, weights2 = None):
        nt, nlats, nlons = var1.shape
        if weights1 is None: weights1 = ones((nlats, nlons)) # weights
        if weights2 is None: weights2 = ones((nlats, nlons))
        area = self.area(lats, nlats, nlons) # area
        av1 = self.av(var1, agg, lats, weights = weights1)
        av2 = self.av(var2, agg, lats, weights = weights2)
        area1 = self.areas(agg, area, weights1)
        area2 = self.areas(agg, area, weights2)
        sz = len(av1)
        area1 = resize(area1, (nt, sz)).T
        area2 = resize(area2, (nt, sz)).T
        totarea = area1 + area2
        totarea = ma.masked_where(totarea == 0, totarea)
        totav = ma.masked_array(zeros((sz, nt, 3)), mask = ones((sz, nt, 3)))
        totav[:, :, 0] = av1
        totav[:, :, 1] = av2
        totav[:, :, 2] = (area1 * av1 + area2 * av2) / totarea
        return totav

    def sum(self, var, agg, area, weights):
        nt, nlats, nlons = var.shape
        aggvals = self.__uniquevals(agg)
        sz = len(aggvals)
        aselect = ma.zeros((sz, nlats, nlons), dtype = bool)
        for i in range(sz): aselect[i] = (agg == aggvals[i])
        ridx, latidx, lonidx = where(aselect)
        sumv = zeros((sz, nt))
        vartmp = zeros((sz, nlats, nlons))
        for t in range(nt):
            vartmp[ridx, latidx, lonidx] = var[t, latidx, lonidx]  * \
                                           weights[latidx, lonidx] * \
                                           area[latidx, lonidx]    * \
                                           aselect[ridx, latidx, lonidx]
            sumv[:, t] = vartmp.sum(axis = 2).sum(axis = 1)
        return sumv

    def areas(self, agg, area, weights):
        aggvals = self.__uniquevals(agg)
        areas = zeros(len(aggvals))
        for i in range(len(aggvals)):
            areas[i] = (weights * area * (agg == aggvals[i])).sum()
        return areas

    def area(self, lats, nlats, nlons):
        A = 100 * (111.2 / 2) ** 2 * cos(pi * lats / 360)
        A = resize(A, (nlons, nlats)).T
        return A

    def __uniquevals(self, d):
        u = ma.unique(d)
        u = u[~u.mask]
        return u

class SumAverager(Averager):    
    def av(self, var, agg, lats, weights = None):
        _, nlats, nlons = var.shape
        if weights is None: weights = ones((nlats, nlons)) # weights
        area = self.area(lats, nlats, nlons) # area
        avv = self.sum(var, agg, area, weights)
        return avv

class MeanAverager(Averager):
    def av(self, var, agg, lats, weights = None):
        nt, nlats, nlons = var.shape
        if weights is None: weights = ones((nlats, nlons)) # weights
        area = self.area(lats, nlats, nlons) # area
        avv = self.sum(var, agg, area, weights)
        areas = self.areas(agg, area, weights)
        valididx = areas > 0.
        if valididx.sum():
            for t in range(nt):
                avv[valididx, t] /= areas[valididx]
        return avv