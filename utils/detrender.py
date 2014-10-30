import abc
from numpy import zeros, arange, polyfit, isnan, nan, diff
from numpy.ma import masked_array, isMaskedArray, masked_where

class Detrender(object):
    @abc.abstractmethod
    def detrend(self, y): return

    def polytrend(self, y, order):
        ly = len(y)
        x = arange(1, 1 + ly)
        if isMaskedArray(y):
            mask = y.mask
            line = masked_array(zeros(ly), mask = mask)
            if not mask.all():
                x2 = x[~mask]
                y2 = y[~mask]
                if len(x2) < order + 1: # fewer points than order of polynomial
                    line[~mask] = y2
                else:
                    coef = polyfit(x2, y2, order)
                    for i in range(order + 1):
                        line[~mask] += coef[i] * x2 ** (order - i)
        else:
            coef = polyfit(x, y, order)
            line = zeros(len(y))
            for i in range(order + 1):
                line += coef[i] * x ** (order - i)
        return line

class NoDetrender(Detrender):
    def detrend(self, y): return y, zeros(len(y)) # no trend line

class PolyDetrender(Detrender):
    def __init__(self, order): self.order = order

    def detrend(self, y):
        line = self.polytrend(y, self.order)
        dy = y - line
        return dy, line

class MovingAveDetrender(Detrender):
    def detrend(self, y):
        dy = y.copy(); dy[:] = nan
        line = y.copy(); line[:] = nan
        line[2 : -2] = self.__winave(y, 5, 7)
        dy[2 : -2] = y[2 : -2] - line[2 : -2]
        dy = masked_where(isnan(dy), dy)
        line = masked_where(isnan(line), line)
        return dy, line

    def __winave(self, d, nstart, nend):
        # perform n-point window average of data d
        nd = len(d)
        dave = d.copy()
        numramp = (nend - nstart) / 2
        n = nstart - 2
        for i in range(nd):
            if i <= numramp:
                n += 2
            elif i >= nd - numramp:
                n -= 2
            idx0 = max(i - n / 2, 0)
            idx1 = min(i + n / 2 + 1, nd)
            dave[i] = d[idx0 : idx1].mean()
        dave = dave[n / 2 : -(n / 2)] # remove left and right edges
        return dave

class FFDDetrender(Detrender):
    def detrend(self, y):
        dy = y.copy(); dy[:] = nan
        dy[1 :] = diff(y) / y[: -1]
        dy = masked_where(isnan(dy), dy)
        line = zeros(len(y)) # no trend line
        return dy, line

class FFDRemoveDetrender(Detrender):
    def detrend(self, y):
        dy = y.copy(); dy[:] = nan
        line = y.copy(); line[:] = nan
        ffd = diff(y) / y[: -1]
        line[1 :] = self.polytrend(ffd, 1)
        dy[1 :] = ffd - line[1 :]
        dy = masked_where(isnan(dy), dy)
        line = masked_where(isnan(line), line)
        return dy, line

class DetrenderWrapper(Detrender):
    def __init__(self, dt, mp):
        if dt == 'none':
            self.detrender = NoDetrender()
        elif dt == 'lin':
            self.detrender = PolyDetrender(1)
        elif dt == 'quad':
            self.detrender = PolyDetrender(2)
        elif dt == 'ma':
            self.detrender = MovingAveDetrender()
        elif dt == 'ffd':
            self.detrender = FFDDetrender()
        elif dt == 'ffdtr':
            self.detrender = FFDRemoveDetrender()
        else:
            raise Exception('Unrecognized detrend method')
        if not mp in ['true', 'false']:
            raise Exception('Unrecognized mean-preserving method')
        self.mp = mp

    def detrend(self, y):
        dy, line = self.detrender.detrend(y)
        if self.mp == 'true':
            return dy - dy.mean() + y.mean(), line # mean-preserving
        else:
            return dy - dy.mean(), line # decenter