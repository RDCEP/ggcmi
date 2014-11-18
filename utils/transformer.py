import abc
from numpy import sqrt
from scipy.stats import scoreatpercentile
from numpy.ma import median, isMaskedArray

class Transformer(object):
    @abc.abstractmethod
    def transform(self, pfull, p, q): return

class NoTransformer(Transformer):
    def transform(self, pfull, p, q): return pfull.copy() # no transformation

class VSTransformer(Transformer):
    def transform(self, pfull, p, q):
        p2 = pfull.copy()

        p2 -= p2.mean() # decenter

        pv = p.var()
        if pv: p2 *= sqrt(q.var() / pv) # scale

        p2 += q.mean() - p2.mean() # shift

        return p2

class MSTransformer(Transformer):
    def transform(self, pfull, p, q):
        p2 = pfull.copy()

        pm = p.mean()

        return p2 * q.mean() / pm if not isMaskedArray(pm) and pm else p2

class QMTransformer(Transformer):
    def transform(self, pfull, p, q):
        pfull2 = pfull.copy() # copy
        p2     = p.copy()
        q2     = q.copy()

        q2mean = q2.mean() # target mean

        pfull2 -= median(pfull2) # center
        p2     -= median(p2)
        q2     -= median(q2)

        q5  = self.__percentile(q2, 5)
        q25 = self.__percentile(q2, 25)
        q75 = self.__percentile(q2, 75)
        q95 = self.__percentile(q2, 95)

        p25 = self.__percentile(p2, 25) # 25th
        if p25: pfull2[pfull2 < 0] *= q25 / p25

        p75 = self.__percentile(p2, 75) # 75th
        if p75: pfull2[pfull2 > 0] *= q75 / p75

        p5 = self.__percentile(p2, 5) # 5th
        if p5: pfull2[pfull2 < q25] *= q5 / p5

        p95 = self.__percentile(p2, 95) # 95th
        if p95: pfull2[pfull2 > q75] *= q95 / p95

        pfull2 += q2mean - pfull2.mean() # reshift

        return pfull2

    def __percentile(self, p, per):
        p2 = p[~p.mask] if isMaskedArray(p) else p.copy()
        return scoreatpercentile(p2, per)

class TransformerWrapper(object):
    def __init__(self, method):
        if method == 'none':
            self.transformer = NoTransformer()
        elif method == 'variance-scale':
            self.transformer = VSTransformer()
        elif method == 'mean-scale':
            self.transformer = MSTransformer()
        elif method == 'quantile-mapping':
            self.transformer = QMTransformer()
        else:
            raise Exception('Unrecognized transform method')

    def transform(self, pfull, p, q): return self.transformer.transform(pfull, p, q)