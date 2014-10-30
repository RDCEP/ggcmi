import abc
from numpy import sqrt
from scipy.stats import scoreatpercentile
from numpy.ma import median, isMaskedArray

class Transformer(object):
    @abc.abstractmethod
    def transform(self, p, q): return

class NoTransformer(Transformer):
    def transform(self, p, q): return p # no transformation

class VSTransformer(Transformer):
    def transform(self, p, q):
        p -= p.mean() # decenter
        pv = p.var()
        if pv: p *= sqrt(q.var() / pv) # scale
        p += q.mean() - p.mean() # shift
        return p

class MSTransformer(Transformer):
    def transform(self, p, q):
        pm = p.mean()
        return p * q.mean() / pm if not isMaskedArray(pm) and pm else p

class QMTransformer(Transformer):
    def transform(self, p, q):
        p2 = p.copy() # copy
        q2 = q.copy()

        q2mean = q2.mean() # target mean

        p2 -= median(p2) # center
        q2 -= median(q2)

        q5 = self.__percentile(q2, 5)
        q25 = self.__percentile(q2, 25)
        q75 = self.__percentile(q2, 75)
        q95 = self.__percentile(q2, 95)

        p25 = self.__percentile(p2, 25) # 25th
        if p25: p2[p2 < 0] *= q25 / p25

        p75 = self.__percentile(p2, 75) # 75th
        if p75: p2[p2 > 0] *= q75 / p75

        p5 = self.__percentile(p2, 5) # 5th
        if p5: p2[p2 < q25] *= q5 / p5

        p95 = self.__percentile(p2, 95) # 95th
        if p95: p2[p2 > q75] *= q95 / p95

        p2 += q2mean - p2.mean() # reshift

        return p2

    def __percentile(self, p, per):
        p2 = p[~p.mask] if isMaskedArray(p) else p.copy()
        return scoreatpercentile(p2, per)

class TransformerWrapper(object):
    def __init__(self, method):
        if method == 'none':
            self.transformer = NoTransformer()
        elif method == 'VS':
            self.transformer = VSTransformer()
        elif method == 'MS':
            self.transformer = MSTransformer()
        elif method == 'QM':
            self.transformer = QMTransformer()
        else:
            raise Exception('Unrecognized transform method')

    def transform(self, p, q): return self.transformer.transform(p, q)