import abc

class Retrender(object):
    @abc.abstractmethod
    def retrend(self, T, R, o): return

class SumRetrender(Retrender):
    def retrend(self, T, R, o): return T + R

class ScaleRetrender(Retrender):
    def retrend(self, T, R, o): return T * R / R.mean()

class FFDLeftRetrender(Retrender):
    def retrend(self, T, R, o):
        r = T.copy()
        r[0] = o[0] # initial condition = leftmost point
        for i in range(1, len(r)):
            r[i] = r[i - 1] * (1. + T[i] + R[i])
        return r

class FFDRightRetrender(Retrender):
    def retrend(self, T, R, o):
        r = T.copy()
        r[-1] = o[-1] # initial condition = rightmost point
        for i in range(-1, -len(r), -1):
            if T.mask[i] or R.mask[i]:
                r.mask[i - 1] = True # prevent division by masked
            else:
                r[i - 1] = r[i] / (1. + T[i] + R[i])
        return r

class RetrenderWrapper(Retrender):
    def __init__(self, dt, mp):
        if dt == 'ffd':
            self.retrender = FFDLeftRetrender()
        elif mp == 'true':
            self.retrender = ScaleRetrender()
        else:
            self.retrender = SumRetrender()

    def retrend(self, T, R, o): return self.retrender.retrend(T, R, o)