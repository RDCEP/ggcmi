from retrender import RetrenderWrapper
from detrender import DetrenderWrapper
from numpy.ma import masked_array, where
from transformer import TransformerWrapper
from numpy import zeros, ones, intersect1d

class BiasCorrecter(object):
    def __init__(self, dt, mp, cr):
        self.transformer = TransformerWrapper(cr)
        self.retrender   = RetrenderWrapper(dt, mp)
        self.detrender   = DetrenderWrapper(dt, mp)

        self.unrealizable = mp == 'false' and cr == 'MS'

    def correct(self, sim, obs, tsim, tobs):
        """ Corrects simulated sequence sim using observations obs,
            with tsim and tobs arrays of years """
        time = intersect1d(tsim, tobs) # get common times
        if not len(time): raise Exception('No common times')

        odT, T = self.detrender.detrend(obs) # detrend observations
        odT, T = self.__toverlap(odT, tobs, time), self.__toverlap(T, tobs, time)

        sdT = self.detrender.detrend(sim)[0] # detrend simulations
        sdT = self.__toverlap(sdT, tsim, time)

        yhat = masked_array(zeros(len(time)), mask = ones(len(time)))
        R    = yhat.copy()
        if not self.unrealizable and not sdT.mask.all() and not odT.mask.all():
            R2 = self.transformer.transform(sdT, odT)   # transform distribution
            obs_trim = self.__toverlap(obs, tobs, time) # combine
            yhat2 = self.retrender.retrend(T, R2, obs_trim)
            ln = min(len(time), len(yhat2))
            yhat[: ln], R[: ln] = yhat2[: ln], R2[: ln]

        return yhat, R

    def __toverlap(self, x, tx, tref):
        ntime = len(tref)
        xref = masked_array(zeros(ntime), mask = ones(ntime))
        for t in range(ntime):
            tidx = where(tx == tref[t])[0][0]
            xref[t] = x[tidx]
        return xref