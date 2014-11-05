from numpy.ma import masked_array
from retrender import RetrenderWrapper
from detrender import DetrenderWrapper
from transformer import TransformerWrapper
from numpy import zeros, ones, logical_and

class BiasCorrecter(object):
    def __init__(self, dt, mp, cr):
        self.transformer = TransformerWrapper(cr)
        self.retrender   = RetrenderWrapper(dt, mp)
        self.detrender   = DetrenderWrapper(dt, mp)

        self.unrealizable = mp == 'false' and cr == 'MS'

    def correct(self, sim, obs, tsim, tobs):
        """ Corrects simulated sequence sim using observations obs,
            with tsim and tobs arrays of years """
        tmin = max([tsim[0],  tobs[0]])
        tmax = min([tsim[-1], tobs[-1]])

        odT, T = self.detrender.detrend(obs) # detrend observations
        sdT, _ = self.detrender.detrend(sim) # detrend simulations

        dT = masked_array(zeros(len(tsim)), mask = ones(len(tsim)))
        rT = dT.copy()
        if not self.unrealizable and not sdT.mask.all() and not odT.mask.all():
            odT_c = self.__toverlap(odT, tobs, tmin, tmax) # transform distribution
            sdT_c = self.__toverlap(sdT, tsim, tmin, tmax)
            dT    = self.transformer.transform(sdT, sdT_c, odT_c)

            obs_c = self.__toverlap(obs, tobs, tmin, tmax) # combine
            T_c   = self.__toverlap(T,   tobs, tmin, tmax)
            dT_c  = self.__toverlap(dT,  tsim, tmin, tmax)
            rT_c  = self.retrender.retrend(T_c, dT_c, obs_c)
            rT[logical_and(tsim >= tmin, tsim <= tmax)] = rT_c

        return dT, rT

    def __toverlap(self, x, tx, tmin, tmax):
        return x[logical_and(tx >= tmin, tx <= tmax)]