import abc
from numpy import sqrt, zeros, ones, where, double
from numpy.ma import corrcoef, masked, intersect1d, masked_array

class Metrics(object):
    def eval(self, sim, obs, obsbase, time):
        if not sim.mask.all() and not obs.mask.all():
            return self.metric(sim, obs, time, obsbase)
        else:
            return masked

    @abc.abstractmethod
    def metric(self, sim, obs, time, obsbase = None): return

    def tslice(self, x, t, tref):
        tcommon = intersect1d(t, tref)
        nt = len(tcommon)

        xs = masked_array(zeros(nt), mask = ones(nt))
        for i in range(nt):
            tidx = where(t == tcommon[i])[0][0]
            xs[i] = x[tidx]

        return xs

class TSCorr(Metrics):
    def metric(self, sim, obs, time, obsbase = None): return corrcoef(sim, obs)[0, 1]

class VarRatio(Metrics):
    def metric(self, sim, obs, time, obsbase = None):
        sim_mean = sim.mean()
        obs_std  = obs.std()
        if sim_mean * obs_std: # no zeros
            return (sim.std() * obs.mean()) / (sim_mean * obs_std)
        else:
            return masked

class RMSE(Metrics):
    def metric(self, sim, obs, time, obsbase = None):
        rmse = sqrt(((sim - obs) ** 2).mean())
        if obsbase is not None and not obsbase.mask.all():
            norm = obsbase.mean()
            if norm:
                rmse /= norm
        return rmse

class HitRate(Metrics):
    def metric(self, sim, obs, time, obsbase = None):
        st = time[sim < sim.mean() - sim.std()]
        ot = time[obs < obs.mean() - obs.std()]
        tc = intersect1d(ot, st)
        return double(tc.size) / ot.size if ot.size else masked

class RMSEExtreme(Metrics):
    def metric(self, sim, obs, time, obsbase = None):
        st = time[sim < sim.mean() - sim.std()]
        ot = time[obs < obs.mean() - obs.std()]
        tc = intersect1d(ot, st)

        sim_extreme = self.tslice(sim, time, tc)
        obs_extreme = self.tslice(obs, time, tc)

        if sim_extreme.size:
            return sqrt(((obs_extreme - sim_extreme) ** 2).mean())
        else:
            return masked

class BiasExtreme(Metrics):
    def metric(self, sim, obs, time, obsbase = None):
        st = time[sim < sim.mean() - sim.std()]
        ot = time[obs < obs.mean() - obs.std()]
        tc = intersect1d(ot, st)

        sim_extreme = self.tslice(sim, time, tc)
        obs_extreme = self.tslice(obs, time, tc)

        if sim_extreme.size:
            return (obs_extreme - sim_extreme).mean()
        else:
            return masked

class MetricsWrapper(Metrics):
    def __init__(self, method):
        super(MetricsWrapper, self).__init__()

        if method == 'tscorr':
            self.metrics = TSCorr()
        elif method == 'varratio':
            self.metrics = VarRatio()
        elif method == 'rmse':
            self.metrics = RMSE()
        elif method == 'hitrate':
            self.metrics = HitRate()
        elif method == 'rmse_extreme':
            self.metrics = RMSEExtreme()
        elif method == 'bias_extreme':
            self.metrics = BiasExtreme()
        else:
            raise Exception('Unrecognized metric')

    def metric(self, sim, obs, time, obsbase = None): return self.metrics.metric(sim, obs, time, obsbase = obsbase)