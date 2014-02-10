#!/usr/bin/env python

# import modules
import sys, re, abc
import numpy.ma as ma
from os import listdir
from netCDF4 import Dataset as nc
from optparse import OptionParser
from os.path import split, splitext, sep
from scipy.stats import scoreatpercentile
from numpy import ceil, double, intersect1d, where, sqrt, polyfit, diff, nan, isnan, arange, zeros, ones

# DETRENDER CLASSES
class Detrender(object):
    @abc.abstractmethod
    def detrend(self, y): return
    def polytrend(self, y, order):
        ly = len(y)
        x = arange(1, 1 + ly)
        if ma.isMaskedArray(y):
            mask = y.mask
            line = ma.masked_array(zeros(ly), mask = mask)
            if not mask.all():
                x2 = x[~mask]
                y2 = y[~mask]
                if len(x2) == 1: # single point
                    line[~mask] = y2
                else:
                    coef = polyfit(x2, y2, order)
                    for i in range(order + 1):
                        line[~mask] += coef[i] * x2 ** (order - i)
        else:
            coef = polyfit(x, y, order)
            for i in range(order + 1):
                line += coef[i] * x ** (order - i)
        return line
class NoDetrender(Detrender):
    def detrend(self, y): return y, None # no trend line
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
        dy = ma.masked_where(isnan(dy), dy)
        line = ma.masked_where(isnan(line), line)
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
        dy = ma.masked_where(isnan(dy), dy)
        line = None # no trend line
        return dy, line
class FFDRemoveDetrender(Detrender):
    def detrend(self, y):
        dy = y.copy(); dy[:] = nan
        line = y.copy(); line[:] = nan
        ffd = diff(y) / y[: -1]
        line[1 :] = self.polytrend(ffd, 1)
        dy[1 :] = ffd - line[1 :]
        dy = ma.masked_where(isnan(dy), dy)
        line = ma.masked_where(isnan(line), line)
        return dy, line
class DetrenderWrapper(Detrender):
    def __init__(self, method):
        if method == 'none':
            self.dt = NoDetrender()
        elif method == 'lin':
            self.dt = PolyDetrender(1)
        elif method == 'ffd':
            self.dt = FFDDetrender()
        else:
            raise Exception('Recognized detrend method')
    def detrend(self, y): return self.dt.detrend(y)

# TRANSFORMER CLASSES
class Transformer(object):
    @abc.abstractmethod
    def transform(p1, p2): return
class VSTransformer(Transformer):
    def transform(self, p1, p2):
        p1 -= p1.mean() # decenter
        p1v = p1.var()
        if p1v: p1 *= sqrt(p2.var() / p1v) # scale
        p1 += p2.mean() - p1.mean() # shift
        return p1
class QMTransformer(Transformer):
    def transform(self, p, q):
        p2 = p.copy() # copy
        q2 = q.copy()
        
        q2mean = q2.mean() # target mean
        
        p2 -= ma.median(p2) # center
        q2 -= ma.median(q2)

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
        p2 = p[~p.mask] if ma.isMaskedArray(p) else p.copy()
        return scoreatpercentile(p2, per)

# RETRENDER CLASSES
class Retrender(object):
    @abc.abstractmethod
    def retrend(self, T, R, o): return
class SumRetrender(Retrender):
    def retrend(self, T, R, o): return T + R
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

# WRAPPER CLASS
class BiasCorrecter(object):
    def __init__(self, method):
        detrend_method, transform_method = method.split(' ')
        if transform_method == 'VS':
            self.transformer = VSTransformer()
        elif transform_method == 'QM':
            self.transformer = QMTransformer()
        else:
            raise Exception('Unrecognized transform method')
        if detrend_method in ['ffd-left', 'ffd-right']:
            self.detrenderObs = FFDRemoveDetrender()
            self.detrenderSim = FFDRemoveDetrender()
            if detrend_method == 'ffd-left':
                self.retrender = FFDLeftRetrender() 
            else:
                self.retrender = FFDRightRetrender()
        else:
            self.retrender = SumRetrender()
            self.detrenderSim = PolyDetrender(1)
            if detrend_method == 'lin':
                self.detrenderObs = PolyDetrender(1)
            elif detrend_method == 'quad':
                self.detrenderObs = PolyDetrender(2)
            elif detrend_method == 'moving-ave':
                self.detrenderObs = MovingAveDetrender()
            else:
                raise Exception('Unrecognized detrend method')
    
    def correct(self, sim, obs, tsim, tobs):
        """ Corrects simulated sequence sim using observations obs,
            with tsim and tobs arrays of years """
        # get common times
        time = intersect1d(tsim, tobs)
        if not len(time): raise Exception('No common times')
        
        # detrend observations
        odT, T = self.detrenderObs.detrend(obs)
        odT = toverlap(odT, tobs, time); T = toverlap(T, tobs, time)
        
        # detrend simulations
        sdT = self.detrenderSim.detrend(sim)[0]
        sdT = toverlap(sdT, tsim, time)

        # residuals
        if sdT.mask.all() or odT.mask.all(): # if either is totally masked
            R = ma.masked_array(zeros((ntime,)), mask = ones((ntime,)))
            yhat = R.copy()
        else:
            # transform distribution
            R = self.transformer.transform(sdT, odT)
            
            # combine
            obs_trim = toverlap(obs, tobs, time)
            yhat = self.retrender.retrend(T, R, obs_trim)
        
        return yhat, R

# AGGREGATED BIAS-CORRECTED FILE CLASS
class AggBiasCorrectedFile(object):
    def __init__(self, filename, gadm, time, scen, detrend_methods, corr_methods):
        # create file
        self.filename = filename
        f = nc(filename, 'w', format = 'NETCDF4_CLASSIC')
        
        # create gadm
        f.createDimension('gadm0_index', len(gadm))
        gadmvar = f.createVariable('gadm0_index', 'i4', ('gadm0_index',))
        gadmvar[:] = gadm
        gadmvar.units = 'GADM L0 index'
        gadmvar.long_name = '253 countries'
        
        # create time
        f.createDimension('time', len(time))
        timevar = f.createVariable('time', 'i4', ('time',))
        timevar[:] = time - time[0]
        timevar.units = 'years since {:d}-01-01'.format(int(time[0]))
        timevar.long_name = 'time'
        
        # create scen
        f.createDimension('scen', len(scen))
        scenvar = f.createVariable('scen', 'i4', ('scen',))
        scenvar[:] = range(1, len(scen) + 1)
        scenvar.units = 'mapping'
        scenvar.long_name = ', '.join(scen)
        
        # create detrend
        f.createDimension('detrend_methods', 3)
        detrendvar = f.createVariable('detrend_methods', 'i4', ('detrend_methods',))
        detrendvar[:] = range(1, len(detrend_methods) + 1)
        detrendvar.units = 'mapping'
        detrendvar.long_name = ', '.join(detrend_methods)
             
        # create corr_methods
        f.createDimension('corr_methods', len(corr_methods))
        corr_methodsvar = f.createVariable('corr_methods', 'i4', ('corr_methods',))
        corr_methodsvar[:] = range(1, len(corr_methods) + 1)
        corr_methodsvar.units = 'mapping'
        corr_methodsvar.long_name = ', '.join(corr_methods)

        # close file
        f.close()
    
    def append(self, varname, var, dims, units, longname):
        # append to file
        f = nc(self.filename, 'a')
        yieldvar = f.createVariable(varname, 'f4', dims, zlib = True, shuffle = False, complevel = 9, fill_value = 1e20)
        yieldvar[:] = var
        yieldvar.units = units
        yieldvar.long_name = longname
        f.close()

# HELPER FUNCTIONS
def short2long(cropin):
    # convert GGCMI crop short name to FAOSTAT full crop name, if available
    short_list = ['mai', 'whe', 'soy', 'ric', 'sor', 'mil', \
                  'mgr', 'sug', 'bar', 'oat', 'rap', \
                  'rye', 'sgb']
    long_list = ['maize', 'wheat', 'soybeans', 'rice', 'sorghum', 'millet', \
                 '', 'sugarcane', 'barley', '', 'rapeseed', \
                 '', 'sugar_beet']
    cropout = long_list[short_list.index(cropin)] if cropin in short_list else ''
    return cropout
def toverlap(x, tx, tref):
    ntime = len(tref)
    xref = ma.masked_array(zeros((ntime,)), mask = ones((ntime,)))
    for t in range(ntime):
        tidx = where(tx == tref[t])[0][0]
        xref[t] = x[tidx]
    return xref

# ==================
# SCRIPT STARTS HERE
# ==================

parser = OptionParser()
parser.add_option("-b", "--batch", dest = "batch", default = "1", type = "int",
                  help = "Batch to process")
parser.add_option("-n", "--numbatches", dest = "num_batches", default = "64", type = "int",
                  help = "Total number of batches")
parser.add_option("-d", "--dir", dest = "dir", default = "", type = "string",
                  help = "Directory in which to perform bias correction")
parser.add_option("-r", "--ref", dest = "ref", default = "", type = "string",
                  help = "Reference data netcdf file", metavar = "FILE")
parser.add_option("-o", "--outdir", dest = "outdir", default = "", type = "string",
                  help = "Output directory to save results")
options, args = parser.parse_args()

print 'RUNNING BIASCORRECT'
print '==================='

rootdir = options.dir # root directory
if rootdir[-1] == sep: rootdir = rootdir[: -1] # remove final separator

fileslist = listdir(rootdir)
files = []
for f in fileslist:
    if f.endswith('.nc4') and not f.startswith('faostat') and not 'biascorr' in f:
        files.append(rootdir + sep + f)
nfiles = len(files)

batch = options.batch # find out start and end indices for batch
numbatches = options.num_batches
bz = int(ceil(double(nfiles) / numbatches))
si = bz * (batch - 1)
ei = nfiles if batch == numbatches else min(si + bz, nfiles)

if si >= nfiles: # no work for processor to do
    print 'No jobs for processor to perform. Exiting . . .'
    sys.exit()

# pull reference data
f_ref = nc(options.ref)
gadm_ref = f_ref.variables['gadm0_index'][:]
time_ref = f_ref.variables['time'][:]
time_ref_units = f_ref.variables['time'].units
detrend_methods = f_ref.variables['detrend'].long_name.split(', ')
didx = detrend_methods.index('none')
f_ref.close()

# get reference time
yr0_ref = int(re.findall(r'\d+', time_ref_units)[0])
time_ref = time_ref + yr0_ref

print 'START INDEX = ' + str(si) + ', END IDX = ' + str(ei) + ' . . .'
for f in files[si : ei]:
    print 'PROCESSING FILE: ' + f + ' . . .'
    
    # pull input data
    f_in = nc(f)
    gadm_in = f_in.variables['gadm0_index'][:]
    time_in = f_in.variables['time'][:]
    time_in_units = f_in.variables['time'].units
    scen = f_in.variables['scen'].long_name.split(', ')
    irr_vals = f_in.variables['irr'].long_name.split(', ')
    yield_in = f_in.variables['yield_gadm0'][:]
    sum_idx = irr_vals.index('sum')
    yield_in = yield_in[:, :, :, sum_idx]
    f_in.close()
    
    # pull crop name from file name
    crop_short = split(f)[1].split('_')[3]
    crop_long = short2long(crop_short)
    if crop_long == '':
        print 'Crop not recognized. Exiting . . .'
        sys.exit(1)
    
    # pull reference crop data
    f_ref = nc(options.ref)
    yield_ref = f_ref.variables['yield_' + crop_long][:, :, didx]
    f_ref.close()
    
    # get time
    yr0_in = int(re.findall(r'\d+', time_in_units)[0])
    time_in = time_in + yr0_in - 1
    
    # find common gadm indices
    gadm = intersect1d(gadm_in, gadm_ref)
    ngadm = len(gadm)
    if not ngadm:
        print 'No common gadm indices. Exiting . . .'
        sys.exit(1)
    nscen = len(scen)
    yield_sim_common = ma.ones((ngadm, len(time_in), nscen), fill_value = 1e20)
    yield_ref_common = ma.ones((ngadm, len(time_ref),), fill_value = 1e20)
    for i in range(ngadm):
        gidx_in = where(gadm_in == gadm[i])[0][0]
        yield_sim_common[i] = yield_in[gidx_in]
        gidx_ref = where(gadm_ref == gadm[i])[0][0]
        yield_ref_common[i] = yield_ref[gidx_ref]
    
    detrend_methods = ['none', 'lin', 'ffd']
    corr_methods = ['lin VS', 'quad VS', 'moving-ave VS', 'ffd-left VS', 'ffd-right VS', \
                    'lin QM', 'quad QM', 'moving-ave QM', 'ffd-left QM', 'ffd-right QM']
    
    # get common times
    time = intersect1d(time_in, time_ref)
    ntime = len(time)
    
    # detrend yield
    yield_sim = ma.zeros((ngadm, ntime, nscen, len(detrend_methods)))
    for d in range(len(detrend_methods)):
        dt = DetrenderWrapper(detrend_methods[d])
        for g in range(ngadm):
            for s in range(nscen):
                yield_dt = dt.detrend(yield_sim_common[g, :, s])[0]
                yield_sim[g, :, s, d] = toverlap(yield_dt, time_in, time)
                if yield_sim[g, :, s, d].mask[-1]:
                    yield_sim[g, 1 :, s, d] = yield_sim[g, : -1, s, d] # shift by one
                    yield_sim[g, 0, s, d] = ma.masked
        
    # correct yield
    yield_biascorr = ma.zeros((ngadm, ntime, nscen, len(corr_methods)))    
    yield_retrend = ma.zeros((ngadm, ntime, nscen, len(corr_methods)))
    yield_biascorr.mask = ones(yield_biascorr.shape) # mask all
    yield_retrend.mask = ones(yield_retrend.shape)
    for c in range(len(corr_methods)):
        # iterate over correction methods
        bc = BiasCorrecter(corr_methods[c])
        for g in range(ngadm):
            # iterate over countries
            yref = yield_ref_common[g]
            if yref.mask.all():
                # blank data
                continue
            for s in range(nscen):
                # iterate over scenarios
                ysim = yield_sim_common[g, :, s]
                if not ysim.mask.all(): # if not blank
                    time_sim = time_in.copy()
                    if ysim.mask[-1]: time_sim += 1 # shift by one year
                    # compute correction
                    yhat = ma.masked_array(zeros((ntime,)), mask = ones((ntime,)))
                    R = yhat.copy()
                    retr, var = bc.correct(ysim, yref, time_sim, time_ref)
                    ln = min(ntime, len(retr))
                    yhat[: ln] = retr[: ln]
                    R[: ln] = var[: ln]
                    # save
                    if ysim.mask[-1]:
                        yield_biascorr[g, 1 :, s, c] = R[: -1] # remove last value
                        yield_retrend[g, 1 :, s, c] = yhat[: -1]
                    else:
                        yield_biascorr[g, :, s, c] = R
                        yield_retrend[g, :, s, c] = yhat
    
    # create file
    fn = options.outdir + sep + splitext(split(f)[1])[0] + '.biascorr.nc4'
    fout = AggBiasCorrectedFile(fn, gadm, time, scen, detrend_methods, corr_methods)            
                                        
    # append to file
    fout.append('yield_sim', yield_sim, ('gadm0_index', 'time', 'scen', 'detrend_methods',), 't ha-1 yr-1', 'average detrended gadm yield')
    fout.append('yield_biascorr', yield_biascorr, ('gadm0_index', 'time', 'scen', 'corr_methods',), 't ha-1 yr-1', 'average bias-corrected gadm yield')
    fout.append('yield_retrend', yield_retrend, ('gadm0_index', 'time', 'scen', 'corr_methods',), 't ha-1 yr-1', 'average retrended gadm yield')