#!/usr/bin/env python

# import modules
import warnings
import sys, re, abc
import numpy.ma as ma
from os import listdir
from netCDF4 import Dataset as nc
from optparse import OptionParser
from os.path import split, splitext, sep
from scipy.stats import scoreatpercentile
from numpy import ceil, double, intersect1d, where, sqrt, polyfit, diff, nan, isnan, arange, zeros, ones

# warnings.filterwarnings('error')

# DETRENDER CLASSES
class Detrender(object):
    @abc.abstractmethod
    def detrend(self, y): return
class NoDetrender(Detrender):
    def detrend(self, y): return y, None # no trend line
class LinDetrender(Detrender):
    def detrend(self, y):
        x = arange(1, 1 + len(y))
        coef = polyfit(x, y, 1)
        line = coef[0] * x + coef[1]
        dy = y - line
        return dy, line
class QuadDetrender(Detrender):
    def detrend(self, y):
        x = arange(1, 1 + len(y))
        coef = polyfit(x, y, 2)
        line = coef[0] * x ** 2 + coef[1] * x + coef[2]
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
        dy = y.copy()
        dy[:] = nan
        dy[1 :] = diff(y) / y[: -1]
        dy = ma.masked_where(isnan(dy), dy)
        line = None # no trend line
        return dy, line
class FFDRemoveDetrender(Detrender):
    def detrend(self, y):
        dy = y.copy(); dy[:] = nan
        line = y.copy(); line[:] = nan
        x = arange(2, 1 + len(y))
        ffd = diff(y) / y[: -1]
        coef = polyfit(x, ffd, 1)
        line[1 :] = coef[0] * x + coef[1]
        dy[1 :] = ffd - line[1 :]
        dy = ma.masked_where(isnan(dy), dy)
        line = ma.masked_where(isnan(line), line)
        return dy, line
class DetrenderWrapper(Detrender):
    def __init__(self, method):
        if method == 'none':
            self.dt = NoDetrender()
        elif method == 'lin':
            self.dt = LinDetrender()
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
        p1 *= sqrt(p2.var() / p1.var()) # scale
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
        p2[p2 < 0] *= q25 / p25
        
        p75 = self.__percentile(p2, 75) # 75th
        p2[p2 > 0] *= q75 / p75

        p5 = self.__percentile(p2, 5) # 5th
        p2[p2 < q25] *= q5 / p5
        
        p95 = self.__percentile(p2, 95) # 95th
        p2[p2 > q75] *= q95 / p95

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
        r = zeros((len(T),))
        r[0] = o[0] # initial condition = leftmost point
        # print r[0]
        # print T
        # print R
        for i in range(1, len(r)):
            # try:
                # print i, r[i-1], T[i], R[i]
                r[i] = r[i - 1] * (1. + T[i] + R[i])
            # except:
            #     print T
            #     print R
        return r
class FFDRightRetrender(Retrender):
    def retrend(self, T, R, o):
        r = zeros((len(T),))
        r[-1] = o[-1] # initial condition = rightmost point
        for i in range(-1, -len(r), -1):
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
            self.detrenderSim = LinDetrender()
            if detrend_method == 'lin':
                self.detrenderObs = LinDetrender()
            elif detrend_method == 'quad':
                self.detrenderObs = QuadDetrender()
            elif detrend_method == 'moving-ave':
                self.detrenderObs = MovingAveDetrender()
            else:
                raise Exception('Unrecognized detrend method')
    
    def correct(self, sim, obs):
        """ Corrects simulated sequence sim using observations obs """
        # detrend observations
        odT, T = self.detrenderObs.detrend(obs)
        
        # detrend simualtions
        sdT = self.detrenderSim.detrend(sim)[0]

        if isnan(odT).sum() == len(odT):
            print odT
            print obs
            raise Exception('blah')
        if isnan(sdT).sum() == len(sdT):
            print sdT
            print sim
            raise Exception('blah2')

        # residuals
        # print 'sdT =', sdT
        # print 'odT =', odT
        if sdT.mask.all() or odT.mask.all(): # if either is totally masked
            R = ma.masked_array(zeros((len(T),)), mask = ones((len(T),)))
            yhat = R.copy()
        else:
            # transform distribution
            R = self.transformer.transform(sdT, odT)
            
            # combine
            yhat = self.retrender.retrend(T, R, obs)
        
        return yhat, T, R

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
        timevar[:] = time
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

# HELPER FUNCTION
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
    
    # find common times and gadm indices
    gadm = intersect1d(gadm_in, gadm_ref)
    ngadm = len(gadm)
    if not ngadm:
        print 'No common gadm indices. Exiting . . .'
        sys.exit(1)
    time = intersect1d(time_in, time_ref)
    ntime = len(time)
    if not ntime:
        print 'No common times. Exiting . . .'
        sys.exit(1)
    nscen = len(scen)
    yield_sim_common = ma.ones((ngadm, ntime, nscen), fill_value = 1e20)
    yield_ref_common = ma.ones((ngadm, ntime,), fill_value = 1e20)
    for i in range(ngadm):
        for j in range(ntime):
            gidx_in = where(gadm_in == gadm[i])[0][0]
            tidx_in = where(time_in == time[j])[0][0]
            yield_sim_common[i, j] = yield_in[gidx_in, tidx_in]
            gidx_ref = where(gadm_ref == gadm[i])[0][0]
            tidx_ref = where(time_ref == time[j])[0][0]
            yield_ref_common[i, j] = yield_ref[gidx_ref, tidx_ref]
    
    detrend_methods = ['none', 'lin', 'ffd']
    corr_methods = ['lin VS', 'quad VS', 'moving-ave VS', 'ffd-left VS', 'ffd-right VS', \
                    'lin QM', 'quad QM', 'moving-ave QM', 'ffd-left QM', 'ffd-right QM']
    
    # detrend yield
    yield_sim = ma.zeros((ngadm, ntime, nscen, len(detrend_methods))) 
    for d in range(len(detrend_methods)):
        dt = DetrenderWrapper(detrend_methods[d])
        for g in range(ngadm):
            for s in range(nscen):
                yield_sim[g, :, s, d] = dt.detrend(yield_sim_common[g, :, s])[0]
        
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
                    # compute correction
                    # print corr_methods[c]
                    # print 'ysim =', ysim
                    # print 'yref =', yref
                    yhat, T, R = bc.correct(ysim, yref)
                    # save
                    yield_biascorr[g, :, s, c] = R
                    yield_retrend[g, :, s, c] = yhat
    
    # create file
    fn = options.outdir + sep + splitext(split(f)[1])[0] + '.biascorr.nc4'
    fout = AggBiasCorrectedFile(fn, gadm, time, scen, detrend_methods, corr_methods)            
                                        
    # append to file
    fout.append('yield_sim', yield_sim, ('gadm0_index', 'time', 'scen', 'detrend_methods',), 't ha-1 yr-1', 'average detrended gadm yield')
    fout.append('yield_biascorr', yield_biascorr, ('gadm0_index', 'time', 'scen', 'corr_methods',), 't ha-1 yr-1', 'average bias-corrected gadm yield')
    fout.append('yield_retrend', yield_retrend, ('gadm0_index', 'time', 'scen', 'corr_methods',), 't ha-1 yr-1', 'average retrended gadm yield')