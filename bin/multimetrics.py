#!/usr/bin/env python

import re, sys
from os import sep, listdir
from itertools import product
from optparse import OptionParser
from netCDF4 import Dataset as nc
from numpy.ma import masked_array, corrcoef, masked
from numpy import array, unique, where, ones, zeros, sqrt, intersect1d, ceil, double

def createnc(filename, gadmindices, scen, detrmethods, corrmethods, timenames, climate, crop, gadmunits, gadmlong):
    with nc(filename, 'w', format = 'NETCDF4_CLASSIC') as f:
        f.createDimension('gadm0_index', len(gadmindices))
        gadmvar = f.createVariable('gadm0_index', 'i4', 'gadm0_index')
        gadmvar[:] = gadmindices
        gadmvar.units = gadmunits
        gadmvar.long_name = gadmlong
        f.createDimension('scen', len(scen))
        scenvar = f.createVariable('scen', 'i4', 'scen')
        scenvar[:] = range(1, len(scen) + 1)
        scenvar.units = 'mapping'
        scenvar.long_name = ', '.join(scen)
        f.createDimension('detr_methods', len(detrmethods))
        detrvar = f.createVariable('detr_methods', 'i4', 'detr_methods')
        detrvar[:] = range(1, len(detrmethods) + 1)
        detrvar.units = 'mapping'
        detrvar.long_name = ', '.join(detrmethods)
        f.createDimension('corr_methods', len(corrmethods))
        corrvar = f.createVariable('corr_methods', 'i4', 'corr_methods')
        corrvar[:] = range(1, len(corrmethods) + 1)
        corrvar.units = 'mapping'
        corrvar.long_name = ', '.join(corrmethods)
        f.createDimension('time_range', len(timenames))
        timevar = f.createVariable('time_range', 'i4', 'time_range')
        timevar[:] = range(1, len(timenames) + 1)
        timevar.units = 'mapping'
        timevar.long_name = ', '.join(timenames)
        f.createDimension('climate', len(climate))
        climatevar = f.createVariable('climate', 'i4', 'climate')
        climatevar[:] = range(1, len(climate) + 1)
        climatevar.units = 'mapping'
        climatevar.long_name = ', '.join(climate)
        f.createDimension('crop', len(crop))
        cropvar = f.createVariable('crop', 'i4', 'crop')
        cropvar[:] = range(1, len(crop) + 1)
        cropvar.units = 'mapping'
        cropvar.long_name = ', '.join(crop)

def computeTSCorr(a, b):
    return corrcoef(a, b)[0, 1]
def computeVarRatio(a, b):
    if not b.mask.all():
        a_mean = a.mean()
        b_std  = b.std()
        if a_mean * b_std: # no zeros
            return (a.std() * b.mean() / a_mean * b_std)
        else:
            return masked
    else:
        return masked
def computeRMSE(a, b, c = None):
    rmse = sqrt(((a - b) ** 2).mean())
    if c is not None:
        if not c.mask.all():
            norm = c.mean()
            if norm:
                rmse /= norm
    return rmse

def restricttime(d, t, tref):
    tcommon = intersect1d(t, tref)
    ntime = len(tcommon)
    sh = array(d.shape)
    sh[0] = ntime
    dout = masked_array(zeros(sh), mask = ones(sh))
    for i in range(ntime):
        tidx = where(t == tcommon[i])[0][0]
        dout[i] = d[tidx]
    return dout

def short2long(cropin):
    long_map = {'mai': 'maize', 'whe': 'wheat', 'soy': 'soybeans', 'ric': 'rice', \
                'sor': 'sorghum', 'mil': 'millet', 'mgr': '', 'sug': 'sugarcane', \
                'bar': 'barley', 'oat': '', 'rap': 'rapeseed', 'rye': '', \
                'sgb': 'sugar_beet'}
    if cropin in long_map:
        return long_map[cropin]
    else:
        raise Exception('Crop not recognized')

parser = OptionParser()
parser.add_option("-b", "--batch", dest = "batch", default = "1", type = "int",
                  help = "Batch to process")
parser.add_option("-n", "--numbatches", dest = "num_batches", default = "64", type = "int",
                  help = "Total number of batches")
parser.add_option("-d", "--dir", dest = "dir", default = "", type = "string",
                  help = "Directory in which to perform multimetrics evaluation")
parser.add_option("-r", "--ref", dest = "ref", default = "", type = "string",
                  help = "Reference data netcdf file", metavar = "FILE")
parser.add_option("-o", "--outdir", dest = "outdir", default = "", type = "string",
                  help = "Output directory to save results")
options, args = parser.parse_args()

batch, numbatches = options.batch, options.num_batches
indir, outdir = options.dir, options.outdir
reffile = options.ref

files = listdir(indir)

models  = unique(array([f.split('_')[0] for f in files]))
nmodels = len(models)

bz = int(ceil(double(nmodels) / numbatches))
si = bz * (batch - 1)
ei = nmodels if batch == numbatches else min(si + bz, nmodels)
if si >= nmodels:
    print 'No jobs for processor to perform. Exiting . . .'
    sys.exit()

timenames = ['full', '1980-2001', '1980-2009']

for m in models[si : ei]:
    modelfiles = [f for f in files if f.startswith(m)]

    climate = unique(array([f.split('_')[1] for f in modelfiles]))
    crops   = unique(array([f.split('_')[3] for f in modelfiles]))

    with nc(indir + sep + modelfiles[0]) as f:
        gadmindices = f.variables['gadm0_index'][:]
        gadmunits   = f.variables['gadm0_index'].units
        gadmlong    = f.variables['gadm0_index'].long_name
        detrmethodssim = f.variables['detr_methods'].long_name.split(', ')
        corrmethodssim = f.variables['corr_methods'].long_name.split(', ')
    ndetr, ncorr = len(detrmethodssim), len(corrmethodssim)

    scen = []
    for f in modelfiles:
        with nc(indir + sep + f) as fm:
            for s in fm.variables['scen'].long_name.split(', '): scen.append(s)
    scen = unique(array(scen))

    outfile = outdir + sep + m + '_metrics.nc4'
    createnc(outfile, gadmindices, scen, detrmethodssim, corrmethodssim, timenames, climate, crops, gadmunits, gadmlong)

    sh = (len(gadmindices), len(scen), ndetr, ncorr, len(timenames), len(climate), len(crops))
    tscorr   = masked_array(zeros(sh), mask = ones(sh))
    varratio = masked_array(zeros(sh), mask = ones(sh))
    rmse     = masked_array(zeros(sh), mask = ones(sh))

    for f in modelfiles:
        climatename  = f.split('_')[1]
        cropname     = f.split('_')[3]
        cropnamelong = short2long(cropname)

        climateidx = where(climate == climatename)[0][0]
        cropidx    = where(crops == cropname)[0][0]

        with nc(indir + sep + f) as fm:
            yieldsim = fm.variables['yield_sim'][:]
            sc = fm.variables['scen'].long_name.split(', ')
            yr0 = int(re.findall(r'\d+', fm.variables['time'].units)[0])
            timesim = fm.variables['time'][:] + yr0
        with nc(reffile) as ref:
            gadmindicesref = ref.variables['gadm0_index'][:]
            yieldref = ref.variables['yield_' + cropnamelong][:]
            corrmethodsref = ref.variables['detrend'].long_name.split(', ')
            yr0 = int(re.findall(r'\d+', ref.variables['time'].units)[0])
            timeref = ref.variables['time'][:] + yr0

        yieldref = yieldref.transpose((1, 0, 2)) # make time first dimension
        yieldsim = yieldsim.transpose((1, 0, 2, 3, 4))

        times = [timesim, range(1980, 2002), range(1980, 2010)]
        refnoneidx = corrmethodsref.index('none')

        for t in range(len(times)):
            tcommon = intersect1d(timesim, intersect1d(timeref, times[t]))
            yieldrefcommon = restricttime(yieldref, timeref, tcommon)
            yieldsimcommon = restricttime(yieldsim, timesim, tcommon)
            for d, c in product(range(ndetr), range(ncorr)):
                refidx = corrmethodsref.index(detrmethodssim[d])
                yieldrefcommondt = yieldrefcommon[:, :, refidx]
                yieldrefcommonnone = yieldrefcommon[:, :, refnoneidx]
                yieldsimcommondt = yieldsimcommon[:, :, :, d, c]
                for g, s in product(range(len(gadmindices)), range(len(sc))):
                    gidx = where(gadmindicesref == gadmindices[g])[0][0]
                    scenidx = where(scen == sc[s])[0][0]
                    dref, drefnone = yieldrefcommondt[:, gidx], yieldrefcommonnone[:, gidx]
                    dsim = yieldsimcommondt[:, g, s]
                    tscorr[g, scenidx, d, c, t, climateidx, cropidx]   = computeTSCorr(dsim, dref)
                    varratio[g, scenidx, d, c, t, climateidx, cropidx] = computeVarRatio(dsim, dref)
                    rmse[g, scenidx, d, c, t, climateidx, cropidx]     = computeRMSE(dsim, dref, drefnone)

    with nc(outfile, 'a') as f:
        tscorrvar = f.createVariable('tscorr', 'f4', ('gadm0_index', 'scen', 'detr_methods', 'corr_methods', 'time_range', 'climate', 'crop'), fill_value = 1e20, zlib = True, complevel = 9)
        tscorrvar[:] = tscorr
        varratiovar = f.createVariable('varratio', 'f4', ('gadm0_index', 'scen', 'detr_methods', 'corr_methods', 'time_range', 'climate', 'crop'), fill_value = 1e20, zlib = True, complevel = 9)
        varratiovar[:] = varratio
        rmsevar = f.createVariable('rmse', 'f4', ('gadm0_index', 'scen', 'detr_methods', 'corr_methods', 'time_range', 'climate', 'crop'), fill_value = 1e20, zlib = True, complevel = 9)
        rmsevar[:] = rmse