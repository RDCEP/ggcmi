#!/usr/bin/env python

# add paths
import os, sys
for p in os.environ['PATH'].split(':'): sys.path.append(p)

# import modules
from re import findall
from itertools import product
from netCDF4 import Dataset as nc
from optparse import OptionParser
from numpy.ma import masked_array
from filespecs import BiasCorrectFile
from biascorrecter import BiasCorrecter
from os.path import split, splitext, sep
from numpy import intersect1d, zeros, ones

parser = OptionParser()
parser.add_option("-i", "--infile", dest = "infile", default = "", type = "string",
                  help = "Input aggregated file", metavar = "FILE")
parser.add_option("-r", "--reffile", dest = "reffile", default = "", type = "string",
                  help = "Reference data netcdf file", metavar = "FILE")
parser.add_option("-o", "--outdir", dest = "outdir", default = "", type = "string",
                  help = "Output directory to save results")
options, args = parser.parse_args()

infile  = options.infile
reffile = options.reffile
outdir  = options.outdir

dt = ['none', 'lin', 'quad', 'ma', 'ffdtr'] # methods
mp = ['true', 'false']
cr = ['none', 'VS', 'MS']
ndt, nmp, ncr = len(dt), len(mp), len(cr)

crop = split(infile)[1].split('_')[3] # pull crop name from file name

with nc(reffile) as fref: # pull reference data
    gref       = fref.variables['gadm0'][:]
    tref       = fref.variables['time'][:]
    tref_units = fref.variables['time'].units
    dtidx      = fref.variables['dt'].long_name.split(', ').index('none')
    mpidx      = fref.variables['mp'].long_name.split(', ').index('true')
    yield_ref  = fref.variables['yield_' + crop][:, :, dtidx, mpidx]
tref += int(findall(r'\d+', tref_units)[0]) # get reference time

with nc(infile) as fin: # pull input data
    gin       = fin.variables['gadm0_index'][:]
    tin       = fin.variables['time'][:]
    tin_units = fin.variables['time'].units
    scen      = fin.variables['scen'].long_name.split(', ')
    sum_idx   = fin.variables['irr'].long_name.split(', ').index('sum')
    yield_in  = fin.variables['yield_gadm0'][:, :, :, sum_idx]
tin += int(findall(r'\d+', tin_units)[0]) - 1 # get time

gadm = intersect1d(gin, gref) # find common gadm indices
time = intersect1d(tin, tref) # find common times
ngadm, ntime, nscen = len(gadm), len(time), len(scen)
if not ngadm: raise Exception('No common gadm indices')

yield_sim_common = masked_array(zeros((ngadm, len(tin), nscen)), mask = ones((ngadm, len(tin), nscen)))
yield_ref_common = masked_array(zeros((ngadm, len(tref))), mask = ones((ngadm, len(tref))))
for i in range(ngadm):
    yield_sim_common[i] = yield_in[list(gin).index(gadm[i])]
    yield_ref_common[i] = yield_ref[list(gref).index(gadm[i])]

sh = (ngadm, ntime, nscen, ndt, nmp, ncr)
yield_sim  = masked_array(zeros(sh), mask = ones(sh))
yield_retr = masked_array(zeros(sh), mask = ones(sh))
for g, s in product(range(ngadm), range(nscen)):
    yref, ysim = yield_ref_common[g], yield_sim_common[g, :, s]
    if not yref.mask.all() and not ysim.mask.all():
        for d, m, c in product(range(ndt), range(nmp), range(ncr)):
            bc = BiasCorrecter(dt[d], mp[m], cr[c])
            yhat, R = bc.correct(ysim, yref, tin, tref)
            yield_sim[g, :, s, d, m, c]  = R
            yield_retr[g, :, s, d, m, c] = yhat

fn = outdir + sep + splitext(split(infile)[1])[0] + '.biascorr.nc4' # create file
fout = BiasCorrectFile(fn, gadm, time, scen, dt, mp, cr)            
fout.append('yield_detrend', yield_sim,  ('gadm0', 'time', 'scen', 'dt', 'mp', 'cr'), 't ha-1 yr-1', 'average detrended gadm yield') # append to file
fout.append('yield_retrend', yield_retr, ('gadm0', 'time', 'scen', 'dt', 'mp', 'cr'), 't ha-1 yr-1', 'average bias-corrected gadm yield')