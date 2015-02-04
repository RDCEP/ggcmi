#!/usr/bin/env python

# add paths
import os, sys
for p in os.environ['PATH'].split(':'): sys.path.append(p)

# import modules
from re import findall
from itertools import product
from netCDF4 import Dataset as nc
from optparse import OptionParser
from filespecs import BiasCorrectFile
from biascorrecter import BiasCorrecter
from os.path import split, splitext, sep
from numpy.ma import masked_array, reshape
from numpy import intersect1d, zeros, ones

parser = OptionParser()
parser.add_option("-i", "--infile", dest = "infile", default = "", type = "string",
                  help = "Input aggregated file", metavar = "FILE")
parser.add_option("-r", "--reffile", dest = "reffile", default = "", type = "string",
                  help = "Reference data netcdf file", metavar = "FILE")
parser.add_option("-a", "--agglvl", dest = "agglvl", default = "gadm0", type = "string",
                  help = "Aggregation level (e.g., gadm0, fpu, kg)")
parser.add_option("-o", "--outdir", dest = "outdir", default = "", type = "string",
                  help = "Output directory to save results")
options, args = parser.parse_args()

infile  = options.infile
reffile = options.reffile
agglvl  = options.agglvl
outdir  = options.outdir

dt = ['none'] # methods
mp = ['true']
cr = ['mean-scale']
ndt, nmp, ncr = len(dt), len(mp), len(cr)

scen = ['default'] # one scenario

crop = split(infile)[1].split('_')[5] # pull crop name from file name

with nc(reffile) as fref: # pull reference data
    aref        = fref.variables[agglvl][:]
    aggunits    = fref.variables[agglvl].units
    agglongname = fref.variables[agglvl].long_name 
    tref        = fref.variables['time'][:]
    tref_units  = fref.variables['time'].units
    dtidx       = fref.variables['dt'].long_name.split(', ').index('none')
    mpidx       = fref.variables['mp'].long_name.split(', ').index('true')

    var = 'yield_' + crop
    if var in fref.variables:
        yield_ref = fref.variables[var][:, :, dtidx, mpidx]
    else:
        print 'Crop %s unavailable in reference file %s. Exiting . . .' % (crop, reffile)
        sys.exit()

with nc(infile) as fin: # pull input data
    ain       = fin.variables[agglvl][:]
    tin       = fin.variables['time'][:]
    tin_units = fin.variables['time'].units
    sum_idx   = fin.variables['irr'].long_name.split(', ').index('sum')
    yield_in = fin.variables['yield_' + crop + '_' + agglvl][:, :, sum_idx]
    yield_in = reshape(yield_in, (len(tin), len(ain), 1))

tref += int(findall(r'\d+', tref_units)[0]) # get reference time
tin  += int(findall(r'\d+', tin_units)[0])  # get simulation time

aggs = intersect1d(ain, aref) # find common gadm indices
naggs, ntime, nscen = len(aggs), len(tin), len(scen)
if not naggs: raise Exception('No common aggregates')

yield_sim_common = masked_array(zeros((naggs, len(tin), nscen)), mask = ones((naggs, len(tin), nscen)))
yield_ref_common = masked_array(zeros((naggs, len(tref))), mask = ones((naggs, len(tref))))
for i in range(naggs):
    yield_sim_common[i] = yield_in[:, list(ain).index(aggs[i])]
    yield_ref_common[i] = yield_ref[list(aref).index(aggs[i])]

sh = (naggs, ntime, nscen, ndt, nmp, ncr)
yield_detr = masked_array(zeros(sh), mask = ones(sh))
yield_retr = masked_array(zeros(sh), mask = ones(sh))
for g, s in product(range(naggs), range(nscen)):
    yref, ysim = yield_ref_common[g], yield_sim_common[g, :, s]
    if not yref.mask.all() and not ysim.mask.all():
        for d, m, c in product(range(ndt), range(nmp), range(ncr)):
            bc = BiasCorrecter(dt[d], mp[m], cr[c])
            detr, retr = bc.correct(ysim, yref, tin, tref)
            yield_detr[g, :, s, d, m, c] = detr
            yield_retr[g, :, s, d, m, c] = retr

fn = outdir + sep + splitext(split(infile)[1])[0] + '.biascorr.nc4' # create file
fout = BiasCorrectFile(fn, aggs, agglvl, aggunits, agglongname, tin, scen, dt, mp, cr)
fout.append('yield_detrend', yield_detr, (agglvl, 'time', 'scen', 'dt', 'mp', 'cr'), 't ha-1 yr-1', 'average detrended yield') # append to file
fout.append('yield_retrend', yield_retr, (agglvl, 'time', 'scen', 'dt', 'mp', 'cr'), 't ha-1 yr-1', 'average retrended yield')