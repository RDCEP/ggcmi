#!/usr/bin/env python

# add paths
import os, sys
for p in os.environ['PATH'].split(':'): sys.path.append(p)

# import modules
from re import findall
from os.path import split
from itertools import product
from numpy.ma import masked_array
from optparse import OptionParser
from netCDF4 import Dataset as nc
from metrics import MetricsWrapper
from filespecs import MultimetricsEnsembleFile
from numpy import where, ones, zeros, arange, logical_and

parser = OptionParser()
parser.add_option("-i", "--infile", dest = "infile", default = "", type = "string",
                  help = "Input model ensemble file", metavar = "FILE")
parser.add_option("-r", "--reffile", dest = "reffile", default = "", type = "string",
                  help = "Reference data netcdf file", metavar = "FILE")
parser.add_option("-a", "--agglvl", dest = "agglvl", default = "gadm0", type = "string",
                  help = "Aggregation level (e.g., gadm0, fpu, kg)")
parser.add_option("-m", "--metric", dest = "metric", default = "rmse", type = "string",
                  help = "Metric name")
parser.add_option("-u", "--munits", dest = "munits", default = "t ha-1 yr-1", type = "string",
                  help = "Metric units")
parser.add_option("-l", "--mlongname", dest = "mlongname", default = "root mean squared error", type = "string",
                  help = "Metric long name")
parser.add_option("-o", "--outfile", dest = "outfile", default = "", type = "string",
                  help = "Output file")
options, args = parser.parse_args()

infile    = options.infile
reffile   = options.reffile
agglvl    = options.agglvl
metric    = options.metric
munits    = options.munits
mlongname = options.mlongname
outfile   = options.outfile

tranges = ['full']
ntimes  = len(tranges)

crop = split(infile)[1].split('_')[3]

with nc(reffile) as fref:
    aref        = fref.variables[agglvl][:]
    aggunits    = fref.variables[agglvl].units
    agglongname = fref.variables[agglvl].long_name
    dtref       = fref.variables['dt'].long_name.split(', ')
    mpref       = fref.variables['mp'].long_name.split(', ')
    tref        = fref.variables['time'][:]
    tref_units  = fref.variables['time'].units

    var = 'yield_' + crop
    if var in fref.variables:
        yield_ref = fref.variables[var][:]
    else:
        print 'Crop %s unavailable in reference file %s. Exiting . . .' % (crop, reffile)
        sys.exit()

with nc(infile) as fin:
    ain       = fin.variables[agglvl][:]
    dt        = fin.variables['dt'].long_name.split(', ')
    mp        = fin.variables['mp'].long_name.split(', ')
    cr        = fin.variables['cr'].long_name.split(', ')
    nm        = fin.variables['nm'].size
    wt        = fin.variables['wt'].long_name.split(', ')
    yield_in  = fin.variables['yield_detrend'][:]
    tin       = fin.variables['time'][:]
    tin_units = fin.variables['time'].units

tref += int(findall(r'\d+', tref_units)[0]) # get reference time
tin  += int(findall(r'\d+', tin_units)[0])  # get simulation time

naggs, ndt, nmp, ncr, nwt = len(ain), len(dt), len(mp), len(cr), len(wt)

sh = (naggs, ndt, nmp, ncr, ntimes, nm, nwt)

times = [tin]
dtidx, mpidx = dtref.index('none'), mpref.index('true')

mobj = MetricsWrapper(metric)
mmat = masked_array(zeros(sh), mask = ones(sh))
for t in range(ntimes):
    tmin = max([tin[0],  tref[0],  times[t][0]])
    tmax = min([tin[-1], tref[-1], times[t][-1]])

    yield_refc = yield_ref[:, logical_and(tref >= tmin, tref <= tmax)]
    yield_inc  = yield_in[:,  logical_and(tin  >= tmin, tin  <= tmax)]

    for d, m, c in product(range(ndt), range(nmp), range(ncr)):
        for a, n, w in product(range(naggs), range(nm), range(nwt)):
            refidx = dtref.index(dt[d])
            aidx = where(aref == ain[a])[0][0]

            dref     = yield_refc[aidx, :, refidx, m]
            drefnone = yield_refc[aidx, :, dtidx, mpidx]
            dsim     = yield_inc[a, :, d, m, c, n, w]

            mmat[a, d, m, c, t, n, w] = mobj.eval(dsim, dref, drefnone, arange(tmin, tmax + 1))

fout = MultimetricsEnsembleFile(outfile, ain, agglvl, aggunits, agglongname, tranges, dt, mp, cr, nm, wt)
fout.append(metric, mmat, (agglvl, 'dt', 'mp', 'cr', 'time_range', 'nm', 'wt'), munits, mlongname) # append to file
