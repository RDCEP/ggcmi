#!/usr/bin/env python

# add paths
import os, sys
for p in os.environ['PATH'].split(':'): sys.path.append(p)

# import modules
from os import sep, listdir
from fnmatch import fnmatch
from os.path import basename
from ensembler import Ensembler
from optparse import OptionParser
from filespecs import ModelEnsembleFile

parser = OptionParser()
parser.add_option("-d", "--indir", dest = "indir", default = "", type = "string",
                  help = "Directory of files to ensemble")
parser.add_option("-m", "--metricsdir", dest = "metricsdir", default = "", type = "string",
                  help = "Directory of metric files")
parser.add_option("-a", "--agglvl", dest = "agglvl", default = "gadm0", type = "string",
                  help = "Aggregation level (e.g., gadm0, fpu, kg)")
parser.add_option("-w", "--weather", dest = "weather", default = "agmerra", type = "string",
                  help = "Weather dataset")
parser.add_option("-c", "--crop", dest = "crop", default = "mai", type = "string",
                  help = "Crop name")
parser.add_option("--metric", dest = "metric", default = "rmse", type = "string",
                  help = "Metric to weight data with (default = rmse)")
parser.add_option("-o", "--outdir", dest = "outdir", default = "", type = "string",
                  help = "Output directory to save results")
options, args = parser.parse_args()

indir      = options.indir
metricsdir = options.metricsdir
agglvl     = options.agglvl
weather    = options.weather
crop       = options.crop
metric     = options.metric
outdir     = options.outdir

files = listdir(indir)
bcfiles = [f for f in files if fnmatch(f, '*_%s_*_%s_*' % (weather, crop))]
bcfiles = [indir + sep + f for f in bcfiles]

files = listdir(metricsdir)
mmfiles = [basename(f).replace('biascorr', 'multimetrics') for f in bcfiles]
mmfiles = [metricsdir + sep + f for f in mmfiles]

nmodels = len(bcfiles)

if not nmodels:
    print 'No files for climate %s, crop %s. Exiting . . .' % (weather, crop)
    sys.exit()

ensembler = Ensembler(bcfiles, mmfiles, agglvl, metric)

yield_detr, yield_retr, model_order, model_weights = ensembler.average()

tmin        = ensembler.tmin
tmax        = ensembler.tmax
time        = ensembler.time
aggs        = ensembler.aggs
aggunits    = ensembler.aggunits
agglongname = ensembler.agglongname
dt          = ensembler.dt
mp          = ensembler.mp
cr          = ensembler.cr
nm          = ensembler.nm
models      = ensembler.models
metricunits = ensembler.metricunits

outfile = outdir + sep + '%s_%s_hist_%s_annual_%d_%d.ensemble.nc4' % (metric, weather, crop, tmin, tmax)
fout = ModelEnsembleFile(outfile, metric, aggs, agglvl, aggunits, agglongname, time, dt, mp, cr, nm)

fout.append('yield_detrend', yield_detr,    (agglvl, 'time', 'dt', 'mp', 'cr', 'nm', 'wt'), 't ha-1 yr-1', 'average ensemble detrended yield')
fout.append('yield_retrend', yield_retr,    (agglvl, 'time', 'dt', 'mp', 'cr', 'nm', 'wt'), 't ha-1 yr-1', 'average ensemble retrended yield')
fout.append('model_order',   model_order,   (agglvl, 'dt', 'mp', 'cr', 'nm'),                'mapping',    ', '.join(models))
fout.append('model_weights', model_weights, (agglvl, 'dt', 'mp', 'cr', 'nm'),                metricunits,  metric)