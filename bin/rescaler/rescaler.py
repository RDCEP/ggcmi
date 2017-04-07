#!/usr/bin/env python

import os, sys
for p in os.environ['PATH'].split(':'): sys.path.append(p)

from re import findall
from netCDF4 import Dataset
from argparse import ArgumentParser
from os.path import basename, splitext, isfile
from numpy.ma import masked_array, masked_where
from numpy import ones, zeros, where, isnan, arange, nanmean, isfinite
from filespecs import RescaledFile
import numpy
import yaml
import warnings

def is_float(fl):
    if not isinstance(fl, numpy.floating):
        return False
    if not isfinite(fl):
        return False
    return True

parser = ArgumentParser()
parser.add_argument("-a", "--agglvl", required=True, help="Aggregation level (e.g., gadm0, fpu, kg)")
parser.add_argument("-b", "--bcfile", required=True, help="Bias corrected file")
parser.add_argument("-i", "--irfile", required=True, help="Input irrigated spatial file")
parser.add_argument("-m", "--mkfile", required=True, help="Mask file")
parser.add_argument("-o", "--outfile", required=True, help="Output rescaled spatial file")
parser.add_argument("-p", "--params", required=True, help="YAML Parameter File")
parser.add_argument("-r", "--rffile", required=True, help="Input rainfed spatial file")
parser.add_argument("-v", "--vname", required=True, help="Variable name")
parser.add_argument("-w", "--wtfile", required=True, help="Weights file")
args = parser.parse_args()
print args

agglvl = args.agglvl
bcfile = args.bcfile
irfile = args.irfile
mkfile = args.mkfile
outfile = args.outfile
rffile = args.rffile
vname = args.vname
wtfile = args.wtfile

params = yaml.load(open(args.params, 'r'))
dt = params['rescaler']['dt']
mp = params['rescaler']['mp']
cr = params['rescaler']['cr']
scen = params['rescaler']['scen']

fsplit = basename(irfile).split('_')
if len(fsplit) == 10:
    scen = fsplit[3]
else:
    scen = '_'.join([fsplit[i] for i in [3, 5]])

# load weights
with Dataset(wtfile) as f:
    lats = f.variables['lat'][:]
    lons = f.variables['lon'][:]
    wir = f.variables['irrigated'][:]
    wrf = f.variables['rainfed'][:]

# load mask
with Dataset(mkfile) as f:
    aggmap = f.variables[agglvl][:]

# load bias corrected file
with Dataset(bcfile) as f:
    agglevels = f.variables[agglvl][:]
    time = f.variables['time'][:]
    tunits = f.variables['time'].units
    snames = f.variables['scen'].long_name.split(', ')
    if scen not in snames:
        print "Cannot find variable %s in %s. Exiting" % (scen, bcfile)
        sys.exit(0)
    else:
        sidx = snames.index(scen)
    detrend = f.variables['yield_detrend'][:]
    try:
        dtidx = f.variables['dt'].long_name.split(', ').index(dt)
        dtnone = f.variables['dt'].long_name.split(', ').index('none')
        mpidx = f.variables['mp'].long_name.split(', ').index(mp)
        cridx = f.variables['cr'].long_name.split(', ').index(cr)
        scenidx = f.variables['scen'].long_name.split(', ').index(scen)
    except ValueError as e:
        print "%s is missing a required dimension (%s, %s, %s, or %s)" % (bcfile, dtidx, 'none', mpidx, cridx)
        print e
        sys.exit(0)

# Get time in years
time += int(findall(r'\d+', tunits)[0])

y1, y2 = [int(y) for y in findall(r'\d+', splitext(basename(irfile))[0])[-2:]]
ftime = arange(y1, y2 + 1)

with Dataset(irfile) as f:
    var = f.variables[vname]
    # Happens with epic-test
    if len(var) != len(ftime):
        ftime = f.variables['time'][:] + \
                int(findall(r'\d+', f.variables['time'].units)[0]) - 1

# Time bounds
tmin = max(time[0], ftime[0])
tmax = min(time[-1], ftime[-1])

# Time indexes
btidx0 = where(time == tmin)[0][0]
btidx1 = where(time == tmax)[0][0] + 1
ftidx0 = where(ftime == tmin)[0][0]
ftidx1 = where(ftime == tmax)[0][0] + 1

time = time[btidx0:btidx1]
detrend = masked_where(isnan(detrend), detrend)

# Get dimensions
(nlats, nlons) = aggmap.shape
nt = len(time)
nirr = 3

varr = masked_array(zeros((nt, nlats, nlons, nirr)), mask=ones((nt, nlats, nlons, nirr)))

with Dataset(irfile) as f:
    var = f.variables[vname]
    if 'units' in var.ncattrs():
        units = var.units
    else:
        units = ''
    if 'long_name' in var.ncattrs():
        lname = var.long_name
    else:
        lname = ''
    varr[:, :, :, 0] = var[ftidx0:ftidx1]

if isfile(rffile):
    with Dataset(rffile) as f:
        varr[:, :, :, 1] = f.variables[vname][ftidx0:ftidx1]

varr[isnan(varr)] = 0.
correction_data = numpy.full((len(lats), len(lons)), 1.e+20)
warnings.filterwarnings('ignore')

for i in range(len(agglevels)):
    latidxs, lonidxs = where(aggmap == agglevels[i])
    detrend_ma = nanmean(detrend[i, :, scenidx, dtidx, mpidx, cridx])
    detrend_none = nanmean(detrend[i, :, scenidx, dtnone, mpidx, cridx])
    correction = detrend_ma / detrend_none

    if is_float(correction):
        print "level=%s, detrend_ma=%s, detrend_none=%s, correction=%s" % (agglevels[i], detrend_ma, detrend_none,
                                                                           correction)
        correction_data[latidxs, lonidxs] = correction
        varr[:, latidxs, lonidxs, 0:2] *= correction
    varr[:, latidxs, lonidxs, 2] = varr[:, latidxs, lonidxs, 0] + varr[:, latidxs, lonidxs, 1]

fout = RescaledFile(outfile, time, lats, lons, ['ir', 'rf', 'sum'])
fout.append(vname, varr, ('time', 'lat', 'lon', 'irr'), units, lname)
