#!/usr/bin/env python

from argparse import ArgumentParser
from glob import glob
from netCDF4 import Dataset
import numpy
import os
import shutil
import sys

# parse inputs
parser = ArgumentParser("GGCMI Aggregator Merger")
parser.add_argument("--crop", required=True, help="Crop short name")
parser.add_argument("--indir", required=True, help="Input directory")
parser.add_argument("--model", required=True, help="Model name")
parser.add_argument("--output", required=True, help="Output file")
parser.add_argument("--adaptation", required=True, help="Adaptation level (A0 or A1)")
parser.add_argument("--carbon_levels", default="C360,C510,C660,C810", help="Comma separated list of carbon levels")
parser.add_argument("--temp_levels", default="T-1,T0,T1,T2,T3,T4,T6", help="Comma separated list of temperature levels")
parser.add_argument("--water_levels", default="W-50,W-30,W-20,W-10,W0,W10,W20,W30,Winf", help="Comma separated list of water levels")
parser.add_argument("--nitrogen_levels", default="N10,N60,N200,NNA", help="Comma separated list of nitrogen levels")
args = parser.parse_args()

adaptation = args.adaptation
crop = args.crop
indir = args.indir
model = args.model
output = args.output

carbon_levels = args.carbon_levels.split(',')
temp_levels = args.temp_levels.split(',')
water_levels = args.water_levels.split(',')
nitrogen_levels = args.nitrogen_levels.split(',')

fill_value = 1.e+20
files = glob(os.path.join(args.indir, "%s_%s*%s.nc4" % (model, crop, adaptation)))
if not files:
    print "No files to process"
    sys.exit(0)

dim_names = list()
dim_sizes = list()
variables = list()

# Gather dimension and variable info from a sample file
with Dataset(files[0]) as sample:
    for dimension in sample.dimensions:
        dim_names.append(dimension)
        dim_sizes.append(sample.dimensions[dimension].size)
    for variable in sample.variables:
        if variable not in ["global", "gadm0", "time", "scen", "irr"]:
            variables.append(str(variable))

# Use sample data as output
shutil.copyfile(files[0], output)
variables_to_delete = list(set(variables) - set(dim_names))
new_variables = list(set(variables) - set(variables_to_delete))
os.system("ncks -h -O -x -v %s %s %s.tmp" % (",".join(variables_to_delete), output, output))
shutil.move("%s.tmp" % output, output)

# gadm0/global, time, scen, irr, c, t, w, n
data = dict()
for v in variables:
    data[v] = numpy.full((dim_sizes[0], dim_sizes[1], dim_sizes[2], dim_sizes[3], len(carbon_levels), len(temp_levels),
                          len(water_levels), len(nitrogen_levels)), fill_value)

for filename in files:
    print "Reading file %s" % filename
    bname = os.path.basename(filename)
    splits = bname.split('_')
    try:
        cidx = carbon_levels.index(splits[2])
        tidx = temp_levels.index(splits[3])
        widx = water_levels.index(splits[4])
        nidx = nitrogen_levels.index(splits[5])
    except ValueError:
        print "Ignoring file %s. Incorrectly formatted file, or invalid c/t/w/n level" % bname
        continue
    nc = Dataset(filename)
    for v in variables:
        try:
            data[v][:, :, :, :, cidx, tidx, widx, nidx] = nc.variables[v][:]
        except KeyError:
            print "Warning: File %s is missing variable %s" % (filename, v)

# Save output
with Dataset(output, "a") as o:
    o.createDimension("c", size=len(carbon_levels))
    o.createDimension("t", size=len(temp_levels))
    o.createDimension("w", size=len(water_levels))
    o.createDimension("n", size=len(nitrogen_levels))

    cvar = o.createVariable("c", "i4", ("c",))
    cvar.units = "mapping"
    cvar.long_name = ", ".join(carbon_levels)
    cvar[:] = numpy.arange(len(carbon_levels))

    tvar = o.createVariable("t", "i4", ("t",))
    tvar.units = "mapping"
    tvar.long_name = ", ".join(temp_levels)
    tvar[:] = numpy.arange(len(temp_levels))

    wvar = o.createVariable("w", "i4", ("w",))
    wvar.units = "mapping"
    wvar.long_name = ", ".join(water_levels)
    wvar[:] = numpy.arange(len(water_levels))

    nvar = o.createVariable("n", "i4", ("n",))
    nvar.units = "mapping"
    nvar.long_name = ", ".join(nitrogen_levels)
    nvar[:] = numpy.arange(len(nitrogen_levels))

    for v in variables:
        newvar = o.createVariable(v, "f4", (dim_names[0], dim_names[1], dim_names[2], dim_names[3],
                                  "c", "t", "w", "n"), zlib=True, complevel=9, fill_value=fill_value)
        newvar[:] = data[v][:]

# Remove intermediate files
for f in files:
    os.remove(f)
