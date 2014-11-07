#!/bin/bash

root=/project/joshuaelliott/ggcmi/reference

# Header
echo infile vrlist outdir

# iizumi
for f in $root/iizumi/30min/*.nc4; do
   echo $f yield05,yield50,yield95,statyield $root/iizumi/detrended
done

# ray
for v in mai ric soy whe; do
   echo $root/ray/weighted/${v}_weight_ray_1961-2008.nc4 yield_$v $root/ray/detrended
done
