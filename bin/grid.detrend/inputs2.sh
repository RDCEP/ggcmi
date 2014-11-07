#!/bin/bash

root=/project/joshuaelliott/ggcmi/reference

# Header
echo indir infile

# iizumi
for f in $root/iizumi/30min/*.nc4; do
   echo $root/iizumi/detrended $(basename $f)
done

# ray
for v in mai ric soy whe; do
   echo $root/ray/detrended ${v}_weight_ray_1961-2008.nc4
done
