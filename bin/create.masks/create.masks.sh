#!/bin/bash

raydir=/project/joshuaelliott/ggcmi/reference/ray/weighted
iizumidir=/project/joshuaelliott/ggcmi/reference/iizumi/30min
wtsdir=/project/ggcmi/AgMIP.output/processed/masks/weight
iizumioutdir=/project/ggcmi/AgMIP.output/processed/masks/weight
rayoutdir=/project/joshuaelliott/ggcmi/reference/ray/masks

wtsnames=(maize rice soy wheat)

# ray
refnames=(mai ric soy whe)
for ((i = 0; i < ${#refnames[@]}; i++)); do
   wn=${wtsnames[$i]}
   rn=${refnames[$i]}
   ./create.masks.py -i $raydir/${rn}_weight_ray_1961-2008.nc4 -w $wtsdir/$wn.nc4 -v area_${rn} -y 1961 -o $rayoutdir/${wn}.ray.nc4
done

# iizumi
refnames=(maize_major rice_major soybean wheat)
for ((i = 0; i < ${#refnames[@]}; i++)); do
   wn=${wtsnames[$i]}
   rn=${refnames[$i]}
   ./create.masks.py -i $iizumidir/iizumi.2013JAN29.${rn}.1982-2006.30min.nc4 -w $wtsdir/$wn.nc4 -v area -y 1982 -o $iizumioutdir/${wn}.iizumi.nc4
done
