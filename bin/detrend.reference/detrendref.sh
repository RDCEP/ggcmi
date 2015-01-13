#!/bin/bash

indir=/project/joshuaelliott/ggcmi/bin/agg.reference
outdir=/project/ggcmi/AgMIP.input/other.inputs/reference

for r in ray iizumi; do
   if [ $r = ray ]; then
      yrs=1961-2008
   else
      yrs=1982-2006
   fi
   for a in fpu kg global; do
      ./detrendref.py -i $indir/${r}.${yrs}.${a}.nc4 -a $a -c mai,ric,soy,whe -o $outdir/$r/${r}.${yrs}.${a}.2.nc4
   done
done
