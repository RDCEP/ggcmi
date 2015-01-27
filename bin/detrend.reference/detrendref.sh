#!/bin/bash

indir=/project/joshuaelliott/ggcmi/bin/agg.reference
outdir=/project/ggcmi/AgMIP.input/other.inputs/reference

PATH=$PATH:/project/joshuaelliott/ggcmi/utils

for r in ray iizumi; do
   if [ $r = ray ]; then
      yrs=1961-2008
   else
      yrs=1982-2006
   fi
   for a in gadm0 fpu kg global; do
      for w in fixed dynamic; do
         if [ $w = fixed ]; then
            ifile=$indir/${r}.${yrs}.${a}.fixed.nc4
            ofile=$outdir/$r/${r}.${yrs}.${a}.fixed.nc4
         else
            ifile=$indir/${r}.${yrs}.${a}.${r}.nc4
            ofile=$outdir/$r/${r}.${yrs}.${a}.${r}.nc4
         fi
         ./detrendref.py -i $ifile -a $a -c mai,ric,soy,whe -o $ofile
      done
   done
done
