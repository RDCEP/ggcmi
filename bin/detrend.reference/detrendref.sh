#!/bin/bash

PATH=$PATH:/project/joshuaelliott/ggcmi/utils

indir=/project/ggcmi/AgMIP.input/other.inputs/reference

for r in ray iizumi; do
   if [ $r = ray ]; then
      yrs=1961-2008
   else
      yrs=1982-2006
   fi
   for a in gadm0 fpu kg global; do
      for w in mirca ray iizumi spam; do
         if [ $w = mirca ] || [ $w = iizumi ] || [ $w = spam ]; then
            ifile=$indir/$r/aggs/${r}.${yrs}.${a}.fixed_${w}_mask.nc4
            ofile=$indir/$r/${r}.${yrs}.${a}.fixed_${w}_mask.nc4
         else
            ifile=$indir/$r/aggs/${r}.${yrs}.${a}.dynamic_${w}_mask.nc4
            ofile=$indir/$r/${r}.${yrs}.${a}.dynamic_${w}_mask.nc4
         fi
         ./detrendref.py -i $ifile -a $a -c mai,ric,soy,whe -o $ofile
      done
   done
done
