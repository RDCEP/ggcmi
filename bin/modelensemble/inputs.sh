#!/bin/bash

bcdir=/project/ggcmi/AgMIP.output/processed/biascorr
mmdir=/project/ggcmi/AgMIP.output/processed/multimetrics
outdir=/project/ggcmi/AgMIP.output/processed/modelensemble

# Header
echo indir metricsdir agglvl outdir

for a in gadm0 fpu kg global; do
   for ref in faostat ray iizumi; do   
      for area in fixed ray iizumi; do
         if [ ! -d $bcdir/$a/$ref/$area ] || [ ! -d $mmdir/$a/$ref/$area ]; then
            continue
         fi
         mkdir -p $outdir/$a/$ref/$area
         echo $bcdir/$a/$ref/$area $mmdir/$a/$ref/$area $a $outdir/$a/$ref/$area
      done
   done
done
