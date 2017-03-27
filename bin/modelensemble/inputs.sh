#!/bin/bash

bcdir=/project/ggcmi/phase1.final/processed.simple/biascorr
mmdir=/project/ggcmi/phase1.final/processed.simple/multimetrics
outdir=/project/ggcmi/phase1.final/processed.simple/modelensemble
aggmasks="gadm0 global"
refs="faostat"
areas="fixed_mirca_mask"

# Header
echo indir metricsdir agglvl outdir

for aggmask in $aggmasks; do
   for ref in $refs; do
      for area in $areas; do
         if [ ! -d $bcdir/$aggmask/$ref/$area ] || [ ! -d $mmdir/$aggmask/$ref/$area ]; then
            continue
         fi
         mkdir -p $outdir/$aggmask/$ref/$area
         echo $bcdir/$aggmask/$ref/$area $mmdir/$aggmask/$ref/$area $aggmask $outdir/$aggmask/$ref/$area
      done
   done
done
