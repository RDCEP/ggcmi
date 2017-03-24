#!/bin/bash

indir=/project/ggcmi/phase1.final/processed.simple/biascorr
outdir=/project/ggcmi/phase1.final/processed.simple/multimetrics

aggmasks=(gadm0 global)
areas=(fixed_mirca_mask)
refshort=(faostat)
reflong=(faostat.1961-2012)

# Header
echo inputfile reffile agglvl outdir

for aggmask in ${aggmasks[@]}; do
   for ((ref = 0; ref < ${#refshort[@]}; ref++)); do
      rs=${refshort[$ref]}
      rl=${reflong[$ref]}
      if [ $rs = ray ]; then
         refdir=/project/joshuaelliott/ggcmi/reference
      else
         refdir=/project/ggcmi/AgMIP.input/other.inputs/reference
      fi
      for area in ${areas[@]}; do
         reffile=$refdir/$rs/$rl.$aggmask.$area.nc4
         if [ ! -f $reffile ] || [ ! -d $indir/$aggmask/$rs/$area ]; then
            continue
         fi
         mkdir -p $outdir/$aggmask/$rs/$area
         for file in $indir/$aggmask/$rs/$area/*; do
            echo $file $reffile $aggmask $outdir/$aggmask/$rs/$area
         done
      done
   done
done
