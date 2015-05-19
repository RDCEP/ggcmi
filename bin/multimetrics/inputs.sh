#!/bin/bash

indir=/project/ggcmi/AgMIP.output/processed/biascorr
outdir=/project/ggcmi/AgMIP.output/processed/multimetrics

refs=(faostat ray iizumi)
refl=(faostat.1961-2012 ray.1961-2008 iizumi.1982-2006)

# Header
echo indir reffile agglvl outdir

for a in gadm0 fpu kg global; do
   for ((ref = 0; ref < ${#refs[@]}; ref++)); do
      rs=${refs[$ref]}
      rl=${refl[$ref]}
      if [ $rs = ray ]; then
         refdir=/project/joshuaelliott/ggcmi/reference
      else
         refdir=/project/ggcmi/AgMIP.input/other.inputs/reference
      fi
      for area in fixed_mirca_mask dynamic_ray_mask fixed_iizumi_mask fixed_spam_mask; do
         reffile=$refdir/$rs/$rl.$a.$area.nc4
         if [ ! -f $reffile ] || [ ! -d $indir/$a/$rs/$area ]; then
            continue
         fi
         mkdir -p $outdir/$a/$rs/$area
         echo $indir/$a/$rs/$area $reffile $a $outdir/$a/$rs/$area
      done
   done
done
