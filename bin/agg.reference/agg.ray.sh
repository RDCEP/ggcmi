#!/bin/bash

# directory to reference data
refdir=/project/joshuaelliott/ggcmi/reference/ray/weighted

# directory to weight files
wtsdir=/project/ggcmi/AgMIP.output/processed/masks/weight

# directory to aggregation files
mskdir=/project/ggcmi/AgMIP.output/processed/masks/aggr

# directory to save output
outdir=/project/joshuaelliott/ggcmi/bin/agg.reference

# crops to process
cpshort=(mai ric soy whe)
cpfull=(maize rice soy wheat)

# masks to process
msks=(fpu kg global)

for ((i = 0; i < ${#msks[@]}; i++)); do # masks
   # mask name
   msk=${msks[$i]}

   # filename
   outfile=$outdir/ray.1961-2008.${msk}.nc4

   for ((j = 0; j < ${#cpshort[@]}; j++)); do # crops
      # crop short name
      cs=${cpshort[$j]}

      # crop full name
      cf=${cpfull[$j]}

      # aggregate
      echo Aggregating crop $cf to $msk level . . .
      ./agg.single.py -i $refdir/${cs}_weight_ray_1961-2008.nc4:yield_${cs} \
                      -a $mskdir/${msk}.mask.nc4:$msk                       \
                      -t mean                                               \
                      -w $wtsdir/${cf}.nc4:sum                              \
                      -o temp.nc4

      # clean data
      ncrename -O -h -v ${msk}_index,$msk temp.nc4 temp.nc4
      ncrename -O -h -v yield_${cs}_${msk},yield_${cs} temp.nc4 temp.nc4
      ncap2 -O -h -s "time=time-1" temp.nc4 temp.nc4
      ncatted -O -h -a units,time,m,c,"years since 1961-01-01" temp.nc4 temp.nc4

      if [ $j = 0 ]; then
         mv temp.nc4 $outfile
      else
         ncks -h -A temp.nc4 $outfile
         rm temp.nc4
      fi
   done

   nccopy -d9 -k4 $outfile $outfile.2
   mv $outfile.2 $outfile
done
