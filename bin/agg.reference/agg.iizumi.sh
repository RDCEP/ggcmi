#!/bin/bash

# directory to reference data
refdir=/project/joshuaelliott/ggcmi/reference/iizumi/30min

# directory to weight files
wtsdir=/project/ggcmi/AgMIP.output/processed/masks/weight

# directory to aggregation files
mskdir=/project/ggcmi/AgMIP.output/processed/masks/aggr

# directory to save output
outdir=/project/joshuaelliott/ggcmi/bin/agg.reference

# crops to process
cpshort=(mai ric soy whe)
cpfull=(maize rice soy wheat)
cpref=(maize_major rice_major soybean wheat)

# masks to process
msks=(fpu kg global)

for ((i = 0; i < ${#msks[@]}; i++)); do # masks
   # mask name
   msk=${msks[$i]}

   # filename
   outfile=$outdir/iizumi.1982-2006.${msk}.nc4

   for ((j = 0; j < ${#cpshort[@]}; j++)); do # crops
      # crop short name
      cs=${cpshort[$j]}

      # crop full name
      cf=${cpfull[$j]}

      # crop reference name
      cr=${cpref[$j]}

      # aggregate
      echo Aggregating crop $cf to $msk level . . .
      ./agg.single.py -i $refdir/iizumi.2013JAN29.${cr}.1982-2006.30min.nc4:yield50 \
                      -a $mskdir/${msk}.mask.nc4:$msk                               \
                      -t mean                                                       \
                      -w $wtsdir/${cf}.nc4:sum                                      \
                      -o temp.nc4

      # clean data
      ncrename -O -h -v ${msk}_index,$msk temp.nc4 temp.nc4
      ncrename -O -h -v yield50_${msk},yield_${cs} temp.nc4 temp.nc4
      ncap2 -O -h -s "time[time]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}" temp.nc4 temp.nc4
      ncatted -O -h -a units,time,c,c,"years since 1982-01-01" temp.nc4 temp.nc4

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
