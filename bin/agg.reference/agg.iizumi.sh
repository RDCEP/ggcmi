#!/bin/bash

PATH=$PATH:/project/joshuaelliott/ggcmi/utils

# directory to reference data
refdir=/project/joshuaelliott/ggcmi/reference/iizumi/30min

# directory to weight files
wtsdir=/project/ggcmi/AgMIP.output/processed/masks/weight
wtsdirray=/project/joshuaelliott/ggcmi/reference/ray/masks

# directory to aggregation files
mskdir=/project/ggcmi/AgMIP.output/processed/masks/aggr

# directory to save output
outdir=/project/ggcmi/AgMIP.input/other.inputs/reference/iizumi/aggs

# crops to process
cpshort=(mai ric soy whe)
cpfull=(maize rice soy wheat)
cpref=(maize_major rice_major soybean wheat)

# aggregation masks to process
amsks=(gadm0 fpu kg global)

# weight masks to process
wmsks=(mirca iizumi ray spam)

for ((i = 0; i < ${#amsks[@]}; i++)); do # aggregation masks
   # aggregation mask name
   amsk=${amsks[$i]}

   for ((j = 0; j < ${#wmsks[@]}; j++)); do # weight masks
      # weight mask name
      wmsk=${wmsks[$j]}

      # filename
      if [ $wmsk = mirca ] || [ $wmsk = iizumi ] || [ $wmsk = spam ]; then
         outfile=$outdir/iizumi.1982-2006.${amsk}.fixed_${wmsk}_mask.nc4
      else
         outfile=$outdir/iizumi.1982-2006.${amsk}.dynamic_${wmsk}_mask.nc4
      fi

      for ((k = 0; k < ${#cpshort[@]}; k++)); do # crops
         # crop short name
         cs=${cpshort[$k]}

         # crop full name
         cf=${cpfull[$k]}

         # crop reference name
         cr=${cpref[$k]}

         # weight mask file
         if [ $wmsk = mirca ]; then
            wfile=$wtsdir/${cf}.nc4
         elif [ $wmsk = iizumi ]; then
            wfile=$wtsdir/${cf}.iizumi.nc4
         elif [ $wmsk = ray ]; then
            wfile=$wtsdirray/${cf}.ray.1982-2006.nc4
         else
            wfile=$wtsdir/${cf}.spam.nc4
         fi

         # aggregate
         echo Aggregating crop $cf to $amsk level with $wmsk weights . . .
         ./agg.single.py -i $refdir/iizumi.2013JAN29.${cr}.1982-2006.30min.nc4:yield50 \
                         -a $mskdir/${amsk}.mask.nc4:$amsk                             \
                         -t mean                                                       \
                         -w $wfile:sum                                                 \
                         -o temp.nc4

         # clean data
         ncrename -O -h -v ${amsk}_index,$amsk temp.nc4 temp.nc4
         ncrename -O -h -v yield50_${amsk},yield_${cs} temp.nc4 temp.nc4
         ncap2 -O -h -s "time[time]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}" temp.nc4 temp.nc4
         ncatted -O -h -a units,time,c,c,"years since 1982-01-01" temp.nc4 temp.nc4

         if [ $k = 0 ]; then
            mv temp.nc4 $outfile
         else
            ncks -h -A temp.nc4 $outfile
            rm temp.nc4
         fi
      done

      nccopy -d9 -k4 $outfile $outfile.2
      mv $outfile.2 $outfile
   done
done
