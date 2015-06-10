#!/bin/bash

PATH=$PATH:/project/joshuaelliott/ggcmi/utils

# directory to reference data
refdir=/project/joshuaelliott/ggcmi/reference/ray/weighted

# directory to weight files
wtsdir=/project/ggcmi/AgMIP.output/processed/masks/weight
wtsdirray=/project/joshuaelliott/ggcmi/reference/ray/masks

# directory to aggregation files
mskdir=/project/ggcmi/AgMIP.output/processed/masks/aggr

# directory to save output
outdir=/project/joshuaelliott/ggcmi/reference/ray/aggs

# crops to process
cpshort=(mai ric soy whe)
cpfull=(maize rice soy wheat)

# aggregation masks to process
amsks=(gadm0 fpu kg global)

# weight masks to process
wmsks=(mirca ray iizumi spam)

for ((i = 0; i < ${#amsks[@]}; i++)); do # aggregation masks
   # aggregation mask name
   amsk=${amsks[$i]}

   for ((j = 0; j < ${#wmsks[@]}; j++)); do # weight masks
      # weight mask name
      wmsk=${wmsks[$j]}

      # filename
      if [ $wmsk = mirca ] || [ $wmsk = iizumi ] || [ $wmsk = spam ]; then
         outfile=$outdir/ray.1961-2008.${amsk}.fixed_${wmsk}_mask.nc4
      else
         outfile=$outdir/ray.1961-2008.${amsk}.dynamic_${wmsk}_mask.nc4
      fi

      for ((k = 0; k < ${#cpshort[@]}; k++)); do # crops
         # crop short name
         cs=${cpshort[$k]}

         # crop full name
         cf=${cpfull[$k]}

         # weight mask file
         if [ $wmsk = mirca ]; then
            wfile=$wtsdir/${cf}.nc4
            infile=${cs}_weight_ray_1961-2008.nc4
         elif [ $wmsk = ray ]; then
            wfile=$wtsdirray/${cf}.ray.nc4
            infile=${cs}_weight_ray_1961-2008.nc4
         elif [ $wmsk = iizumi ]; then
            wfile=$wtsdir/${cf}.iizumi.nc4
            infile=${cs}_weight_ray_1982-2006.nc4
         else
            wfile=$wtsdir/${cf}.spam.nc4
            infile=${cs}_weight_ray_1961-2008.nc4
         fi

         # aggregate
         echo Aggregating crop $cf to $amsk level with $wmsk weights . . .
         ./agg.single.py -i $refdir/$infile:yield_${cs}    \
                         -a $mskdir/${amsk}.mask.nc4:$amsk \
                         -t mean                           \
                         -w $wfile:sum                     \
                         -o temp.nc4

         # clean data
         ncrename -O -h -v ${amsk}_index,$amsk temp.nc4 temp.nc4
         ncrename -O -h -v yield_${cs}_${amsk},yield_${cs} temp.nc4 temp.nc4

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
