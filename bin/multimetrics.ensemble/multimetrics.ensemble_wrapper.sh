#!/bin/bash

infile=$1
reffile=$2
agglvl=$3
outdir=$4

metrics=(rmse)
munits=("t ha-1 yr-1")
mlongnames=("root mean squared error")
outfile=$outdir/$(basename $infile)
outfile=${outfile/ensemble/multimetrics}

for ((i = 0; i < ${#metrics[@]}; i++)); do
   multimetrics.ensemble.py -i $infile -r $reffile -a $agglvl -m "${metrics[$i]}" -u "${munits[$i]}" -l "${mlongnames[$i]}" -o tmp.nc
   if [ $i = 0 ]; then
      mv tmp.nc $outfile
   else
      ncks -h -A tmp.nc $outfile
      rm tmp.nc
   fi
done

exit 0
