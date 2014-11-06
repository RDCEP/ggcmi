#!/bin/bash

infile=$1
reffile=$2
agglvl=$3
outdir=$4

metrics=(tscorr varratio rmse hitrate rmse_extreme bias_extreme)
munits=("" "" "t ha-1 yr-1" "" "t ha-1 yr-1" "t ha-1 yr-1")
mlongnames=("time series correlation" "variance ratio" "root mean squared error" "hit rate" "RMSE of extreme events" "bias of extreme events")

outfile=$outdir/$(basename $infile)
outfile=${outfile/ensemble/multimetrics}

for ((i = 0; i < ${#metrics[@]}; i++)); do
   /project/joshuaelliott/ggcmi/bin/multimetrics.ensemble/multimetrics.ensemble.py -i $infile -r $reffile -a $agglvl -m "${metrics[$i]}" -u "${munits[$i]}" -l "${mlongnames[$i]}" -o tmp.nc
   if [ i = 0 ]; then
      mv tmp.nc $outfile
   else
      ncks -h -A tmp.nc $outfile
      rm tmp.nc
   fi
done

exit 0
