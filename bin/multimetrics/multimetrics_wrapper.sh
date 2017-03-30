#!/bin/bash

infile=$1
reffile=$2
agglvl=$3
outdir=$4
params=$5

# Source common wrapper functions
COMMONDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../common" && pwd )"
source $COMMONDIR/common_inputs.sh

metrics=( $( get_param metrics ) )
metrics_units=( $( get_param metrics_units ) )
metrics_longnames=( $( get_param metrics_longnames ) )

outfile=$outdir/$(basename $infile)
outfile=${outfile/biascorr/multimetrics}

for ((i = 0; i < ${#metrics[@]}; i++)); do
   multimetrics.py -i $infile -r $reffile -a $agglvl -m "${metrics[$i]}" -u "${metrics_units[$i]}" -l "${metrics_longnames[$i]}" -o tmp.nc -p $params
   if [ $i = 0 ]; then
      mv tmp.nc $outfile
   else
      ncks -h -A tmp.nc $outfile
      rm tmp.nc
   fi
done

exit 0
