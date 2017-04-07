#!/bin/bash

# Parameter file
params=$1
if [ -z "$params" ]; then
    echo "Usage: $0 <params>"
    exit 1
fi

COMMONDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../common" && pwd )"
source $COMMONDIR/common_inputs.sh

model_directory=$( get_param model_directory )
biascorr_directory=$( get_param biascorr_directory )
weight_directory=$( get_param weight_directory )
rescaled_directory=$( get_param rescaled_directory )
models=( $( get_param models ) )
climates=$( get_param climates )
crops_long=( $(get_param crops_long) )
crops_short=( $(get_param crops_short) )
area_long=$( area_to_long $( get_param rescaler:area ))
reference=$( get_param rescaler:reference )
level=$( get_param rescaler:level )

# Header
echo irfile rffile bcfile wtfile vname outfile

for model in ${models[@]}; do
  for climate in $climates; do
    for ((i = 0; i < ${#crops_long[@]}; i++)); do
      for v in yield; do
        irfiles=($(ls ${model_directory}/$model/$climate/${crops_long[$i]}/*firr*$v*${crops_short[$i]}* 2>/dev/null))
        rffiles=($(ls ${model_directory}/$model/$climate/${crops_long[$i]}/*noirr*$v*${crops_short[$i]}* 2>/dev/null))
        bcfile=$(ls ${biascorr_directory}/$level/$reference/$area_long/${model,,}*_${climate,,}*${crops_short[$i]}* 2>/dev/null)
        wtfile=${weight_directory}/${crops_long[$i]}.nc4
        if [ ${#irfiles[@]} -ge 1 ] && [ $bcfile ]; then # irrigated and bias-corrected file exist
          for ((j = 0; j < ${#irfiles[@]}; j++)); do
            if [ ${#rffiles[@]} = 0 ]; then
              rffile=None # no rainfed file
            else
              rffile=${rffiles[$j]}
            fi
            mkdir -p ${rescaled_directory}/$level/$reference/$area_long
            outfile=${rescaled_directory}/$level/$reference/$area_long/$(basename ${irfiles[$j]//_firr})
            outfile=${outfile/nc4/rescaled.nc4}
            echo ${irfiles[$j]} $rffile $bcfile $wtfile ${v}_${crops_short[$i]} $outfile
          done
        fi
      done
    done
  done
done
