#!/bin/bash

# Parameter file
params=$1
if [ -z "$params" ]; then
    echo "Usage: $0 <params>"
    exit 1
fi

# Import common functions
COMMONDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../common" && pwd )"
source $COMMONDIR/common_inputs.sh

biascorr_directory=$( get_param biascorr_directory )
multimetrics_directory=$( get_param multimetrics_directory )
aggregation_levels=$( get_param aggregation_levels )
areas=( $( get_param areas ) )
reference_short=( $( get_param reference_short ) )
reference_long=( $( get_param reference_long ) )

# Header
echo inputfile reffile agglvl outdir

for aggmask in ${aggregation_levels[@]}; do
   for ((ref = 0; ref < ${#reference_short[@]}; ref++)); do
      rs=${reference_short[$ref]}
      rl=${reference_long[$ref]}
      if [ $rs = ray ]; then
         refdir=/project/joshuaelliott/ggcmi/reference
      else
         refdir=/project/ggcmi/AgMIP.input/other.inputs/reference
      fi
      for area in ${areas[@]}; do
         area=$( area_to_long $area )

         reffile=$refdir/$rs/$rl.$aggmask.${area}.nc4
         if [ ! -f $reffile ] || [ ! -d ${biascorr_directory}/$aggmask/$rs/${area} ]; then
            continue
         fi
         mkdir -p ${multimetrics_directory}/$aggmask/$rs/${area}
         for file in ${biascorr_directory}/$aggmask/$rs/${area}/*; do
            echo $file $reffile $aggmask ${multimetrics_directory}/$aggmask/$rs/${area}
         done
      done
   done
done
