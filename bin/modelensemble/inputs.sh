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
modelensemble_directory=$( get_param modelensemble_directory )
aggregation_levels=$( get_param aggregation_levels )
reference_short=$( get_param reference_short )
areas=$( get_param areas )

# Header
echo indir metricsdir agglvl outdir

for aggmask in $aggregation_levels; do
   for ref in $reference_short; do
      for area in $areas; do
         area=$( area_to_long $area )
         if [ ! -d ${biascorr_directory}/$aggmask/$ref/$area ] || [ ! -d ${multimetrics_directory}/$aggmask/$ref/$area ]; then
            continue
         fi
         mkdir -p ${modelensemble_directory}/$aggmask/$ref/$area
         echo ${biascorr_directory}/$aggmask/$ref/$area ${multimetrics_directory}/$aggmask/$ref/$area $aggmask ${modelensemble_directory}/$aggmask/$ref/$area
      done
   done
done
