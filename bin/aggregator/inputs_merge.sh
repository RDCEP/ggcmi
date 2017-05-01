#!/bin/bash

# Parameter file
params=$1
if [ -z "$params" ]; then
    echo "Usage: $0 <params>"
    exit 1
fi

COMMONDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../common" && pwd )"
source $COMMONDIR/common_inputs.sh

# Read parameters
agg_directory=$( get_param agg_directory )
agg_levels=$( get_param aggregation_levels )
models=$( get_param models )
crops_short=$( get_param crops_short )
adaptation_levels=$( get_param adaptation_levels )
areas=$( get_param areas )
clevs=$( get_param co2_levels | sed s/' '/,/g )
tlevs=$( get_param temperature_levels | sed s/' '/,/g )
wlevs=$( get_param precipitation_levels | sed s/' '/,/g )
nlevs=$( get_param nitrogen_levels | sed s/' '/,/g )

# Header
echo indir model crop clevs tlevs wlevs nlevs adaptation output

for agg_level in $agg_levels; do
    for area in $areas; do
        area=$( area_to_long $area )
        for model in $models; do
            model=$( echo $model | tr '[:upper:]' '[:lower:]' )
            for crop in $crops_short; do
                for adaptation in $adaptation_levels; do
                    indir=$agg_directory/$agg_level/$area
                    if [ ! -d $indir ] || [ $( ls $indir | wc -l ) = 0 ]; then
                        echo Directory $indir empty or does not exist
                        continue
                    fi
		    output=$agg_directory/$agg_level/$area/${model}_${crop}_${adaptation}.nc4
                    echo $indir $model $crop $clevs $tlevs $wlevs $nlevs $adaptation $output
                done
            done
        done
    done
done
