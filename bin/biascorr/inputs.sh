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

agg_directory=$( get_param agg_directory )
biascorr_directory=$( get_param biascorr_directory )
areas=( $( get_param areas ) )
aggregation_levels=( $( get_param aggregation_levels ) )
reference_short=( $( get_param reference_short ) )
reference_long=( $( get_param reference_long ) )
reference_directory=$( get_param reference_directory )
ray_reference_directory=$( get_param ray_reference_directory )

# Header
echo inputfile reffile agglvl outdir

for aggmask in ${aggregation_levels[@]}; do
    for ((ref = 0; ref < ${#reference_short[@]}; ref++)); do
        rs=${reference_short[$ref]}
        rl=${reference_long[$ref]}
        if [ $rs = ray ]; then
            refdir=$ray_reference_directory
        else
            refdir=$reference_directory
        fi

        for area in ${areas[@]}; do
            area=$( area_to_long $area )
            reffile=$refdir/$rs/${rl}.${aggmask}.${area}.nc4

            if [ ! -f $reffile ] || [ ! -d $agg_directory/$aggmask/${area} ]; then
                echo Cant find $reffile or $agg_directory/$aggmask/${area} is not a directory
                continue
            fi

            mkdir -p $biascorr_directory/$aggmask/$rs/${area}
            for file in $agg_directory/$aggmask/${area}/*; do
                echo $file $reffile $aggmask $biascorr_directory/$aggmask/$rs/${area}
            done
        done
    done
done
