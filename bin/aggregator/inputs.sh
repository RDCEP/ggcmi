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
aggr_mask_directory=$( get_param aggr_mask_directory )
growing_season_directory=$( get_param growing_season_directory )
model_directory=$( get_param model_directory )
agg_directory=$( get_param agg_directory )
ray_mask_directory=$( get_param ray_mask_directory )
weight_directory=$( get_param weight_directory )
aggregation_levels=( $( get_param aggregation_levels ) )
areas=( $( get_param areas) )
crops_long=( $( get_param crops_long ) )
crops_short=( $( get_param crops_short ) )
models=( $( get_param models ) )
climates=( $( get_param climates ) )
climate_years=( $( get_param climate_years ) )
variables=$( get_param variables )
adaptation_levels=$( get_param adaptation_levels )
co2_levels=$( get_param co2_levels )
temperature_levels=$( get_param temperature_levels )
precipitation_levels=$( get_param precipitation_levels )
nitrogen_levels=$( get_param nitrogen_levels )

# Header
echo indir crop lufile agg gsfile co2 temperature precip nitrogen adaptation outfile

for model in ${models[@]}; do
    indir=$model_directory/$model/phase2
    if [ ! -d $indir ] || [ $(ls $indir | wc -l) = 0 ]; then
        continue
    fi
    for ((c = 0; c < ${#crops_long[@]}; c++)); do
        croplong=${crops_long[$c]}
        cropshort=${crops_short[$c]}
        gsfile=$growing_season_directory/${cropshort}_growing_season_dates.nc4
        for variable in $variables; do
            for aggmask in ${aggregation_levels[@]}; do
                aggmaskfile=$aggr_mask_directory/$aggmask.mask.nc4
                for area in ${areas[@]}; do
                    if [ $area = mirca ]; then
                        wfile=$weight_directory/$croplong.nc4
                    elif [ $area = ray ]; then
                        wfile=$ray_directory/$croplong.$area.nc4
                    else
                        wfile=$weight_directory/$croplong.$area.nc4
                    fi
                    if [ ! -f $wfile ]; then
                        continue
                    fi
                    if [ $area = mirca ] || [ $area = iizumi ] || [ $area = spam ]; then
                        outdir=$agg_directory/$aggmask/fixed_${area}_mask
                    else
                        outdir=$agg_directory/$aggmask/dynamic_${area}_mask
                    fi
                    mkdir -p $outdir
	            for co2 in $co2_levels; do
                        for temperature in $temperature_levels; do
                            for precip in $precipitation_levels; do
                                for nitrogen in $nitrogen_levels; do
                                    for adaptation in $adaptation_levels; do
                                        ofile=$outdir/${model,,}_${cropshort}_${co2}_${temperature}_${precip}_${nitrogen}_${adaptation}.nc4
                                        echo $indir $cropshort $wfile $aggmaskfile:$aggmask $gsfile $co2 $temperature $precip $nitrogen $adaptation $ofile
                                    done
                                done
                            done
                        done
                    done
                done
            done
        done
    done
done
