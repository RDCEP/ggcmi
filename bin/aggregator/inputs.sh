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

# Header
echo indir crop lufile agg gsfile year_start outfile

for model in ${models[@]}; do
   for ((w = 0; w < ${#climates[@]}; w++)); do
      weath=${climates[$w]}
      wyear=${climate_years[$w]}
      start_year=$( echo $wyear | cut -d'_' -f1 )
      for ((c = 0; c < ${#crops_long[@]}; c++)); do
         croplong=${crops_long[$c]}
         cropshort=${crops_short[$c]}
         indir=$model_directory/$model/$weath/$croplong
         if [ ! -d $indir ] || [ $(ls $indir | wc -l) = 0 ]; then
            continue
         fi
         gsfile=$growing_season_directory/${cropshort}_growing_season_dates.nc4
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
               ofile=$outdir/${model,,}_${weath,,}_hist_${cropshort}_annual_${wyear}.nc4
               mkdir -p $outdir
               echo $indir $cropshort $wfile $aggmaskfile:$aggmask $gsfile $start_year $ofile
            done
         done
      done
   done
done
