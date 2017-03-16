#!/bin/bash

# Directory structures
processed=/project/ggcmi/phase1.final/processed.simple
aggrdir=$processed/masks/aggr
growdir=/project/ggcmi/AgMIP.input/other.inputs/agmip_growing_season.1.23.processed
inputdir=/project/ggcmi/phase1.final
outroot=$processed/aggs
raydir=/project/joshuaelliott/ggcmi/reference/ray/masks
weightdir=$processed/masks/weight

# Processing options
aggmasks=(gadm0 fpu kg global)
areas=(mirca)
cropslong=(maize wheat rice soy sorghum millet)
cropsshort=(mai whe ric soy sor mil)
models=(pDSSAT pAPSIM LPJ-GUESS LPJmL PEGASUS GEPIC EPIC-IIASA EPIC-Boku EPIC-test CGMS-WOFOST CLM-Crop EPIC-TAMU ORCHIDEE ORCHIDEE-crop PEPIC PRYSBI2)
weaths=(AgMERRA AgCFSR CFSR ERAI GRASP WATCH WFDEI.CRU WFDEI.GPCC Princeton)
wyears=(1980_2010 1980_2010 1980_2010 1979_2010 1961_2010 1958_2001 1979_2012 1979_2010 1948_2008)

# Header
echo indir crop lufile agg gsfile outfile

for model in ${models[@]}; do
   for ((w = 0; w < ${#weaths[@]}; w++)); do
      weath=${weaths[$w]}
      wyear=${wyears[$w]}
      for ((c = 0; c < ${#cropslong[@]}; c++)); do
         croplong=${cropslong[$c]}
         cropshort=${cropsshort[$c]}
         indir=$inputdir/$model/$weath/$croplong
         if [ ! -d $indir ] || [ $(ls $indir | wc -l) = 0 ]; then
            continue
         fi
         gsfile=$growdir/${cropshort}_growing_season_dates.nc4
         for aggmask in ${aggmasks[@]}; do
            aggmaskfile=$aggrdir/$aggmask.mask.nc4
            for area in ${areas[@]}; do
               if [ $area = mirca ]; then
                  wfile=$weightdir/$croplong.nc4
               elif [ $area = ray ]; then
                  wfile=$raydir/$croplong.$area.nc4
               else
                  wfile=$weightdir/$croplong.$area.nc4
               fi
               if [ ! -f $wfile ]; then
                  continue
               fi
               if [ $area = mirca ] || [ $area = iizumi ] || [ $area = spam ]; then
                  outdir=$outroot/$aggmask/fixed_${area}_mask
               else
                  outdir=$outroot/$aggmask/dynamic_${area}_mask
               fi
               ofile=$outdir/${model,,}_${weath,,}_hist_${cropshort}_annual_${wyear}.nc4
               mkdir -p $outdir
               echo $indir $cropshort $wfile $aggmaskfile:$aggmask $gsfile $ofile
            done
         done
      done
   done
done
