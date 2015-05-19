#!/bin/bash

idir=/project/ggcmi/AgMIP.output
wdir=/project/ggcmi/AgMIP.output/processed/masks/weight
adir=/project/ggcmi/AgMIP.output/processed/masks/aggr
gdir=/project/ggcmi/AgMIP.input/other.inputs/agmip_growing_season.1.23.processed
odir=/project/ggcmi/AgMIP.output/processed/aggs

mods=(pDSSAT pAPSIM LPJ-GUESS LPJmL PEGASUS GEPIC EPIC-IIASA EPIC-Boku EPIC-test CGMS-WOFOST CLM-Crop EPIC-TAMU ORCHIDEE ORCHIDEE-crop PEPIC PRYSBI2)

cropl=(maize wheat rice soy sorghum millet)
crops=(mai whe ric soy sor mil)

weaths=(AgMERRA AgCFSR CFSR ERAI GRASP WATCH WFDEI.CRU WFDEI.GPCC Princeton)
wyears=(1980_2010 1980_2010 1980_2010 1979_2010 1961_2010 1958_2001 1979_2012 1979_2010 1948_2008)

# Header
echo indir crop lufile agg gsfile outfile

for m in ${mods[@]}; do
   for ((w = 0; w < ${#weaths[@]}; w++)); do
      weath=${weaths[$w]}
      wyear=${wyears[$w]}
      for ((c = 0; c < ${#cropl[@]}; c++)); do
         cl=${cropl[$c]}
         cs=${crops[$c]}
         indir=$idir/$m/$weath/$cl
         if [ ! -d $indir ] || [ $(ls $indir | wc -l) = 0 ]; then
            continue
         fi
         gsfile=$gdir/${cs}_growing_season_dates.nc4
         for a in gadm0 fpu kg global; do
            afile=$adir/$a.mask.nc4
            for area in mirca ray iizumi spam; do
               if [ $area = mirca ]; then
                  wfile=$wdir/$cl.nc4
               else
                  wfile=$wdir/$cl.$area.nc4
               fi
               if [ ! -f $wfile ]; then
                  continue
               fi
               if [ $area = mirca ] || [ $area = iizumi ] || [ $area = spam ]; then
                  outdir=$odir/$a/fixed_${area}_mask
               else
                  outdir=$odir/$a/dynamic_${area}_mask
               fi
               ofile=$outdir/${m,,}_${weath,,}_hist_${cs}_annual_${wyear}.nc4
               if [ ! -f $ofile ]; then
                  mkdir -p $outdir
                  echo $indir $cs $wfile $afile:$a $gsfile $ofile
               fi
            done
         done
      done
   done
done
