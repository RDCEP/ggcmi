#!/bin/bash

indir=/project/ggcmi/AgMIP.input/other.inputs/AGMIP_GROWING_SEASON.HARM.version1.23
outdir=/project/ggcmi/AgMIP.input/other.inputs/agmip_growing_season.1.23.processed

clong=(Maize Wheat Rice Soybeans Sorghum Millet)
cshort=(mai whe ric soy sor mil)

for ((i = 0; i < ${#clong[@]}; i++)); do
   cl=${clong[$i]}
   cs=${cshort[$i]}
   echo Processing $cl . . .

   # irrigated
   ncecat -O -h -u irr $indir/${cl}_ir_growing_season_dates_v1.23.nc4 ir.nc4
   ncap2 -O -h -s "irr[irr]=0" ir.nc4 ir.nc4

   # rainfed
   ncecat -O -h -u irr $indir/${cl}_rf_growing_season_dates_v1.23.nc4 rf.nc4
   ncap2 -O -h -s "irr[irr]=1" rf.nc4 rf.nc4

   # combined
   ncrcat -h ir.nc4 rf.nc4 combined.nc4
   ncrename -O -h -v "harvest day",harvest_day combined.nc4 combined.nc4
   ncrename -O -h -v "planting day",planting_day combined.nc4 combined.nc4
   ncrename -O -h -v "growing season length",growing_season_length combined.nc4 combined.nc4
   ncatted -O -h -a long_name,irr,c,c,"ir, rf" combined.nc4 combined.nc4
   ncpdq -O -h -a -lat combined.nc4 combined.nc4
   nccopy -d9 -k4 combined.nc4 $outdir/${cs}_growing_season_dates.nc4

   rm ir.nc4 rf.nc4 combined.nc4
done
