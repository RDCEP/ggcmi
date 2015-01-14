#!/bin/bash

indir=/project/ggcmi/AgMIP.output/processed/masks/weight
mskdir=/project/ggcmi/AgMIP.output/processed/masks/aggr
outdir=/project/ggcmi/AgMIP.output/processed/masks/weight/aggs

for c in maize millet rice sorghum soy wheat; do
   if=$indir/$c.nc4
   cp $if .
   ncecat -O -h -u time $c.nc4 $c.nc4
   ncap2 -O -h -s "time[time]=1" $c.nc4 $c.nc4
   for a in fpu kg; do
      ./agg.single.py -i $c.nc4:rainfed,irrigated -a $mskdir/${a}.mask.nc4:$a -t sum -o temp.nc4
      ncwa -O -h -a time temp.nc4 temp.nc4
      ncks -O -h -x -v time temp.nc4 temp.nc4
      nccopy -d9 -k4 temp.nc4 $outdir/${c}.${a}.nc4   
      rm temp.nc4
   done
   rm $c.nc4
done
