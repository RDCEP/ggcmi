#!/bin/bash

infiles=(maize.resamp.nc4 millet.resamp.nc4 rice.resamp.nc4 sorghum.resamp.nc4 soybean.resamp.nc4 wheat.resamp.nc4)

agfiles=(maize.agg.nc4 millet.agg.nc4 rice.agg.nc4 sorghum.agg.nc4 soybean.agg.nc4 wheat.agg.nc4)

for ((i = 0; i < ${#infiles[@]}; i++)); do
  if=${infiles[$i]}
  of=${agfiles[$i]}
  ncecat -O -h -u time $if $if
  ncap2 -O -h -s "time[time]=1" $if $if
  ./agg.single.py -i $if:irrigated,rainfed -a fpu.kg.mask.nc4:fpu,kg -t sum -o $of
  ncwa -O -h -a time $of $of
  ncks -O -h -x -v time $of $of
  nccopy -d9 -k4 $of $of.2
  mv $of.2 $of
done
