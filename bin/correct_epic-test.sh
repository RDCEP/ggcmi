#!/bin/bash

for f in /project/ggcmi/AgMIP.output/EPIC-test/WFDEI.GPCC/maize/*; do
   hasspinup=$(ncdump -h $f | grep "time = 96")
   if [ "$hasspinup" ]; then
      ncks -O -h -d time,64,94 $f $f
      ncap2 -O -h -s "time=time-64" $f $f
      ncatted -O -h -a units,time,m,c,"growing seasons since 1979-01-01" $f $f
   fi
done
