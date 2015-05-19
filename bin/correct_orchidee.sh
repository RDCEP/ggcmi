#!/bin/bash

dir=/project/ggcmi/AgMIP.output/ORCHIDEE/WFDEI.GPCC

for c in maize wheat sunflower; do
   for f in `ls $dir/$c/*`; do
      ncpdq -O -h -a lat,lon $f $f
   done
done
