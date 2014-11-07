#!/bin/bash

iizumi=/project/joshuaelliott/ggcmi/reference/iizumi/detrended
ray=/project/joshuaelliott/ggcmi/reference/ray/detrended

for f in `ls /project/joshuaelliott/ggcmi/reference/iizumi/30min/*.nc4`; do
   fn=$(basename $f)
   echo $fn
   ncrcat -O -h $iizumi/$fn* $iizumi/$fn
   ncpdq -O -h -a time,lat $iizumi/$fn $iizumi/$fn
   ncpdq -O -h -a lat,lon $iizumi/$fn $iizumi/$fn
done

rm $iizumi/*.nc4.*
for d in {1..8}; do rm -r $iizumi/$d; done

for f in `ls /project/joshuaelliott/ggcmi/reference/ray/weighted/*.nc4`; do
   fn=$(basename $f)
   echo $fn
   ncrcat -O -h $ray/$fn* $ray/$fn
   ncpdq -O -h -a time,lat $ray/$fn $ray/$fn
   ncpdq -O -h -a lat,lon $ray/$fn $ray/$fn
done

rm $ray/*.nc4.*
for d in {1..8}; do rm -r $ray/$d; done
