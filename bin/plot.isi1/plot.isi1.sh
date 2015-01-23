#!/bin/bash

plot=$1
crop=$2

mfile=/project/ggcmi/isi1.metrics.all/metrics.nc4
afile=/project/ggcmi/AgMIP.output/processed/masks/aggr/fpu.mask.nc4
sfile=/project/joshuaelliott/data/masks/fpu_shp_files/fpu_wgs84/fpu_polyg_wgs84_3_to_dieter

wdir=/project/ggcmi/AgMIP.output/processed/masks/weight
odir=/project/ggcmi/isi1.plots

if [ $plot = map ]; then
   /project/joshuaelliott/ggcmi/bin/plot.isi1/map.isi1.py -i $mfile -c $crop -a $afile -s $sfile -w $wdir/$crop.nc4 -p 0.1 -m $odir/plots/${plot}_${crop}.png -n $odir/files/${plot}_${crop}.nc4
elif [ $plot = box ]; then
   /project/joshuaelliott/ggcmi/bin/plot.isi1/box.isi1.py -i $mfile -c $crop -b $odir/plots/${plot}_${crop}.png -n $odir/files/${plot}_${crop}.nc4
else
   echo Unknown plotting option $plot . . .
   exit
fi
