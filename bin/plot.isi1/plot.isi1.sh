#!/bin/bash

plot=$1
crop=$2

mfile=/project/ggcmi/isi1/isi1.metrics.all/metrics.nc4
afile=/project/ggcmi/AgMIP.output/processed/masks/aggr/fpu.mask.nc4
hffile=/project/ggcmi/AgMIP.output/processed/masks/weight/aggs/all.fpu.nc4
hgfile=/project/ggcmi/AgMIP.output/processed/masks/weight/aggs/all.global.nc4
sfile=/project/joshuaelliott/data/masks/fpu_shp_files/fpu_wgs84/fpu_polyg_wgs84_3_to_dieter

wdir=/project/ggcmi/AgMIP.output/processed/masks/weight
odir=/project/ggcmi/isi1/isi1.plots

if [ $plot = blmap ]; then
   /project/joshuaelliott/ggcmi/bin/plot.isi1/blmap.isi1.py -i $mfile                          \
                                                            -c $crop                           \
                                                            -a $afile                          \
                                                            -r $hffile                         \
                                                            -s $sfile                          \
                                                            -w $wdir/$crop.nc4                 \
                                                            -p 0.1                             \
                                                            -m $odir/plots/${plot}_${crop}.png \
                                                            -n $odir/files/${plot}_${crop}.nc4
elif [ $plot = dymap ]; then
   for var in delta_yield_26 delta_yield_85; do
      /project/joshuaelliott/ggcmi/bin/plot.isi1/dymap.isi1.py -i $mfile                                 \
                                                               -c $crop                                  \
                                                               -a $afile                                 \
                                                               -r $hffile                                \
                                                               -s $sfile                                 \
                                                               -w $wdir/$crop.nc4                        \
                                                               -p 0.1                                    \
                                                               -v $var                                   \
                                                               -m $odir/plots/${plot}_${var}_${crop}.png \
                                                               -n $odir/files/${plot}_${var}_${crop}.nc4
   done
elif [ $plot = cmap ]; then
   for var in delta_yield_26 delta_yield_85; do
      /project/joshuaelliott/ggcmi/bin/plot.isi1/cmap.isi1.py -i $mfile                                 \
                                                              -c $crop                                  \
                                                              -a $afile                                 \
                                                              -r $hffile                                \
                                                              -s $sfile                                 \
                                                              -w $wdir/$crop.nc4                        \
                                                              -p 0.1                                    \
                                                              -v $var                                   \
                                                              -m $odir/plots/${plot}_${var}_${crop}.png \
                                                              -n $odir/files/${plot}_${var}_${crop}.nc4
   done
elif [ $plot = box ]; then
   for var in delta_yield_26 delta_yield_85; do
      /project/joshuaelliott/ggcmi/bin/plot.isi1/box.isi1.py -i $mfile                                 \
                                                             -c $crop                                  \
                                                             -r $hgfile                                \
                                                             -v $var                                   \
                                                             -b $odir/plots/${plot}_${var}_${crop}.png \
                                                             -n $odir/files/${plot}_${var}_${crop}.nc4
   done
else
   echo Unknown plotting option $plot . . .
   exit
fi
