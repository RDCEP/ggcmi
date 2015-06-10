#!/bin/bash

PATH=$PATH:/project/joshuaelliott/ggcmi/utils

function shortnames {
    fname=$1
    if [ $fname = maize ]; then
        sname="mai"
    elif [ $fname = wheat ]; then
        sname="whe"
    elif [ $fname = soy ]; then
        sname="soy"
    elif [ $fname = rice ]; then
        sname="ric"
    else
        sname=""
    fi
    echo $sname
}

mod=$1
g=$2
cf=$3
co=$4
r=$5
var=$6
year=$7

idir=/project/ggcmi/isi1/processed/isi1.long.clean/$mod/$g/$cf/$r/$co
odir=/project/ggcmi/isi1/processed/isi1.long.agg.ray.slices/$mod/$g/$cf/$r/$co/$year

wdir=/project/joshuaelliott/ggcmi/reference/ray/masks
mfile=/project/ggcmi/AgMIP.output/processed/masks/aggr/gadm0.mask.nc4

ls $idir/*${var}* >/dev/null 2>&1

if [ $? = 0 ] && [ $(ls $idir/*${var}* | wc -l) = 1 ]; then
   f=$(ls $idir/*${var}*)
   cs=$(shortnames $cf)
   mkdir -p $odir
   /project/joshuaelliott/ggcmi/bin/agg.isi1/agg.out.py -i $f:${var}_${cs} -w $wdir/${cf}.ray.${year}.nc4 -a $mfile -l time -n 10 -y ${var}_${cs} -o $odir/$(basename $f) 
else
   echo No data for $mod, $g, $cf, $co, $r, $var . . .
   exit
fi
