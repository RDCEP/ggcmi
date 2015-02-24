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

idir=/project/ggcmi/isi1/isi1.long.clean/$mod/$g/$cf/$r/$co
odir=/project/ggcmi/isi1/isi1.long.agg/$mod/$g/$cf/$r/$co

wdir=/project/ggcmi/AgMIP.output/processed/masks/weight
mfile=/project/ggcmi/AgMIP.output/processed/masks/aggr/fpu.gadm.global.mask.nc4

if [ -d $idir ] && [ $(ls $idir | wc -l) = 1 ]; then
   f=$idir/$(ls $idir)
   cs=$(shortnames $cf)
   mkdir -p $odir
   /project/joshuaelliott/ggcmi/bin/agg.isi1/agg.out.py -i $f:yield_${cs} -w $wdir/$cf.nc4 -a $mfile -l time -n 10 -y yield_${cs} -o $odir/$(basename $f) 
else
   echo No data for $mod, $g, $cf, $co, $r . . .
   exit
fi
