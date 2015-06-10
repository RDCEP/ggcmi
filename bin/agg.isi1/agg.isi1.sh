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

idir=/project/ggcmi/isi1/processed/isi1.clean/$mod/$g/$cf/$r/$co
odir=/project/ggcmi/isi1/processed/isi1.agg

wdir=/project/ggcmi/AgMIP.output/processed/masks/weight
mfile=/project/ggcmi/AgMIP.output/processed/masks/aggr/fpu.global.mask.nc4

cs=$(shortnames $cf)
f=$idir/${mod,,}_${g,,}_ssp2_${co}_yield_${cs}_annual_1980_2099.nc4

if [ -f $f ]; then
   fdir=$odir/$mod/$g/$cf/$r/$co
   mkdir -p $fdir
   /project/joshuaelliott/ggcmi/bin/agg.isi1/agg.out.py -i $f:yield_${cs} -w $wdir/$cf.nc4 -a $mfile -l time -n 10 -y yield_${cs} -o $fdir/$(basename $f) 
else
   echo No data for $mod, $g, $cf, $co, $r . . .
   exit
fi
