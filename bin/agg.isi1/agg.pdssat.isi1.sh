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

g=$1
cf=$2
co=$3
r=$4

mod=pDSSAT

idir=/project/ggcmi/isi1/isi1.pdssat.clean/$mod/$g/$cf/$r/$co
odir=/project/ggcmi/isi1/isi1.pdssat.agg

wdir=/project/ggcmi/AgMIP.output/processed/masks/weight
mfile=/project/ggcmi/AgMIP.output/processed/masks/aggr/fpu.gadm.global.mask.nc4

cs=$(shortnames $cf)
f=$idir/${mod,,}_${g,,}_ssp2_${co}_yield_${cs}_annual_1951_2099.nc4

if [ -f $f ]; then
   fdir=$odir/$mod/$g/$cf/$r/$co
   mkdir -p $fdir
   /project/joshuaelliott/ggcmi/bin/agg.isi1/agg.out.py -i $f:yield_${cs} -w $wdir/$cf.nc4 -a $mfile -l time -n 10 -y yield_${cs} -o $fdir/$(basename $f) 
else
   echo No data for $mod, $g, $cf, $co, $r . . .
   exit
fi
