#!/bin/bash

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

rcp26dir=/project/ggcmi/isi1/processed/isi1.biascorr/$mod/$g/$cf/rcp2p6/$co
rcp85dir=/project/ggcmi/isi1/processed/isi1.biascorr/$mod/$g/$cf/rcp8p5/$co
odir=/project/ggcmi/isi1/processed/isi1.metrics/$mod/$g/$cf/$co

cs=$(shortnames $cf)

rcp26file=$rcp26dir/${mod,,}_${g,,}_ssp2_${co}_yield_${cs}_annual_1980_2099.nc4
rcp85file=$rcp85dir/$(basename $rcp26file)

if [ -f $rcp26file ] && [ -f $rcp85file ]; then
   mkdir -p $odir
   /project/joshuaelliott/ggcmi/bin/metrics.isi1/metrics.isi1.py --rcp26file $rcp26file --rcp85file $rcp85file -v yield_${cs} -o $odir/$(basename $rcp26file) 
else
   echo No data for $mod, $g, $cf, $co . . .
   exit
fi
