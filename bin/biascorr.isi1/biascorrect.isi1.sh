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

idir=/project/ggcmi/isi1/processed/isi1.agg/$mod/$g/$cf/$r/$co
odir=/project/ggcmi/isi1/processed/isi1.biascorr/$mod/$g/$cf/$r/$co

cs=$(shortnames $cf)
if=$idir/${mod,,}_${g,,}_ssp2_${co}_yield_${cs}_annual_1980_2099.nc4
of=$odir/$(basename $if)

fpufile=/project/ggcmi/AgMIP.input/other.inputs/reference/iizumi/iizumi.1982-2006.fpu.fixed_iizumi_mask.nc4
globalfile=/project/ggcmi/AgMIP.input/other.inputs/reference/iizumi/iizumi.1982-2006.global.fixed_iizumi_mask.nc4

if [ -f $if ]; then
   mkdir -p $odir 

   # fpu
   /project/joshuaelliott/ggcmi/bin/biascorr.isi1/biascorrect.isi1.py -i $if -r $fpufile -a fpu -o $odir
   mv $of $of.2
   ncrename -O -h -v yield_detrend,yield_${cs}_fpu $of.2 $of.2
   ncks -O -h -x -v yield_retrend $of.2 $of.2

   # global
   /project/joshuaelliott/ggcmi/bin/biascorr.isi1/biascorrect.isi1.py -i $if -r $globalfile -a global -o $odir
   ncrename -O -h -v yield_detrend,yield_${cs}_global $of $of
   ncks -O -h -x -v yield_retrend $of $of

   # append
   ncks -h -A $of.2 $of
   rm $of.2
else
   echo No data for $mod, $g, $cf, $co, $r . . .
   exit
fi
