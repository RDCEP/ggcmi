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

idir=/project/ggcmi/isi1/isi1.long.agg/$mod/$g/$cf/$r/$co
odir=/project/ggcmi/isi1/isi1.long.biascorr/$mod/$g/$cf/$r/$co

fpufile=/project/ggcmi/AgMIP.input/other.inputs/reference/iizumi/iizumi.1982-2006.fpu.fixed_iizumi_mask.nc4
gadmfile=/project/ggcmi/AgMIP.input/other.inputs/reference/iizumi/iizumi.1982-2006.gadm0.fixed_iizumi_mask.nc4
globalfile=/project/ggcmi/AgMIP.input/other.inputs/reference/iizumi/iizumi.1982-2006.global.fixed_iizumi_mask.nc4

if [ -d $idir ] && [ $(ls $idir | wc -l) = 1 ]; then
   mkdir -p $odir

   cs=$(shortnames $cf)
   if=$idir/$(ls $idir)
   of=$odir/$(basename $if)

   # fpu
   /project/joshuaelliott/ggcmi/bin/biascorr.isi1/biascorrect.isi1.py -i $if -r $fpufile -a fpu -o $odir
   mv $of $of.2
   ncrename -O -h -v yield_detrend,yield_${cs}_fpu $of.2 $of.2
   ncks -O -h -x -v yield_retrend $of.2 $of.2

   # gadm
   /project/joshuaelliott/ggcmi/bin/biascorr.isi1/biascorrect.isi1.py -i $if -r $gadmfile -a gadm0 -o $odir
   mv $of $of.3
   ncrename -O -h -v yield_detrend,yield_${cs}_gadm0 $of.3 $of.3
   ncks -O -h -x -v yield_retrend $of.3 $of.3

   # global
   /project/joshuaelliott/ggcmi/bin/biascorr.isi1/biascorrect.isi1.py -i $if -r $globalfile -a global -o $odir
   ncrename -O -h -v yield_detrend,yield_${cs}_global $of $of
   ncks -O -h -x -v yield_retrend $of $of

   # append
   ncks -h -A $of.2 $of
   ncks -h -A $of.3 $of
   rm $of.2 $of.3
else
   echo No data for $mod, $g, $cf, $co, $r . . .
   exit
fi
