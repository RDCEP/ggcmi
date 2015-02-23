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

g=$1
cf=$2
co=$3
r=$4

rootdir=/project/ggcmi/isi1/raw
outdir=/project/ggcmi/isi1/isi1.pdssat.clean

mod=pDSSAT # only pDSSAT

irrs=(firr noirr)

if [ $co = noco2 ] && [ $g != HadGEM2-ES ]; then
    echo "noco2 data only available for HadGEM2-ES"
    exit
fi

cs=$(shortnames $cf)

for ((i = 0; i < ${#irrs[@]}; i++)); do # irrigation
    irr=${irrs[$i]}

    hdir=$rootdir/$mod/$g/hist/ssp2/$co/$irr/$cf
    fdir=$rootdir/$mod/$g/$r/ssp2/$co/$irr/$cf
    if [ ! -d $hdir ] || [ ! -d $fdir ]; then
        echo Skipping $mod, $g, $cf, $irr, $co, $r . . .
        exit
    fi
    echo Processing $mod, $g, $cf, $irr, $co, $r . . .

    # historical
    if [ $g = HadGEM2-ES ]; then
        ncrcat -h $hdir/*yield* hist.nc4
        ncap2 -O -h -s 'time(:)={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53}' hist.nc4 hist.nc4
    else
        ncrcat -h $hdir/*yield* hist.nc4
        ncap2 -O -h -s 'time(:)={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54}' hist.nc4 hist.nc4
    fi

    ncatted -O -h -a units,time,m,c,"years since 1951-01-01" hist.nc4 hist.nc4
    ncecat -O -h -u irr hist.nc4 hist.nc4 &> /dev/null
    ncap2 -O -h -s "irr[irr]=$i" hist.nc4 hist.nc4
    ncpdq -O -h -a time,irr hist.nc4 hist.nc4

    odir=$outdir/$mod/$g/$cf/$r/$co
    fn=$odir/${mod,,}_${g,,}_ssp2_${co}_yield_${cs}_annual_1951_2099.nc4
    mkdir -p $odir

    # concatenate
    if [ $g = HadGEM2-ES ]; then
        for f in $fdir/*yield*; do
            ncks -O -h --mk_rec_dim time $f $f
        done
        ncrcat -h $fdir/*yield* rcp.final.nc4
        ncap2 -O -h -s 'time(:)={54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148}' rcp.final.nc4 rcp.final.nc4
    else
        ncrcat -h $fdir/*yield* rcp.final.nc4
        ncap2 -O -h -s 'time(:)={55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148}' rcp.final.nc4 rcp.final.nc4
    fi

    ncatted -O -h -a units,time,m,c,"years since 1951-01-01" rcp.final.nc4 rcp.final.nc4
    ncecat -O -h -u irr rcp.final.nc4 rcp.final.nc4 &> /dev/null
    ncap2 -O -h -s "irr[irr]=$i" rcp.final.nc4 rcp.final.nc4
    ncpdq -O -h -a time,irr rcp.final.nc4 rcp.final.nc4
    ncrcat -O -h hist.nc4 rcp.final.nc4 rcp.final.nc4
    ncpdq -O -h -a irr,time rcp.final.nc4 rcp.final.nc4
    ncatted -O -h -a long_name,irr,c,c,"ir, rf" rcp.final.nc4 rcp.final.nc4

    if [ $i = 0 ]; then
        mv rcp.final.nc4 $fn
    else
        ncrcat -O -h $fn rcp.final.nc4 $fn
        nccopy -d9 -k4 $fn $fn.2
        mv $fn.2 $fn
    fi
    rm -f rcp* hist.nc4
done
