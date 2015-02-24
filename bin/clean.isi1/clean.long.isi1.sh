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
r=$5

if [ $mod != pDSSAT ] && [ $mod != GEPIC ]; then
    echo Unrecognized model $mod. Exiting . . .
fi

rootdir=/project/ggcmi/isi1/raw
outdir=/project/ggcmi/isi1/isi1.long.clean

irrs=(firr noirr)

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

    # make time record dimension in some files
    if [ $g = HadGEM2-ES ]; then
        for f in $fdir/*yield*; do
            ncks -O -h --mk_rec_dim time $f $f
        done
    fi

    # time ranges
    if [ $mod = pDSSAT ]; then
       syear=1951
       if [ $g = HadGEM2-ES ]; then
            timeh=$(seq -s, 0 53)
            timef=$(seq -s, 54 148)
        else
            timeh=$(seq -s, 0 54)
            timef=$(seq -s, 55 148)
        fi
    elif [ $mod = GEPIC ]; then
        syear=1971
        if [ $g = HadGEM2-ES ]; then
            timeh=$(seq -s, 0 33)
            timef=$(seq -s, 34 128)
        else
            timeh=$(seq -s, 0 34)
            timef=$(seq -s, 35 128)
        fi
    fi

    # historical
    ncrcat -h $hdir/*yield* hist.nc4
    ncap2 -O -h -s 'time(:)={${timeh%?}}' hist.nc4 hist.nc4
    ncatted -O -h -a units,time,m,c,"years since ${syear}-01-01" hist.nc4 hist.nc4
    ncecat -O -h -u irr hist.nc4 hist.nc4 &> /dev/null
    ncap2 -O -h -s "irr[irr]=$i" hist.nc4 hist.nc4
    ncpdq -O -h -a time,irr hist.nc4 hist.nc4

    odir=$outdir/$mod/$g/$cf/$r/$co
    fn=$odir/${mod,,}_${g,,}_ssp2_${co}_yield_${cs}_annual_${syear}_2099.nc4
    mkdir -p $odir

    # future
    ncrcat -h $fdir/*yield* rcp.final.nc4
    ncap2 -O -h -s 'time(:)={${timef%?}}' rcp.final.nc4 rcp.final.nc4
    ncatted -O -h -a units,time,m,c,"years since ${syear}-01-01" rcp.final.nc4 rcp.final.nc4
    ncecat -O -h -u irr rcp.final.nc4 rcp.final.nc4 &> /dev/null
    ncap2 -O -h -s "irr[irr]=$i" rcp.final.nc4 rcp.final.nc4
    ncpdq -O -h -a time,irr rcp.final.nc4 rcp.final.nc4

    # concatenate
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
