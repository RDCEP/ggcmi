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
var=$6

rootdir=/project/ggcmi/isi1/raw
outdir=/project/ggcmi/isi1/isi1.long.clean

epicfill=/project/ggcmi/isi1/bin/epic.fill.2066_2068

irrs=(firr noirr)

cs=$(shortnames $cf)

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
elif [ $mod = GEPIC ] || [ $mod = LPJ-GUESS ] || [ $mod = LPJmL ] || [ $mod = PEGASUS ] || [ $mod = IMAGE_LEITAP ]; then
    if [ $mod = LPJmL ] && [ $var = yield ]; then # yield extends to 1951 for LPJmL
        syear=1951
        if [ $g = HadGEM2-ES ]; then
            timeh=$(seq -s, 0 53)
            timef=$(seq -s, 54 148)
        else
            timeh=$(seq -s, 0 54)
            timef=$(seq -s, 55 148)
        fi
    else
        syear=1971
        if [ $g = HadGEM2-ES ]; then
            timeh=$(seq -s, 0 33)
            timef=$(seq -s, 34 128)
        else
            timeh=$(seq -s, 0 34)
            timef=$(seq -s, 35 128)
        fi
    fi
elif [ $mod = EPIC ]; then
    syear=1980
    timeh=$(seq -s, 0 24)
    timef=$(seq -s, 25 119)
else
    echo Unrecognized model $mod
    exit
fi

for ((i = 0; i < ${#irrs[@]}; i++)); do # irrigation
    irr=${irrs[$i]}

    hdir=$rootdir/$mod/$g/hist/ssp2/$co/$irr/$cf
    fdir=$rootdir/$mod/$g/$r/ssp2/$co/$irr/$cf

    ls $hdir/*${var}* >/dev/null 2>&1
    hasvar1=$?
    ls $fdir/*${var}* >/dev/null 2>&1
    hasvar2=$?

    if  [ $hasvar1 != 0 ] || [ $hasvar2 != 0 ]; then
        echo Skipping $mod, $g, $cf, $irr, $co, $r, $var . . .
        exit
    fi
    echo Processing $mod, $g, $cf, $irr, $co, $r, $var . . .

    # make time record dimension in some files
    if [ $g = HadGEM2-ES ]; then
        for f in $fdir/*${var}*; do
            ncks -O -h --mk_rec_dim time $f $f
        done
    fi

    odir=$outdir/$mod/$g/$cf/$r/$co
    fn=$odir/${mod,,}_${g,,}_ssp2_${co}_${var}_${cs}_annual_${syear}_2099.nc4
    mkdir -p $odir

    # historical
    if [ $mod = LPJmL ]; then
        ncrcat -h $hdir/*${var}* hist.nc4
    else
        ncrcat -h $hdir/*${var}*.nc4 hist.nc4 # locate nc4 specifically
    fi
    if [ $mod = EPIC ]; then
        ncks -O -h -d time,0,24 hist.nc4 hist.nc4 # select first 25 years
    fi
    ncap2 -O -h -s "time(:)={$timeh}" hist.nc4 hist.nc4
    ncatted -O -h -a units,time,m,c,"years since ${syear}-01-01" hist.nc4 hist.nc4
    ncecat -O -h -u irr hist.nc4 hist.nc4 &> /dev/null
    ncap2 -O -h -s "irr[irr]=$i" hist.nc4 hist.nc4
    ncpdq -O -h -a time,irr hist.nc4 hist.nc4

    # future
    if [ $mod = EPIC ]; then
        f1=$fdir/${mod,,}_${g,,}_${r}_ssp2_${co}_${irr}_${var}_${cs}_annual_2005_2035.nc4
        f2=$fdir/${mod,,}_${g,,}_${r}_ssp2_${co}_${irr}_${var}_${cs}_annual_2035_2065.nc4            
        f3=$fdir/${mod,,}_${g,,}_${r}_ssp2_${co}_${irr}_${var}_${cs}_annual_2069_2099.nc4
        ncks -O -h -d time,1,30 $f2 $f2.2
        cp $epicfill.$cf.nc4 epicfill.tmp.nc4
        if [ $var != yield ]; then
            ncrename -O -h -v yield_${cs},${var}_${cs} epicfill.tmp.nc4 epicfill.tmp.nc4
        fi
        ncrcat -h $f1 $f2.2 epicfill.tmp.nc4 $f3 rcp.final.nc4
        rm $f2.2 epicfill.tmp.nc4
    elif [ $mod = LPJmL ]; then
        ncrcat -h $fdir/*${var}* rcp.final.nc4
    else
        ncrcat -h $fdir/*${var}*.nc4 rcp.final.nc4 # locate nc4 specifically
    fi
    ncap2 -O -h -s "time(:)={$timef}" rcp.final.nc4 rcp.final.nc4
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
