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

indir=/project/ggcmi/isi1/raw/LPJmL-extra
outdir=/project/ggcmi/isi1/raw/LPJmL

run=hist

for m in GFDL-ESM2M HadGEM2-ES IPSL-CM5A-LR MIROC-ESM-CHEM NorESM1-M; do
    for c in co2 noco2; do
        for i in firr noirr; do
            for crop in maize wheat soy rice; do
                crops=$(shortnames $crop)

                odir=$outdir/$m/$run/ssp2/$c/$i/$crop
                if [ ! -d $odir ]; then
                    echo Skipping . . .
                    continue
                fi

                f1=$(ls $indir/*${m,,}*rcp2p6*_${c}*${i}*${crops}*1951*)
                f2=$(ls $indir/*${m,,}*rcp2p6*_${c}*${i}*${crops}*1961*)

                f1new=${f1/_rcp2p6/}
                f1new=${f1new%?}
                f2new=${f2/_rcp2p6/}
                f2new=${f2new%?}

                echo "$odir corresponds to $f1 ($f1new), $f2 ($f2new)"
                cp $f1 $odir/$(basename $f1new)
                cp $f2 $odir/$(basename $f2new)
            done
        done
    done
done