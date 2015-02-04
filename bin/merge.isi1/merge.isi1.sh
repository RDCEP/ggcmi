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

models=(EPIC GEPIC IMAGE_LEITAP LPJ-GUESS LPJmL pDSSAT PEGASUS)
gcms=(GFDL-ESM2M HadGEM2-ES IPSL-CM5A-LR MIROC-ESM-CHEM NorESM1-M)
crops=(maize wheat soy rice)
co2s=(co2 noco2)

indir=/project/ggcmi/isi1/isi1.metrics

outfile=/project/ggcmi/isi1/isi1.metrics.all/metrics.nc4

firstfile=true
for m in ${models[@]}; do
    for g in ${gcms[@]}; do
        for c in ${crops[@]}; do
            cs=$(shortnames $c)
            for co in ${co2s[@]}; do
                ml=${m,,}
                gl=${g,,}

                infile=$indir/$m/$g/$c/$co/${ml}_${gl}_ssp2_${co}_yield_${cs}_annual_1980_2099.nc4
                if [ ! -f $infile ]; then
                    continue
                fi

                if [ $firstfile = true ]; then
                    cp $infile $outfile
                    firstfile=false
                else
                    ncks -h -A $infile $outfile
                fi

                ncrename -O -h -v beta_fpu,beta_fpu_${ml}_${gl}_${c}_${co} $outfile $outfile
                ncrename -O -h -v beta_global,beta_global_${ml}_${gl}_${c}_${co} $outfile $outfile
                ncrename -O -h -v lambda_fpu,lambda_fpu_${ml}_${gl}_${c}_${co} $outfile $outfile
                ncrename -O -h -v lambda_global,lambda_global_${ml}_${gl}_${c}_${co} $outfile $outfile
            done
        done
    done
done

nccopy $outfile $outfile.2
mv $outfile.2 $outfile