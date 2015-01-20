#!/bin/bash

mod=$1

rootdir=/project/ggcmi/isi1
outdir=/project/ggcmi/isi1.clean
tempdir=/project/ggcmi/temp

epicfill=/project/ggcmi/isi1/epic.fill.2066_2068

gcms=(GFDL-ESM2M HadGEM2-ES IPSL-CM5A-LR MIROC-ESM-CHEM NorESM1-M)
crpf=(maize wheat soy rice)
crps=(mai whe soy ric)
rcps=(rcp2p6 rcp8p5)
irrs=(firr noirr)
co2s=(co2 noco2)

for g in ${gcms[@]}; do # GCM
    for ((c = 0; c < ${#crpf[@]}; c++)); do # crop
        cf=${crpf[$c]}
        cs=${crps[$c]}

        for co in ${co2s[@]}; do # co2
            if [ $co = noco2 ] && [ $g != HadGEM2-ES ]; then
                continue
            fi

            for ((i = 0; i < ${#irrs[@]}; i++)); do # irrigation
                irr=${irrs[$i]}

                # historical
                hdir=$rootdir/$mod/$g/hist/ssp2/$co/$irr/$cf
                if [ ! -d $hdir ]; then
                    echo Skipping $mod, $g, $cf, $irr, $co . . .
                fi
                if [ $mod = EPIC ]; then
                    hist=$hdir/${mod,,}_${g,,}_hist_ssp2_${co}_${irr}_yield_${cs}_annual_1980_2010.nc4
                    ncks -h -d time,0,24 $hist $tempdir/hist.nc4
                    ncap2 -O -h -s "time=time-79" $tempdir/hist.nc4 $tempdir/hist.nc4
                elif [ $mod = GEPIC ] || [ $mod = IMAGE_LEITAP ] || [ $mod = LPJ-GUESS ] || [ $mod = LPJmL ] || [ $mod = PEGASUS ]; then
                    if [ $g = HadGEM2-ES ]; then
                        ncrcat -h $hdir/*yield* $tempdir/hist.nc4
                        ncks -O -h -d time,9,33 $tempdir/hist.nc4 $tempdir/hist.nc4
                        ncap2 -O -h -s 'time[$time]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}' $tempdir/hist.nc4 $tempdir/hist.nc4
                    else
                        ncrcat -h $hdir/*yield* $tempdir/hist.nc4
                        ncks -O -h -d time,9,34 $tempdir/hist.nc4 $tempdir/hist.nc4
                        ncap2 -O -h -s 'time[$time]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25}' $tempdir/hist.nc4 $tempdir/hist.nc4
                    fi
                elif [ $mod = pDSSAT ]; then
                    if [ $g = HadGEM2-ES ]; then
                        ncrcat -h $hdir/*yield* $tempdir/hist.nc4
                        ncks -O -h -d time,29,53 $tempdir/hist.nc4 $tempdir/hist.nc4
                        ncap2 -O -h -s 'time[$time]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}' $tempdir/hist.nc4 $tempdir/hist.nc4
                    else
                        ncrcat -h $hdir/*yield* $tempdir/hist.nc4
                        ncks -O -h -d time,29,54 $tempdir/hist.nc4 $tempdir/hist.nc4
                        ncap2 -O -h -s 'time[$time]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25}' $tempdir/hist.nc4 $tempdir/hist.nc4
                    fi
                fi

                ncatted -O -h -a units,time,m,c,"years since 1980-01-01" $tempdir/hist.nc4 $tempdir/hist.nc4
                ncecat -O -h -u irr $tempdir/hist.nc4 $tempdir/hist.nc4 &> /dev/null
                ncap2 -O -h -s "irr[irr]=$i" $tempdir/hist.nc4 $tempdir/hist.nc4
                ncpdq -O -h -a time,irr $tempdir/hist.nc4 $tempdir/hist.nc4

                for r in ${rcps[@]}; do # RCP
                    echo Processing $mod, $g, $cf, $irr, $co, $r . . .

                    fdir=$rootdir/$mod/$g/$r/ssp2/$co/$irr/$cf
                    odir=$outdir/$mod/$g/$cf/$r/$co
                    fn=$odir/${mod,,}_${g,,}_ssp2_${co}_yield_${cs}_annual_1980_2099.nc4
                    mkdir -p $odir

                    # concatenate
                    if [ $mod = EPIC ]; then
                        # prepare files
                        f1=$fdir/${mod,,}_${g,,}_${r}_ssp2_${co}_${irr}_yield_${cs}_annual_2005_2035.nc4
                        f2=$fdir/${mod,,}_${g,,}_${r}_ssp2_${co}_${irr}_yield_${cs}_annual_2035_2065.nc4            
                        f3=$fdir/${mod,,}_${g,,}_${r}_ssp2_${co}_${irr}_yield_${cs}_annual_2069_2099.nc4
                        ncap2 -h -s "time=time-79" $f1 $tempdir/rcp1.nc4
                        ncatted -O -h -a units,time,m,c,"years since 1980-01-01" $tempdir/rcp1.nc4 $tempdir/rcp1.nc4
                        ncks -h -d time,1,30 $f2 $tempdir/rcp2.nc4
                        ncap2 -O -h -s "time=time-79" $tempdir/rcp2.nc4 $tempdir/rcp2.nc4
                        ncatted -O -h -a units,time,m,c,"years since 1980-01-01" $tempdir/rcp2.nc4 $tempdir/rcp2.nc4
                        ncap2 -h -s "time=time-79" $f3 $tempdir/rcp3.nc4
                        ncatted -O -h -a units,time,m,c,"years since 1980-01-01" $tempdir/rcp3.nc4 $tempdir/rcp3.nc4
                        ncrcat -h $tempdir/rcp1.nc4 $tempdir/rcp2.nc4 $epicfill.$cf.nc4 $tempdir/rcp3.nc4 $tempdir/rcp.final.nc4
                    elif [ $mod = GEPIC ] || [ $mod = IMAGE_LEITAP ] || [ $mod = LPJ-GUESS ] || [ $mod = LPJmL ] || [ $mod = pDSSAT ] || [ $mod = PEGASUS ]; then
                        if [ $g = HadGEM2-ES ]; then
                            for f in $fdir/*yield*; do
                                ncks -O -h --mk_rec_dim time $f $f
                            done
                            ncrcat -h $fdir/*yield* $tempdir/rcp.final.nc4
                            ncap2 -O -h -s 'time[$time]={25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119}' $tempdir/rcp.final.nc4 $tempdir/rcp.final.nc4
                        else
                            ncrcat -h $fdir/*yield* $tempdir/rcp.final.nc4
                            ncap2 -O -h -s 'time[$time]={26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119}' $tempdir/rcp.final.nc4 $tempdir/rcp.final.nc4
                        fi
                    fi

                    ncatted -O -h -a units,time,m,c,"years since 1980-01-01" $tempdir/rcp.final.nc4 $tempdir/rcp.final.nc4
                    ncecat -O -h -u irr $tempdir/rcp.final.nc4 $tempdir/rcp.final.nc4 &> /dev/null
                    ncap2 -O -h -s "irr[irr]=$i" $tempdir/rcp.final.nc4 $tempdir/rcp.final.nc4
                    ncpdq -O -h -a time,irr $tempdir/rcp.final.nc4 $tempdir/rcp.final.nc4
                    ncrcat -O -h $tempdir/hist.nc4 $tempdir/rcp.final.nc4 $tempdir/rcp.final.nc4
                    ncpdq -O -h -a irr,time $tempdir/rcp.final.nc4 $tempdir/rcp.final.nc4
                    ncatted -O -h -a long_name,irr,c,c,"ir, rf" $tempdir/rcp.final.nc4 $tempdir/rcp.final.nc4

                    if [ $i = 0 ]; then
                        mv $tempdir/rcp.final.nc4 $fn
                    else
                        ncrcat -O -h $fn $tempdir/rcp.final.nc4 $fn
                        nccopy -d9 -k4 $fn $fn.2
                        mv $fn.2 $fn
                    fi
                    rm -f $tempdir/rcp*
                done # RCP

                rm -f $tempdir/hist.nc4
            done # irrigation
        done # co2
    done # crop
done # GCM