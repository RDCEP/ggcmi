#!/bin/bash

root=/project/ggcmi

clong=(maize soy wheat rice sorghum millet)
cshort=(mai soy whe ric sor mil)

# Header
echo irfile rffile bcfile wtfile vname outfile

for m in pDSSAT pAPSIM GEPIC PEGASUS LPJmL EPIC-IIASA EPIC-Boku LPJ-GUESS; do
  for w in AgMERRA AgCFSR WFDEI.GPCC WFDEI.CRU Princeton GRASP WATCH ERAI CFSR; do
    for ((i = 0; i < ${#clong[@]}; i++)); do
      for v in yield; do
        irfiles=($(ls $root/AgMIP.output/$m/$w/${clong[$i]}/*firr*$v*${cshort[$i]}*                          2> /dev/null))
        rffiles=($(ls $root/AgMIP.output/$m/$w/${clong[$i]}/*noirr*$v*${cshort[$i]}*                         2> /dev/null))
        bcfile=$(ls $root/AgMIP.output/processed/biascorr/gadm0/faostat/fixed_mask/${m}*_${w}*${cshort[$i]}* 2> /dev/null)
        wtfile=$root/AgMIP.output/processed/masks/weight/${clong[$i]}.nc4
        if [ $bcfile ]; then # bias-corrected file exists
          for ((j = 0; j < ${#irfiles[@]}; j++)); do
            mkdir -p $root/AgMIP.output/processed/rescaled/gadm0/faostat/fixed_mask
            outfile=$root/AgMIP.output/processed/rescaled/gadm0/faostat/fixed_mask/$(basename ${irfiles[$j]//_firr})
            outfile=${outfile/nc4/rescaled.nc4}
            echo ${irfiles[$j]} ${rffiles[$j]} $bcfile $wtfile ${v}_${cshort[$i]} $outfile
          done
        fi
      done
    done
  done
done
