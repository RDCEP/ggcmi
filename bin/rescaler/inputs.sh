#!/bin/bash

root=/project/ggcmi

mods=(pDSSAT pAPSIM LPJ-GUESS LPJmL PEGASUS GEPIC EPIC-IIASA EPIC-Boku CGMS-WOFOST CLM-Crop EPIC-TAMU ORCHIDEE ORCHIDEE-crop PEPIC PRYSBI2)

clong=(maize soy wheat rice sorghum millet)
cshort=(mai soy whe ric sor mil)

# Header
echo irfile rffile bcfile wtfile vname outfile

for m in ${mods[@]}; do
  for w in AgMERRA AgCFSR WFDEI.GPCC WFDEI.CRU Princeton GRASP WATCH ERAI CFSR; do
    for ((i = 0; i < ${#clong[@]}; i++)); do
      for v in yield; do
        irfiles=($(ls $root/AgMIP.output/$m/$w/${clong[$i]}/*firr*$v*${cshort[$i]}*                                    2> /dev/null))
        rffiles=($(ls $root/AgMIP.output/$m/$w/${clong[$i]}/*noirr*$v*${cshort[$i]}*                                   2> /dev/null))
        bcfile=$(ls $root/AgMIP.output/processed/biascorr/gadm0/faostat/fixed_mirca_mask/${m,,}*_${w,,}*${cshort[$i]}* 2> /dev/null)
        wtfile=$root/AgMIP.output/processed/masks/weight/${clong[$i]}.nc4
        if [ ${#irfiles[@]} -ge 1 ] && [ $bcfile ]; then # irrigated and bias-corrected file exist
          for ((j = 0; j < ${#irfiles[@]}; j++)); do
            if [ ${#rffiles[@]} = 0 ]; then
              rffile=None # no rainfed file
            else
              rffile=${rffiles[$j]}
            fi
            mkdir -p $root/AgMIP.output/processed/rescaled/gadm0/faostat/fixed_mirca_mask
            outfile=$root/AgMIP.output/processed/rescaled/gadm0/faostat/fixed_mirca_mask/$(basename ${irfiles[$j]//_firr})
            outfile=${outfile/nc4/rescaled.nc4}
            echo ${irfiles[$j]} $rffile $bcfile $wtfile ${v}_${cshort[$i]} $outfile
          done
        fi
      done
    done
  done
done
