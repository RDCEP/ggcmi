#!/bin/bash

indir=/project/ggcmi/phase1.final
outdir=/project/ggcmi/phase1.final/processed.simple
models=(pDSSAT pAPSIM LPJ-GUESS LPJmL PEGASUS GEPIC EPIC-IIASA EPIC-Boku EPIC-test CGMS-WOFOST CLM-Crop EPIC-TAMU ORCHIDEE ORCHIDEE-crop PEPIC PRYSBI2)
weathers="AgMERRA AgCFSR WFDEI.GPCC WFDEI.CRU Princeton GRASP WATCH ERAI CFSR"
cropslong=(maize soy wheat rice sorghum millet)
cropsshort=(mai soy whe ric sor mil)

# Header
echo irfile rffile bcfile wtfile vname outfile

for model in ${models[@]}; do
  for weather in $weathers; do
    for ((i = 0; i < ${#cropslong[@]}; i++)); do
      for v in yield; do
        irfiles=($(ls $indir/$model/$weather/${cropslong[$i]}/*firr*$v*${cropsshort[$i]}* 2>/dev/null))
        rffiles=($(ls $indir/$model/$weather/${cropslong[$i]}/*noirr*$v*${cropsshort[$i]}* 2>/dev/null))
        bcfile=$(ls $outdir/biascorr/gadm0/faostat/fixed_mirca_mask/${model,,}*_${weather,,}*${cropsshort[$i]}* 2>/dev/null)
        wtfile=$outdir/masks/weight/${cropslong[$i]}.nc4
        if [ ${#irfiles[@]} -ge 1 ] && [ $bcfile ]; then # irrigated and bias-corrected file exist
          for ((j = 0; j < ${#irfiles[@]}; j++)); do
            if [ ${#rffiles[@]} = 0 ]; then
              rffile=None # no rainfed file
            else
              rffile=${rffiles[$j]}
            fi
            mkdir -p $outdir/rescaled/gadm0/faostat/fixed_mirca_mask
            outfile=$outdir/rescaled/gadm0/faostat/fixed_mirca_mask/$(basename ${irfiles[$j]//_firr})
            outfile=${outfile/nc4/rescaled.nc4}
            echo ${irfiles[$j]} $rffile $bcfile $wtfile ${v}_${cropsshort[$i]} $outfile
          done
        fi
      done
    done
  done
done
