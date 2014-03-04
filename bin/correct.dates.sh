#!/bin/bash

dir=/project/joshuaelliott/ggcmi/AgMIP.output/LPJ-GUESS
inputdir=/project/joshuaelliott/data/ggcmi/other.inputs/AGMIP_GROWING_SEASON.HARM.version1.23
scens=(default_firr default_noirr harmnon_firr harmnon_noirr fullharm_firr fullharm_noirr)
corrector=/project/joshuaelliott/ggcmi/bin/ggcmi.corrector.py
correctorInput=/project/joshuaelliott/ggcmi/bin/ggcmi.corrector.input.py

for i in `ls $dir`; do # weather
  for j in `ls $dir"/"$i`; do # crop
    dirij=$dir"/"$i"/"$j
    echo "Scanning directory "$dirij" . . ."
    if [[ $j != "others" ]] && [[ `ls $dirij` ]]; then # if not empty
      for s in ${scens[@]}; do
        if [[ `ls $dirij |  grep $s` ]]; then # if has scenario
          files=(`ls $dirij/*$s*`) # all files
          if [[ `echo ${files[@]} | grep "plant-day"` ]] && [[ `echo ${files[@]} | grep "maty-day"` ]]; then
            pfile=$(ls $dirij/*$s"_plant-day"*) # planting file
            mfile=$(ls $dirij/*$s"_maty-day"*) # maturity file
            for f in ${files[@]}; do # run corrector
              echo "Running corrector for "$f" . . ."
              `$corrector -i $f -p $pfile -m $mfile -o $f.new`
            done
            for f in ${files[@]}; do # rename files
              mv $f $f.old
              mv $f.new $f
            done
          else
            OLDIFS=$IFS
            IFS='_'
            scensplit=($s)
            IFS=$OLDIFS
            if [[ ${scensplit[0]} == "fullharm" ]] || [[ ${scensplit[0]} == "harmnon" ]]; then # correct harmonized data
              gfile=$(ls $inputdir/$j"_"${scensplit[1]}*)
              if [[ $gfile ]]; then
                for f in ${files[@]}; do # run corrector
                  echo "Running corrector with harmonized input for "$f" . . ."
                  `$correctorInput -i $f -g $gfile -o $f.new`
                done
                for f in ${files[@]}; do # rename files
                  mv $f $f.old
                  mv $f.new $f
                done
              fi
            fi
          fi
        fi
      done
    fi
  done
done
