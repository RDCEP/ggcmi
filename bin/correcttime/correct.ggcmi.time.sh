#!/bin/bash

cropname() {
  n=$1
  if [ $n = "maize" ]; then
    echo "Maize"
  elif [ $n = "sorghum" ]; then
    echo "Sorghum"
  elif [ $n = "soy" ]; then
    echo "Soybeans"
  elif [ $n = "wheat" ]; then
    echo "Wheat"
  elif [ $n = "rice" ]; then
    echo "Rice"
  elif [ $n = "sugar_beet" ]; then
    echo "Sugar_beet"
  elif [ $n = "rye" ]; then
    echo "Rye"
  elif [ $n = "rapeseed" ]; then
    echo "Rapeseed"
  elif [ $n = "millet" ]; then
    echo "Millet"
  elif [ $n = "barley" ]; then
    echo "Barley"
  elif [ $n = "sugarcane" ]; then
    echo "Sugar_cane"
  else
    echo ""
  fi
}

irtag() {
  t=$1
  if [ $t = "firr" ]; then
    echo "ir"
  elif [ $t = "noirr" ]; then
    echo "rf"
  else
    echo ""
  fi
}

dir=/project/joshuaelliott/ggcmi/AgMIP.output/EPIC-IIASA
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
              crop=$(cropname $j)
              ir=$(irrtag ${scensplit[1]})
              gfile=$(ls $inputdir/$crop"_"$ir*)
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
