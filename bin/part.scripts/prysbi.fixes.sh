#!/bin/bash

crops="rice maize soy wheat"
vars="yield aet biom maty-day"
clim=(AgCFSR AgMERRA GRASP Princeton WATCH WFDEI.GPCC)
YYYY=(1980 1980 1961 1948 1958 1979)

for c in $crops ; do 
  for v in $vars ; do 
    for i in {0..5} ; do
      folder=./${clim[$i]}/$c
      file=`ls $folder | grep $v"_"${c:0:3}`
      echo $folder/$file
      ncrename -v var,$v"_"${c:0:3} $folder/$file
      ncpdq -O -h -a -lat $folder/$file $folder/$file
      ncatted -O -a units,time,o,c,"growing seasons since ${YYYY[$i]}-01-01 00:00:00" $folder/$file
    done
  done
done
