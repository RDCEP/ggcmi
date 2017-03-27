#!/bin/bash

PATH=$PATH:/project/joshuaelliott/ggcmi/utils:$PWD

site=$1
if [ -z "$site" ]; then
    site="sandyb"
fi

swift -tc.file tc.data -sites.file ${site}.xml -config swift.properties modelensemble.swift \
      -w=agmerra,agcfsr,cfsr,erai,grasp,watch,wfdei.cru,wfdei.gpcc,princeton \
      -c=mai,whe,ric,soy,sor,mil \
      -m=rmse

ec=$?
echo Cleaning up, please wait...
sleep 5

if [ "$ec" -eq 0 ]; then
   rm -rf run???
   rm finder.out
fi

echo Done
