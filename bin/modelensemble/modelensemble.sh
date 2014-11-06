#!/bin/bash

PATH=$PATH:/project/joshuaelliott/ggcmi/utils

swift -tc.file tc.data -sites.file midway.xml -config swift.properties modelensemble.swift \
      -w=agmerra,agcfsr,cfsr,erai,grasp,watch,wfdei.cru,wfdei.gpcc,princeton \
      -c=mai,whe,ric,soy,sor,mil \
      -m=rmse,tscorr

if [ $? -eq 0 ]; then
   rm -rf run???
   rm finder.out
fi
