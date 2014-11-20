#!/bin/bash

PATH=$PATH:/project/joshuaelliott/ggcmi/utils

swift -tc.file tc.data -sites.file midway.xml -config swift.properties rescaler.swift \
      -mkfile=/project/ggcmi/AgMIP.output/processed/masks/aggr/gadm0.mask.nc4         \
      -crmthd=2,0,2                                                                   \
      -agglvl=gadm0

if [ $? -eq 0 ]; then
   rm -rf run???
   rm finder.out
fi
