#!/bin/bash

PATH=$PATH:/project/joshuaelliott/ggcmi/utils

swift -tc.file tc.data -sites.file midway.xml -config swift.properties downscaler.swift   \
      -mkfile=/project/joshuaelliott/ggcmi/processed/masks/aggr/gadm0.mask.nc4            \
      -reffile=/project/joshuaelliott/ggcmi/reference/faostat/faostat.1961-2012.gadm0.nc4 \
      -agglvl=gadm0

if [ $? -eq 0 ]; then
   rm -rf run???
   rm finder.out
fi
