#!/bin/bash

PATH=$PATH:/project/joshuaelliott/ggcmi/utils

swift -tc.file tc.data -sites.file midway.xml -config swift.properties downscaler.swift                               \
      -mkfile=/project/ggcmi/AgMIP.output/processed/masks/aggr/gadm0.mask.nc4                                         \
      -reffile=/project/ggcmi/AgMIP.input/other.inputs/reference/faostat/faostat.1961-2012.gadm0.fixed_mirca_mask.nc4 \
      -agglvl=gadm0

if [ $? -eq 0 ]; then
   rm -rf run???
   rm finder.out
fi
