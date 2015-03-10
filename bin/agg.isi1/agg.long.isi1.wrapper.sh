#!/bin/bash

export PATH=$PATH:$PWD

swift -sites.file midway.xml -tc.file tc.data agg.long.isi1.swift       \
      -models=pDSSAT,GEPIC,EPIC,IMAGE_LEITAP,LPJ-GUESS,LPJmL,PEGASUS    \
      -gcms=GFDL-ESM2M,HadGEM2-ES,IPSL-CM5A-LR,MIROC-ESM-CHEM,NorESM1-M \
      -crops=maize,wheat,soy,rice                                       \
      -co2s=co2,noco2                                                   \
      -rcps=rcp2p6,rcp4p5,rcp6p0,rcp8p5                                 \
      -vars=yield,plant-day,maty-day

# Remove run directories if Swift finishes with no errors
if [ $? -eq 0 ]; then
   echo Removing run directory . . .
   rm -rf run???
fi
