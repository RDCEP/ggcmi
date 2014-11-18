#!/bin/bash

export PATH=$PATH:$PWD

swift -sites.file midway.xml -tc.file tc.data aggregator.swift \
      -b=1 -n=1 -d=/project/ggcmi/AgMIP.output \
      -m=pDSSAT,pAPSIM,LPJ-GUESS,LPJmL,PEGASUS,GEPIC,EPIC-IIASA,EPIC-Boku \
      -w=AgMERRA,AgCFSR,CFSR,ERAI,GRASP,WATCH,WFDEI.CRU,WFDEI.GPCC,Princeton \
      -c=maize,wheat,rice,soy,sorghum,millet \
      -i=/project/ggcmi/AgMIP.output/processed/masks/weight/landuse.ir.nc4 \
      -r=/project/ggcmi/AgMIP.output/processed/masks/weight/landuse.rf.nc4 \
      -a=/project/ggcmi/AgMIP.output/processed/masks/aggr/gadm0.mask.nc4,/project/ggcmi/AgMIP.output/processed/masks/aggr/fpu.mask.nc4,/project/ggcmi/AgMIP.output/processed/masks/aggr/kg.mask.nc4 \
      -t=/project/joshuaelliott/ggcmi/bin/aggregator/timestamps.txt \
      -g=/project/ggcmi/AgMIP.input/other.inputs/AGMIP_GROWING_SEASON.HARM.version1.23 \
      -o=/project/ggcmi/AgMIP.output/processed/aggs/gadm0,/project/ggcmi/AgMIP.output/processed/aggs/fpu,/project/ggcmi/AgMIP.output/processed/aggs/kg

# Remove run directories if Swift finishes with no errors
if [ $? -eq 0 ]; then
   echo Removing run directory . . .
   rm -rf run???
fi
