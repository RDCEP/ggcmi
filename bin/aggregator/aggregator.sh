#!/bin/bash

export PATH=$PATH:$PWD

swift -sites.file midway.xml -tc.file tc.data aggregator.swift \
      -b=1 -n=1 -d=/project/joshuaelliott/ggcmi/AgMIP.output \
      -m=pDSSAT,pAPSIM,LPJ-GUESS,LPJmL,PEGASUS,GEPIC,EPIC-IIASA,EPIC-Boku \
      -w=AgMERRA,AgCFSR,CFSR,ERAI,GRASP,WATCH,WFDEI.CRU,WFDEI.GPCC,Princeton \
      -c=maize,wheat,rice,soy,sorghum,millet \
      -i=/project/joshuaelliott/ggcmi/processed/masks/weight/landuse.ir.nc4 \
      -r=/project/joshuaelliott/ggcmi/processed/masks/weight/landuse.rf.nc4 \
      -a=/project/joshuaelliott/ggcmi/processed/masks/aggr/gadm0.mask.nc4,/project/joshuaelliott/ggcmi/processed/masks/aggr/fpu.mask.nc4,/project/joshuaelliott/ggcmi/processed/masks/aggr/kg.mask.nc4 \
      -t=$PWD/timestamps2.txt \
      -g=/project/joshuaelliott/data/ggcmi/other.inputs/AGMIP_GROWING_SEASON.HARM.version1.23 \
      -o=$PWD/processed/aggs/gadm0,$PWD/processed/aggs/fpu,$PWD/processed/aggs/kg

# Remove run directories if Swift finishes with no errors
if [ $? -eq 0 ]; then
   echo Removing run directory . . .
   rm -rf run???
fi
