#!/bin/bash

export PATH=$PATH:$PWD

swift -sites.file midway.xml -tc.file tc.data multimetrics.swift \
      -b=1 \
      -n=1 \
      -d=/project/joshuaelliott/ggcmi/processed/biascorr \
      -m=epic-boku,epic-iiasa,gepic,lpj-guess,lpjml,papsim,pdssat,pegasus \
      -r=/project/joshuaelliott/ggcmi/reference/faostat/faostat.1961-2012.gadm0.nc4 \
      -o=$PWD/processed/multimetrics

# Remove run directories if Swift finishes with no errors
if [ $? -eq 0 ]; then
   echo Removing run directory . . .
   rm -rf run???
fi
