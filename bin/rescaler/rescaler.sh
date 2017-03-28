#!/bin/bash

utils=$( readlink -f ../../utils )
PATH=$PATH:$utils:$PWD
mkfile="/project/ggcmi/phase1.final/processed/masks/aggr/gadm0.mask.nc4"

site=$1
if [ -z "$site" ]; then
    site="sandyb"
fi

swift -tc.file tc.data -sites.file ${site}.xml -config swift.properties rescaler.swift \
      -mkfile=$mkfile \
      -crmthd=2,0,1 \
      -agglvl=gadm0

rc=$?
echo Cleaning up, please wait
sleep 5

if [ "$rc" -eq 0 ]; then
   rm -rf run???
   rm finder.out
fi
