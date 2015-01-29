#!/bin/bash

range="1980-2009"
mp=true
cr="mean-scale"

for crp in mai whe soy ric sor; do # maize, wheat, etc.
   for g in 240 48 11; do # US, China, and Argentina
      ./multimetrics_plot.py -d /project/ggcmi/AgMIP.output/processed/multimetrics/gadm0/faostat/fixed -v tscorr -x model -y climate,AgCFSR,AgMERRA,CFSR,ERAI,GRASP,Princeton,WATCH,WFDEI.CRU,WFDEI.GPCC -c gadm0,$g -o scen -o dt -c mp,$mp -c cr,$cr -c time_range,$range -c crop,$crp --outdir /project/ggcmi/AgMIP.output/processed/plots/heatmapMM -f jpg 
   done
done
