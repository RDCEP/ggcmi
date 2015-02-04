#!/bin/bash

range="1980-2009"
mp=true
cr="mean-scale"

for crp in mai whe soy ric sor; do # maize, wheat, etc.
   for g in 240 48 11; do # US, China, and Argentina
      ./multimetrics_plot.py -d /project/ggcmi/AgMIP.output/processed/multimetrics/gadm0/faostat/fixed_mirca_mask -v tscorr -x model -y climate,agcfsr,agmerra,cfsr,erai,grasp,princeton,watch,wfdei.cru,wfdei.gpcc -c gadm0,$g -o scen -o dt -c mp,$mp -c cr,$cr -c time_range,$range -c crop,$crp --outdir /project/ggcmi/AgMIP.output/processed/plots/heatmapMM -f jpg 
   done
done
