#!/bin/bash

for g in 240 48; do # US, china
  ./multimetrics_plot.py -d /project/ggcmi/AgMIP.output/processed/multimetrics/gadm0/faostat -v tscorr -x model -y climate,agcfsr,agmerra,cfsr,erai,grasp,princeton,watch,wfdei.cru,wfdei.gpcc -c gadm0,$g -o scen -o dt -o mp -c cr,mean-scale -c time_range,full -c crop,mai -o . --anon -f eps
done
