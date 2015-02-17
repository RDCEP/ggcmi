#!/bin/bash

range="1980-2009"
mp=true
cr="mean-scale"
crp=mai
g=240 # US

./multimetrics_plot.py -d /project/ggcmi/AgMIP.output/processed/multimetrics/gadm0/faostat/fixed_mirca_mask -m tscorr -x model -y climate:agcfsr,agmerra,cfsr,erai,grasp,princeton,watch,wfdei.cru,wfdei.gpcc -o scen -c gadm0:$g -c dt:quad -c mp:$mp -c cr:$cr -c time_range:$range -c crop:$crp --outdir . --nm 1 --wt 2 -f png
