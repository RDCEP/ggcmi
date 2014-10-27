#!/bin/bash

# plot u.s.
./multimetrics_plot.py -d /project/joshuaelliott/ggcmi/processed/multimetrics -v tscorr -x model -y climate,agcfsr,agmerra,cfsr,erai,grasp,princeton,watch,wfdei.cru,wfdei.gpcc -c gadm0_index,240 -o scen -c corr_methods,VS -c time_range,full -c crop,mai -o detr_methods

# plot china
./multimetrics_plot.py -d /project/joshuaelliott/ggcmi/processed/multimetrics -v tscorr -x model -y climate,agcfsr,agmerra,cfsr,erai,grasp,princeton,watch,wfdei.cru,wfdei.gpcc -c gadm0_index,48 -o scen -c corr_methods,VS -c time_range,full -c crop,mai -o detr_methods

# plot brazil
./multimetrics_plot.py -d /project/joshuaelliott/ggcmi/processed/multimetrics -v tscorr -x model -y climate,agcfsr,agmerra,cfsr,erai,grasp,princeton,watch,wfdei.cru,wfdei.gpcc -c gadm0_index,32 -o scen -c corr_methods,VS -c time_range,full -c crop,mai -o detr_methods

# plot india
./multimetrics_plot.py -d /project/joshuaelliott/ggcmi/processed/multimetrics -v tscorr -x model -y climate,agcfsr,agmerra,cfsr,erai,grasp,princeton,watch,wfdei.cru,wfdei.gpcc -c gadm0_index,105 -o scen -c corr_methods,VS -c time_range,full -c crop,mai -o detr_methods

# plot mexico
./multimetrics_plot.py -d /project/joshuaelliott/ggcmi/processed/multimetrics -v tscorr -x model -y climate,agcfsr,agmerra,cfsr,erai,grasp,princeton,watch,wfdei.cru,wfdei.gpcc -c gadm0_index,145 -o scen -c corr_methods,VS -c time_range,full -c crop,mai -o detr_methods

# average over top maize production areas
./multimetrics_plot.py -d /project/joshuaelliott/ggcmi/processed/multimetrics -v tscorr -x model -y climate,agcfsr,agmerra,cfsr,erai,grasp,princeton,watch,wfdei.cru,wfdei.gpcc -a gadm0_index,240,48,32,105,145,163,106,225,11,237 -o scen -c corr_methods,VS -c time_range,full -c crop,mai -o detr_methods -w /project/joshuaelliott/ggcmi/reference/faostat/weights.maize.nc4

# move to plot directory
mv *.png /project/joshuaelliott/ggcmi/processed/multimetrics_plots
