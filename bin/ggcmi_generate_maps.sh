#!/bin/bash

#___sample runs_________________ 
#       default for rice:  ./ggcmi_generate_maps.sh ./ crops=2
#       plot rmse for erai climate and crops mai,ric,soy,whe:  ./ggcmi_generate_maps.sh ./ plot_metric='tscorr',climate_products=3,crops=[0,2,4,5]
	    
outputdir=$1
options=$2   #OK to not use options

if [ $outputdir == "./" ]; then outputdir="$(pwd)""/"
fi


#____run IDL__________________
echo "running IDL"

cd /home/glotter/idl_scripts/ggcmi/  #change path to that of idl script

if [ -n "$options" ]; then
    /software/exelis/idl82/bin/idl -e "ggcmi_climate_maps,'$outputdir',$options" &> /dev/null
else
    /software/exelis/idl82/bin/idl -e "ggcmi_climate_maps,'$outputdir'" &> /dev/null
fi


#____convert to jpg___________
echo "printing jpg"

eps_file=$(find $outputdir -mtime -0.25 -type f -name '*.eps')    #find eps file generated in the last 15 seconds
convert -density 300 $eps_file "$eps_file"".jpg"
rm $eps_file

 
#_____options______________________________________________ 
#   plot_metric:      which metric to plot         -->  options: 'tscorr','rmse', 'varratio' (default='tscorr')
#   weighting_metric  metric to use for weighting  -->  options: 'tscorr','rmse' (default: weighting_metric= plot_metric -- for varratio, default is 'tscorr')
#   crops:            crops to plot                -->  options: 0=mai, 1=mil, 2=ric, 3=sor, 4=soy, 5=whe (default: crops=0)
#   climate_products: list of climate products     -->  options: 0=agcfsr, 1=agmerra, 2=cfsr, 3=erai, 4=grasp, 5=princeton, 6=watch, 7=wfdei.cru, 8=wfdei.gpcc  (default: climate_products=[0,1,2,3,4,5,6,7,8])
#   time_range:       time range for metrics       -->  options: 0=full, 1=1980-2001, 2=1980-2009  (default: time_range=0)
#   detr_methods:     detrending method            -->  options: 0=none, 1=lin, 2=quad, 3=ma, 4=ffdtr  (default: detr_methods=0)
#   corr_methods:     variance correction method   -->  options: 0=none, 1=VS                          (default: corr_methods=0)
#   weighted:         ensemble weighting method    -->  options: 0=unweighted, 1=tscorr_weighted       (default: weighted=0)
#   top_models:       how many models to average   -->  options: 0=highest, 1=top 2 highest,2=top 3 highest, ..., -1=all models   (default:top_models=-1)
#   projection:       type of map projection       -->  options: 'aitoff', 'albers', 'azimuthal', 'conic', 'cylindrical', 'gnomic', 'goodeshomolosine', 'hammer', 'lambert, 'mercator', 'miller_cylindrical', 'mollweide', 'orthographic', 'robinson', 'satellite', 'sinusoidal', 'stereographic', 'traverse_mercator' (default: 'robinson')

