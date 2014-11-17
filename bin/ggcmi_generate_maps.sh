#!/bin/bash

#___sample runs_________________ 
#       default for rice:  ./ggcmi_generate_maps.sh ./ crops=2
#       plot rmse for erai climate and crops mai,ric,soy,whe:  ./ggcmi_generate_maps.sh ./ plot_metric=1,climate_products=3,crops=[0,2,4,5]
	    
outputdir=$1
options=$2   #OK to not use options

if [ $outputdir == "./" ]; then outputdir="$(pwd)""/"
fi


#____run IDL__________________
echo "running IDL"

if [ -n "$options" ]; then
    /software/exelis/idl82/bin/idl -e "ggcmi_generate_maps,'$outputdir',$options" &> /dev/null
else
    /software/exelis/idl82/bin/idl -e "ggcmi_generate_maps,'$outputdir'" &> /dev/null
fi


#____convert to jpg___________
echo "printing jpg"

eps_file=$(find $outputdir -mtime -0.25 -type f -name '*.eps')    #find eps file generated in the last 15 seconds
convert -density 300 $eps_file "$eps_file"".jpg"
rm $eps_file

 
#_____options______________________________________________ 
# plot_metric:      which metric to plot         -->  options: 0=tscorr, 1=rmse, 2=varratio, 3=hitrate, 4=rmse_extreme, 5=bias_extreme  (default: plot_metric=0)
# weighting_metric  metric to use for weighting  -->  options: 0=tscorr, 1=rmse (default: weighting_metric=0 ; for rmse, default is weighting_metric=1)
# crops:            crops to plot                -->  options: 0=mai, 1=mil, 2=ric, 3=sor, 4=soy, 5=whe (default: crops=0)
# climate_products: list of climate products     -->  options: 0=agcfsr, 1=agmerra, 2=cfsr, 3=erai, 4=grasp, 5=princeton, 6=watch, 7=wfdei.cru, 8=wfdei.gpcc  (default: climate_products=[0,1,2,3,4,5,6,7,8])
# time_range:       time range for metrics       -->  options: 0=full, 1=1980-2001, 2=1980-2009  (default: time_range=0)
# detr_methods:     detrending method            -->  options: 0=none, 1=lin, 2=quad, 3=ma, 4=ffdtr  (default: detr_methods=0)
# corr_methods:     variance correction method   -->  options: 0=none, 1=VS, 2=MS                          (default: corr_methods=0)
# weighted:         ensemble weighting method    -->  options: 0=unweighted, 1=tscorr_weighted       (default: weighted=0)
# top_models:       how many models to average   -->  options: 0=highest, 1=top 2 highest,2=top 3 highest, ..., -1=all models   (default:top_models=-1)
# mean_preserving:  preserve the model mean      -->  options: 0=true, 1=false    (default: mean_preserving=0)
# projection:       type of map projection       -->  options: 0=aitoff, 1=albers, 2=azimuthal, 3=conic, 4=cylindrical, 5=gnomic, 6=goodeshomolosine, 7=hammer, 8=lambert, 9=mercator, 10=miller_cylindrical, 11=mollweide, 12=orthographic, 13=robinson, 14=satellite, 15=sinusoidal, 16=stereographic, 17=traverse_mercator (default: projection=13)

