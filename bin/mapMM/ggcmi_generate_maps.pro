;this script creates global maps of ggcmi output

pro ggcmi_generate_maps,output_dir,PLOT_METRIC=plot_metric,WEIGHTING_METRIC=weighting_metric,CROPS=crops,CLIMATE_PRODUCTS=climate_products,TIME_RANGE=time_range,DETR_METHODS=detr_methods,CORR_METHODS=corr_methods,WEIGHTED=weighted,TOP_MODELS=top_models,MEAN_PRESERVING=mean_preserving,PROJECTION=projection
;sample runs:
    ;default for rice:  ./ggcmi_generate_maps.sh ./ crops=2
    ;plot rmse for erai climate and crops mai,ric,soy,whe:  ./ggcmi_generate_maps.sh ./ plot_metric=1,climate_products=3,crops=[0,2,4,5]

;plot_metric:      which metric to plot         -->  options: 0=tscorr, 1=rmse, 2=varratio, 3=hitrate, 4=rmse_extreme, 5=bias_extreme  (default: plot_metric=0)
;weighting_metric  metric to use for weighting  -->  options: 0=tscorr, 1=rmse (default: weighting_metric=0 ; for rmse, default is weighting_metric=1)
;crops:            crops to plot                -->  options: 0=mai, 1=mil, 2=ric, 3=sor, 4=soy, 5=whe (default: crops=0)
;climate_products: list of climate products     -->  options: 0=agcfsr, 1=agmerra, 2=cfsr, 3=erai, 4=grasp, 5=princeton, 6=watch, 7=wfdei.cru, 8=wfdei.gpcc  (default: climate_products=[0,1,2,3,4,5,6,7,8])
;time_range:       time range for metrics       -->  options: 0=full, 1=1980-2001, 2=1980-2009  (default: time_range=0)
;detr_methods:     detrending method            -->  options: 0=none, 1=lin, 2=quad, 3=ma, 4=ffdtr  (default: detr_methods=0)
;corr_methods:     variance correction method   -->  options: 0=none, 1=VS, 2=MS                          (default: corr_methods=0)
;weighted:         ensemble weighting method    -->  options: 0=unweighted, 1=tscorr_weighted       (default: weighted=0)
;top_models:       how many models to average   -->  options: 0=highest, 1=top 2 highest,2=top 3 highest, ..., -1=all models   (default:top_models=-1)
;mean_preserving:  preserve the model mean      -->  options: 0=true, 1=false    (default: mean_preserving=0)
;projection:       type of map projection       -->  options: 0=aitoff, 1=albers, 2=azimuthal, 3=conic, 4=cylindrical, 5=gnomic, 6=goodeshomolosine, 7=hammer, 8=lambert, 9=mercator, 10=miller_cylindrical, 11=mollweide, 12=orthographic, 13=robinson, 14=satellite, 15=sinusoidal, 16=stereographic, 17=traverse_mercator (default: projection=13)

;define paths
  path_ggcmi='/project/ggcmi/AgMIP.output/processed/'
  path_plot=output_dir;path_plot='/home/glotter/ag_impacts/plots/ggcmi/maps/'
  
  ;Add the following line to IDL when starting:
    ;!PATH=!PATH+':'+Expand_Path('+/home/glotter/idl_scripts/utilities')

;____define default values___________
  ;define default weighting metric
    if (n_elements(plot_metric) eq 0) then plot_metric= 0 ;'tscorr'
    if (n_elements(weighting_metric) eq 0) then begin
      if (plot_metric eq 1) then weighting_metric=1 ;for 'rmse'
      if (plot_metric ne 1) then weighting_metric=0 ;set all others to 'tscorr'
    endif
  ;define the other variables
    if (n_elements(crops) eq 0) then crops=0 ;['mai','mil','ric','sor','soy','whe']
    if (n_elements(climate_products) eq 0) then climate_products=[0,1,2,3,4,5,6,7,8]  ;['agcfsr','agmerra','cfsr','erai','grasp','princeton','watch','wfdei.cru','wfdei.gpcc']
    if (n_elements(time_range) eq 0) then time_range=0
    if (n_elements(detr_methods) eq 0) then detr_methods=0
    if (n_elements(corr_methods) eq 0) then corr_methods=0
    if (n_elements(weighted) eq 0) then weighted=0
    if (n_elements(top_models) eq 0) then top_models=-1   ;set to -1 instead of 8 because not all files have 8 models included
    if (n_elements(mean_preserving) eq 0) then mean_preserving=0
    if (n_elements(projection) eq 0) then projection=13 ;'robinson'
    
    ;change select variables to strings
       string_arrays={crop:['mai','mil','ric','sor','soy','whe'],clim:['agcfsr','agmerra','cfsr','erai','grasp','princeton','watch','wfdei.cru','wfdei.gpcc'],plot:['tscorr','rmse','varratio','hitrate','rmse_extreme','bias_extreme'],weight:['tscorr','rmse'],proj:['aitoff','albers','azimuthal','conic','cylindrical','gnomic','goodeshomolosine','hammer','lambert','mercator','miller_cylindrical','mollweide','orthographic','robinson','satellite','sinusoidal','stereographic','traverse_mercator']}
       crops= string_arrays.crop(crops)
       climate_products= string_arrays.clim(climate_products)
       plot_metric= string_arrays.plot(plot_metric)
       weighting_metric= string_arrays.weight(weighting_metric)
       projection= string_arrays.proj(projection)
        
;_____________read in country mask_____________________
  id=ncdf_open(path_ggcmi + 'masks/aggr/gadm0.mask.nc4')
  lon_id=ncdf_varid(id,'lon')
  lat_id=ncdf_varid(id,'lat')
  gadm0_mask_id=ncdf_varid(id,'gadm0')
  
  ncdf_varget,id,lon_id,lon
  ncdf_varget,id,lat_id,lat
  ncdf_varget,id,gadm0_mask_id,country_mask
  ncdf_close,id 

  ;______switch missing values to NaN____________
    country_mask=missing_nan(country_mask,-2147483648)


for clim_loop=0,n_elements(climate_products)-1 do begin ;loop over all climate products
for crop_loop=0,n_elements(crops)-1 do begin ;loop over all crops

;____________________read in netcdf file__________________________
  filename=file_search(path_ggcmi + 'multimetrics/gadm0/faostat/' + weighting_metric + '_' + climate_products(clim_loop) + '_hist_' + crops(crop_loop) + '_annual_*.multimetrics.nc4')

  id=ncdf_open(filename)
  country_id=ncdf_varid(id,'gadm0')
  tscorr_id=ncdf_varid(id,'tscorr')
  varratio_id=ncdf_varid(id,'varratio')
  rmse_id=ncdf_varid(id,'rmse')
  hitrate_id=ncdf_varid(id,'hitrate')
  rmse_extreme_id=ncdf_varid(id,'rmse_extreme')
  bias_extreme_id=ncdf_varid(id,'bias_extreme')
  weighted_id=ncdf_varid(id,'wt')
  time_id=ncdf_varid(id,'time_range')
  corr_id=ncdf_varid(id,'cr')
  detr_id=ncdf_varid(id,'dt')
  mean_id=ncdf_varid(id,'mp')

  ncdf_varget,id,country_id,country
  ncdf_varget,id,tscorr_id,tscorr
  ncdf_varget,id,varratio_id,varratio
  ncdf_varget,id,rmse_id,rmse
  ncdf_varget,id,hitrate_id,hitrate
  ncdf_varget,id,rmse_extreme_id,rmse_extreme
  ncdf_varget,id,bias_extreme_id,bias_extreme

  ncdf_attget,id,weighted_id,'long_name',weighted_names_byt
  ncdf_attget,id,time_id,'long_name',time_names_byt
  ncdf_attget,id,corr_id,'long_name',corr_names_byt
  ncdf_attget,id,detr_id,'long_name',detr_names_byt
  ncdf_attget,id,mean_id,'long_name',mean_names_byt
  
  ncdf_close,id
  
  ;______switch missing values to NaN____________
    tscorr=missing_nan(tscorr,1.e20)
    varratio=missing_nan(varratio,1.e20)
    rmse=missing_nan(rmse,1.e20)
    hitrate=missing_nan(hitrate,1.e20)
    rmse_extreme=missing_nan(rmse_extreme,1.e20)
    bias_extreme=missing_nan(bias_extreme,1.e20)

;____define variable names________
  weighted_names=strsplit(string(weighted_names_byt),', ',/extract) ;['unweighted',weighting_metric + 'weighted']
  time_names=strsplit(string(time_names_byt),', ',/extract) ;['full','1980-2001','1980-2009']
  corr_names=strsplit(string(corr_names_byt),', ',/extract) ;['none','VS','MS']
  detr_names=strsplit(string(detr_names_byt),', ',/extract) ;['none','lin','quad','ma','ffdtr']
  mean_names=strsplit(string(mean_names_byt),', ',/extract) ;['true','false']
 
  ;only keep subset that is used for plot
      weighted_names= weighted_names(weighted)
      crop_names= crops  ;already listed as names subset
      climate_names= climate_products  ;already listed as names subset
      time_names= time_names(time_range)
      corr_names= corr_names(corr_methods)
      detr_names= detr_names(detr_methods)
      mean_names= mean_names(mean_preserving)

  ;only keep subset of arrays   (dimensions: (weighted,top_models,corr_methods,detr_methods,time,countries)
      if (plot_metric eq 'tscorr') then metric_long= tscorr(weighted,*,time_range,corr_methods,mean_preserving,detr_methods,*)
      if (plot_metric eq 'varratio') then metric_long= varratio(weighted,*,time_range,corr_methods,mean_preserving,detr_methods,*)
      if (plot_metric eq 'rmse') then metric_long= rmse(weighted,*,time_range,corr_methods,mean_preserving,detr_methods,*)
      if (plot_metric eq 'hitrate') then metric_long= hitrate(weighted,*,time_range,corr_methods,mean_preserving,detr_methods,*)
      if (plot_metric eq 'rmse_extreme') then metric_long= rmse_extreme(weighted,*,time_range,corr_methods,mean_preserving,detr_methods,*)
      if (plot_metric eq 'bias_extreme') then metric_long= bias_extreme(weighted,*,time_range,corr_methods,mean_preserving,detr_methods,*)

   ;mask out top models
     metric=fltarr(n_elements(weighted),n_elements(time_range),n_elements(corr_methods),n_elements(mean_preserving),n_elements(detr_methods),n_elements(country))*!Values.F_NaN
     
       ;find the last index in top_models that is not missing (for each country separately) -- each country/crop/climate has a different number of available crop models
         
         for weight_ind=0,n_elements(weighted)-1 do begin         ;loop over all weights
         for time_ind=0,n_elements(time_range)-1 do begin         ;loop over all time ranges
         for corr_ind=0,n_elements(corr_methods)-1 do begin       ;loop over all correction methods
         for mean_ind=0,n_elements(mean_preserving)-1 do begin    ;loop over all mean preserving options
         for detr_ind=0,n_elements(detr_methods)-1 do begin       ;loop over all detrending methods
         for count_ind=0,n_elements(country)-1 do begin           ;loop over all countries
          
           top_ind=-1
           while ((finite(metric_long(weight_ind,top_ind,time_ind,corr_ind,mean_ind,detr_ind,count_ind)) eq 0) and (abs(top_ind) lt n_elements(metric_long(0,*,0,0,0,0,0)))) do top_ind--     ;keep stepping back until yield is not missing
           
           if (top_models ne -1) then top_ind= top_models   ;set top_ind to top_models if top_models is not -1

           metric(weight_ind,time_ind,corr_ind,mean_ind,detr_ind,count_ind)= metric_long(weight_ind,top_ind,time_ind,corr_ind,mean_ind,detr_ind,count_ind)

         endfor   ;stop looping over countries
         endfor   ;stop looping over detrending methods
         endfor   ;stop looping over mean preserving options
         endfor   ;stop looping over correction methods
         endfor   ;stop looping over time ranges
         endfor   ;stop looping over weights


;_____flatten country data to an array with elements (countries x scenarios)____________
   if (n_elements(metric_flat) eq 0) then metric_flat=[[]]
   if (n_elements(plot_names) eq 0) then plot_names=[[]]
 
     for weight_ind=0,n_elements(weighted)-1 do begin         ;loop over all weights
     for time_ind=0,n_elements(time_range)-1 do begin         ;loop over all time ranges
     for corr_ind=0,n_elements(corr_methods)-1 do begin       ;loop over all correction methods
     for mean_ind=0,n_elements(mean_preserving)-1 do begin    ;loop over all mean preserving options
     for detr_ind=0,n_elements(detr_methods)-1 do begin       ;loop over all detrending methods
            
            metric_flat= [[metric_flat],[reform(metric(weight_ind,time_ind,corr_ind,mean_ind,detr_ind,*))]]    ;cols=countries, rows=panel plots
    
            ;define name for plot title
              plot_name=''
              
              ;define plot names by the variable that is changing (i.e., has multiple elements)
                if (n_elements(climate_names) gt 1) then plot_name= plot_name + ' ' + climate_names(clim_loop)
                if (n_elements(crop_names) gt 1) then plot_name= plot_name + ' ' + crop_names(crop_loop)
                if (n_elements(weighted_names) gt 1) then plot_name= plot_name + ' ' + weighted_names(weight_ind)
                if (n_elements(time_names) gt 1) then plot_name= plot_name + ' ' + time_names(time_ind)
                if (n_elements(corr_names) gt 1) then plot_name= plot_name + ' ' + corr_names(corr_ind)
                if (n_elements(mean_names) gt 1) then plot_name= plot_name + ' ' + mean_names(mean_ind)
                if (n_elements(detr_names) gt 1) then plot_name= plot_name + ' ' + detr_names(detr_ind)
              ;if no element is changing (i.e., only one plot is being produced) then use the climate and crop name
                if (n_elements(plot_name) eq 0) then plot_name= [climate_names(clim_loop) + ' ' + crop_names(crop_loop)]
             
              plot_names= [plot_names,plot_name] 
    
     endfor   ;stop looping over detrending methods
     endfor   ;stop looping over mean preserving options
     endfor   ;stop looping over correction methods
     endfor   ;stop looping over time ranges
     endfor   ;stop looping over weights    

endfor   ;stop looping over crops
endfor   ;stop looping over climates            
          
;_____change country data into lat/lon array____________   
   metric_map= fltarr(n_elements(lon),n_elements(lat),n_elements(metric_flat(0,*)))*!Values.F_NaN    ;lon x lat x num_panels
   for panel=0,n_elements(metric_flat(0,*))-1 do begin   ;loop over all panels
     metric_map_temp= metric_map(*,*,panel)
    
     for k=0,n_elements(metric_flat(*,0))-1 do begin   ;loop over all countries
       w_country= where(country_mask eq country(k))
       metric_map_temp(w_country)= metric_flat(k,panel)
     endfor  
     
     metric_map(*,*,panel)= metric_map_temp

   endfor

;____________________plot_____________________________
;define number of panels
  num_panels= n_elements(metric_flat(0,*))
  num_panel_cols= ceil(sqrt(num_panels))   ;round up so all panels fit
  num_panel_rows= ceil(num_panels/num_panel_cols)

;define file name parameters________________________
  crop_print=''
  climate_print=''
  top_print=''
  weight_print=''
  corr_print=''
  mean_print=''
  detr_print=''
  time_print=''
  if (n_elements(crops) ne 6) then begin
    for q=0,n_elements(crop_names)-1 do begin
      crop_print= crop_print + '_' + crop_names(q)
    endfor
  endif else crop_print= '_all'
  
  if (n_elements(climate_products) ne 9) then begin
    for q=0,n_elements(climate_names)-1 do begin
      climate_print= climate_print + '_' + climate_names(q)
    endfor
  endif else climate_print= '_all'
  
  if (top_models ne -1) then top_print='.top_' + strtrim(top_models,2)
  
  if (weighted ne 0) then begin
    weight_print='.weight'
    for q=0,n_elements(weighted_names)-1 do begin
       weight_print= weight_print + '_' + weighted_names(q)
    endfor
  endif
  
  if (corr_methods ne 0) then begin
    corr_print='.corr'
    for q=0,n_elements(corr_names)-1 do begin
      corr_print= corr_print + '_' + corr_names(q)
    endfor
  endif

  if (mean_preserving ne 0) then begin
    mean_print='.mean'
    for q=0,n_elements(mean_names)-1 do begin
      mean_print= mean_print + '_' + mean_names(q)
    endfor
  endif
  
  if (detr_methods ne 0) then begin
    detr_print='.detr'
    for q=0,n_elements(detr_names)-1 do begin
      detr_print= detr_print + '_' + detr_names(q)
    endfor
  endif
  
  if (time_range ne 0) then begin
    time_print='.time'
    for q=0,n_elements(time_names)-1 do begin
      time_print= time_print + '_' + time_names(q)
    endfor
  endif
;____________________________________

set_plot,'ps'                  ;use to save as postscript file
device,filename=path_plot + plot_metric + '.crop' + crop_print + '.climate' + climate_print + top_print + weight_print + corr_print + mean_print + detr_print + time_print + '.eps'    ;set path to save postscript file
device,/encapsul,/color,set_font='Helvetica',/TT_Font,xsize=25,ysize=20
device,decomposed=0        ;set device to get color tables to work
device,/ISOLATIN1    ;use to get symbols


;Define image settings
!Y.OMargin=[6,0]            ;Set margins for title, colorbar
;!X.OMargin=[0,4]            ;Set margins on sides as well
!P.Multi=[0,num_panel_cols,num_panel_rows,0,0]        ;use for multiple panel plots
!P.background=255           ;Set background to white
!P.Color=0                  ;Set drawn lines to black
!P.Font=1                   ;Set font to TrueType fonts (as set in device)

numlevels=8

;set the levels  
  if (plot_metric eq 'tscorr') then begin
    max_data=1.
    min_data=-1.
  endif else begin
    max_data=max(metric_map,/nan)
    min_data=min(metric_map,/nan)
  endelse

  step=(max_data-min_data)/numlevels
  userlevels=indgen(numlevels)*step + min_data
  usercolors=indgen(numlevels)+1

;______load color table and levels for num days tmax___________________________________
;flop lat array so its increaing
  lat= reverse(lat)
  metric_map= reverse(metric_map,2)

  color_start=1.
  cgloadct,0,NCOLORS=1   ;set gray to level 1
  cgloadct,70,/REVERSE,BOTTOM=color_start,NCOLORS=numlevels
  
  for cells=0,num_panels-1 do begin 

    map_set,name=projection,/grid,/noborder,limit=[-90,-180,90,180],/advance,/isotropic,title=plot_names(cells),charsize=2.,/noerase,color=cgcolor('black')
    map_continents,/FILL_CONTINENTS,/overplot,color=cgcolor('white')
    contour,metric_map(*,*,cells),lon,lat,levels=userlevels,c_colors=usercolors,/overplot,/cell_fill
;    map_byt=bytscl(metric_map(*,*,cells),max=max_data,min=min_data,/nan,top=numlevels-1)+1B   ;bytscl(metric_map(*,*,cells))*numlevels
;    w_miss=where(finite(metric_map(*,*,cells)) eq 0)
;    map_byt(w_miss)= 10   ;set missing values to 10

;    result= Map_Image(map_byt, xstart, ystart, xsize, ysize, LatMin=-90, LatMax=90, LonMin=-180, LonMax=180,max_value=9, missing=cgcolor('white'), Compress=1)
;    TV, result, xstart, ystart, XSize=xsize, YSize=ysize
    ;result= map_patch(map_byt,lon,lat, xstart=x0, ystart=y0, xsize=x1, ysize=y1, Lat0=-90, Lat1=90, Lon0=-180, Lon1=180,max_value=9, missing=cgcolor('lightgray'))   
    ;tv,result,x0,y0,xsize=x1,ysize=y1
    map_continents,/COUNTRIES,/overplot,color=cgcolor('darkgray')
    map_continents,/overplot,color=cgcolor('darkgray')

  endfor

  ;Create Colorbar
    if (plot_metric eq 'tscorr') then colorbar_title= 'correlation'
    if (plot_metric eq 'rmse') then colorbar_title= 'rmse (tonnes/ha)'
    if (plot_metric eq 'varratio') then colorbar_title= 'varratio'
    if (plot_metric eq 'hitrate') then colorbar_title= 'hitrate'
    if (plot_metric eq 'rmse_extreme') then colorbar_title= 'rmse_extreme'
    if (plot_metric eq 'bias_extreme') then colorbar_title= 'bias_extreme'
    
    cgcolorbar,range=[min_data,max_data],position=[0.35,0.09,0.65,0.11],charsize=1.5,/discrete,ncolors=numlevels,divisions=6,bottom=1,format='(F0.1)',title=colorbar_title,color=cgcolor('black')

device,/close


!Y.OMargin=[0,0]
!X.OMargin=[0,0]
!P.Multi=[0,0,0,0,0]
!P.background=0           ;Set background to white
!P.Color=255              ;Set drawn lines to black
!P.Font=-1                ;Set back to default Hershey vector-drawn fonts
set_plot,'x'


end