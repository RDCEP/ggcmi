# GGCMI Processed Output Data Pipeline

![](http://users.rcc.uchicago.edu/~joshuaelliott/images/ScreenShot046.png)

A processed, quality-checked, and time-stamped version of the GGCMI phase 1 archive is available [on the DKRZ server](http://vre1.dkrz.de:8080/thredds/catalog/isi_mipEnhanced/AgMIP.output/processed.140515/catalog.html).

##Aggregation
-------------
Using weightings from the MIRCA global land-cover dataset, aggregate all simulated outputs to various different spatial boundaries:
 
 * gadm0: country-level aggregations using the Global Administrative Areas (http://www.gadm.org/)
 * fpu: contains aggregation to food producing units, river-basins and continents
 * kg: contains a koeppen-geiger aggregation based on http://koeppen-geiger.vu-wien.ac.at/present.htm

As an example, the dimensions of the gadm0 files are:

 * dimensions:
 * time = 31 ; // for agmerra which has 31 years
 * scen = 2 ; // lpjml simulates 2 scenarios (default/harmnon)
 * irr = 3 ; // includes rainfed, irrigated, and the weighted sum
 * gadm0_index = 222 ; // includes data for 222 countries

The file includes all the variables present in the original data listed like yield_gadm0, aet_gadm0, ... etc. 
In addition, all aggregation files include an area variable, e.g. area_gadm0(gadm0_index, time, scen, irr), and global aggregates for all variables are also included in all files (e.g. yield_global(time, scen, irr)).

##Detrending and Bias-correction
--------------------------------
The next step in the pipeline is to detrend the simulated outputs (the irr='sum' versions from the gadm0 files) and compare them to the detrended FAOstat data to remove possible biases and add observed country-level trends. The dimensions of the biascorr files are: 

 * dimensions:
 * gadm0_index = 208 ; // includes the 208 countries in the faostat data
 * time = 31 ; // 31 times for agmerra
 * scen = 4 ; // pdssat has 4 scenarios, "default, fullharm, fullharm_pt, harmnon" 
 * detr_methods = 5 ; // detrending methods are "none, lin, quad, ma, ffdtr" 
 * corr_methods = 2 ; // correction methods for now are just "none, VS"

This uses 5 distinct detrending methods directly comparable to the method used on the FAOstat time-series:

* none (i.e. raw data with no detrending)
* linear regression trend removed (lin)
* quadratic regression trend removed (quad)
* moving average trend removed (ma)
* fraction first differences [( y(t)-y(t-1) )/y(t-1)] of the simulation with a linear trend subtracted (ffdtr)

Detrending effectively eliminates any means biases, but it also may be desirable in some cases to remove distributional biases as well. This is accomplished through a bias correction step defined by the corr_methods dimension. Currently this dimension only has 2 values:

1. None (i.e. no attempt is made to remove distributional biases)
2. Variance scaling (VS; the variance of the simulated distribution is rescaled)

The files have 2 variables:

1. yield_sim(gadm0_index, time, scen, detr_methods, corr_methods) includes the simulated data for all countries, times, scenarios, detrending methods, and correction methods.
2. yield_retrend(gadm0_index, time, scen, detr_methods, corr_methods) includes versions of each of the time-series from yield_sim projected onto the trends derived from FAOstat. In each case, the trends are matched so that we combine the linearly detrended simulations with linear trends from FAOstat and quadratic detrended simulations with quadratic trends from FAOstat and etc. Note that for detr_method='none', retrending has no effect and the raw aggregation simply populates through to the yield_retrend variable. Note also that for the retrended versions of the simulated FFDs (ffdtr and ffdtr with VS correction) the retrending procedure is to take the FAO yield at time zero (Y(0)) and define the retrended sim at time T to be Y(0) \Prod_{t=1..T} (1+ffd(t)).

##Multi-metric
--------------
 
This contains one file for each model that has a 7 dimensional matrix of model performance metrics (currently there's just RMSE, the ratio of simulated to observed coefficient-of-variance, and the time-series correlation, all for FAOstat compared against country-level aggregates). The dimensions are inherited from the biascorr files and depend on how many crops, climate products, and scenarios are simulated. For pDSSAT: 

* dimensions:
* gadm0_index = 208 ; // 208 countries
* scen = 4 ; // as above
* detr_methods = 5 ; // as above
* corr_methods = 2 ; // as above
* time_range = 3 ; // 3 sub-sets of the time-series: "full, 1980-2001, 1980-2009" 
* climate = 9 ; // pDSSAT simulates all 9 possible climate datasets
* crop = 6 ; // pDSSAT simulates 6 crops

Notes:
The 3 time ranges are to facilitate comparability among the climate datasets. If a particular analysis is not intended to evaluate differences among climate datasets, then the value 'full' should be used. 1980-2001 is the period that is done by all 9 climate datasets (because WATCH stops in 2001). 1980-2009 was decided as the main 30-year overlap period for climate model comparisons besides WATCH. Actually thus far we're only looking at mai,soy,whe,ric,sor,and mil so the max of crop is 6.

##Multi-model Ensembles
-----------------------

The next step is to produce multi-model ensemble versions of the output to evaluate, for example, how well the ensemble perform relative to individual models. This step uses the multi-metrics files to produce versions of the biascorr files that aggregate all the models into weighted averages. The dimensions of tscorr_watch_.... are:

* dimensions:
* gadm0_index = 208 ; // 208 countries
* time = 41 ; // 41 years of WATCH data
* detr_methods = 5 ; // 5 detrending methods "none,lin,quad,ma,ffdtr" 
* corr_methods = 2 ; // 2 additional correction methods "none,VS" 
* top_models = 7 ; // number of the top performing models to include in the weighting, from 1..Nmods
* weighted = 2 ; // weight by the metric score or just evenly weight all the n top-performing mods.

With top_models=1, this just takes the best model (according to the chosen metric) for each country, and the weighting dimension has not effect in this case. For now we have chosen to leave out model/country/crop combinations with negative correlations, so some countries might not have all 7 values along the top_models dimension. For ones without, a fill_value '1.e+20f' is there instead. Since the scenarios vary so much between models, we don't include scenario as a separate dimension in the ensemble file. each model only contributes once to a given time-series ensemble, and we just take the "best" version of the model according to the given metric to contribute (so the ensembles mix default and fullharm and harmnon and etc. among different models).

##Calculating Metrics from the Ensembles 
---------------------------------------

(also go into the /modelensemble/ folder)
Now we produce multi-metrics files for the ensemble combinations in the same way we did for the models. This allows us to easily include the ensembles in visualization and analysis products produced from the multi-metric files. These ensemble metric files end up being 8 dimensional:
dimensions:
gadm0_index = 208 ;
detr_methods = 5 ;
corr_methods = 2 ;
time_range = 3 ;
climate = 9 ;
crop = 6 ;
top_models = 7 ;
weighted = 2 ;

Notes:
There is no scenario dimension in these files but instead there is a 'top_models' and a 'weighted' dimension that indicates the different properties from the ensembling step.

##Rescaling to a Grid
---------------------
We include this in the version 1 processed folder, but its not actually used at this point for calculating any metrics. The rescaled folder contains versions of the data at the original 0.5 degree resolution that have been "rescaled" using the yield_retrend variable from the biascorr files so that the gridded data reproduces the country level trends from FAOstat. The idea is that this data will be somewhat directly comparable to the Ray et al and Iizumi et al datasets. The dimensions of the rescaled files are
 
* dimensions:
* time = 31 ; // for agmerra
* lat = 360 ;
* lon = 720 ;
* scen = 4 ; // for pdssat
* irr = 3 ; // includes rainfed, irrigated, and weighted sum

The rainfed and irrigated versions of the rescaled files are nonzero everywhere where there were valid simulation values (i.e. they are NOT weighted or masked by any areas). In the 'sum' version of the irr dimension, the value in each grid cell is weighted based on the fraction rainfed and irrigated, and thus gridcells with neither any rainfed or irrigated land for the given crop (according to the MIRCA mask) simply have null values. Eventually we will produce rescaled versions of the weighted ensemble combinations as well, which should prove interesting for a number of purposes.

