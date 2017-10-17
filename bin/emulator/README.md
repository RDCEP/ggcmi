# GGCMI Emulator

# Software Requirements
  - NCO (required for aggregation)
  - netCDF4
  - numpy

# Usage
```
usage: emulator.py -c co2vals -t tempvals -w watervals -n nitrovals -fit fitfile [-agg aggfile] [-aggmask aggmaskfile] [-aggweight aggweightfile] [-aggoutput aggoutputfile] [-o O]
```
| Argument | Description | Required? | Default |
| -------- | ----------- | --------- | ------- |
| --agg | Aggregation level (eg. gadm0) | No | |
| --aggmask | Aggregation mask file | No | gadm.nc |
| --aggweight | Aggregation weight file | No | |
| --aggoutput | Aggregation output filename | No | output.agg.nc |
| --c | Atmospheric co2 concentration | Yes | |
| --t | Temperature anomolies | Yes | |
| --w | Precipitation anomolies | Yes | |
| --n | Fertilizer levels | Yes | |
| --fit | Fit file | Yes | |
| --output | Output file | No | output.nc |

# Specifying C, T, W, N values
The GGCMI Emulator takes input values for C, T, W, and N levels. Users may specify a static value to be used globally, or they may specify a NetCDF file for values that vary spatially. Values must be within the ranges specified below. If a static value is outside the valid range, an error is thrown before the emulator runs. For NetCDF files, values less than the minimum will be set to the minimum, and values greater than the maximum will be set to the maximum.

| Level | Description | Units | Minimum | Maximum | NetCDF Variable |
| ----- | ----------- | ----- | ------- | ------- | --------------- |
| C | Atmospheric CO2 concentration | PPM | 360 | 810 | |
| T | Temperature anomalies | C | -1 | 6 | temperature |
| W | Precipitation anomalies | % | -50 | 30 | precipitation |
| N | Fertilizer levels | kg/ha | 10 | 200 | fertilizer |

If specifying a NetCDF filename, the file must contain a variable with the name listed above. Data must be stored as floats with lat/lon dimensions. For example:

```
float fertilizer(lat, lon)
```
These files must match the lat/lon/resolution of all other files (see 'Latitudes and Longitudes in NetCDF files')

# Fitfile Format
Fit files are the basic input file for the emulator. If the fit file contains information about rainfed simulations, it should be stored in a variable called 'fitparms'. If the fit file contains information about irrigated simulations, it should be stored in a variable called 'fitparmsirr'. At least one of those variable is required in order to run. The first dimension of the parms variable represents a function dimension, the second dimension is latitude, and the third dimension is longitude.

```
float fitparms(funcs=15, lat=360, lon=720);
float fitparmsirr(funcsirr=10, lat=360, lon=720);
```

The function dimension used by fitparms or fitparmsirr must have a property called long_name that is used to store the equation to be evaluated. For example:

```
  long_name = "1,c,t,w,n,c2,t2,w2,n2,ct,cw,cn,tw,tn,wn";
```
In this example, the first dimension is multiplied by 1, the second dimension is multiplied by the user-supplied carbon level, the sixth dimension is multiplied by c^2, etc. These values are then summed to generate data in the output file.

The fit file may also contain metadata used in the generation of the output file. The global attribute "variable" defines the output variable name. The "units" property describes output variable units, and "long_name" describes the output variable long_name. If these are not specified, the output variable will be named "result" and not contain units or long_name.

# Aggregation
Aggregation allows you to to split up your results by region. The emulator includes gadm data for aggregation by country (gadm0), state (gadm1), and county (gadm2). Users may specify this with the --agg option (--agg gadm0, for example). When using gadm, an aggmask is not required as one is included. 

Users may specify their own aggregation file with the --aggmask flag. The data must be in a NetCDF file with the same lat/lon space as the GGCMI Emulations (see section "Latitudes and Longitudes in NetCDF files"). If the data is stored as a shapefile, the gdal_rasterize command may be used to convert to NetCDF. Next, run the included create_agg_limits.py on the nc file to finalize the creation process:

```
$ /path/to/create_agg_limits.py -i yourfile.nc -o aggmask.nc
```

Aggregation weight files are provided in the masks/aggweight directory for barley, cotton, maize, millet, rapeseed, rice, sorghum, soy, and wheat. The 'area' variable has the dimensions irr, lat, and lon. It represents the harvested area for that crop. These area variable must have an irr dimension. The first irr dimension represents rainfed data, and the second irr dimension represents irrigated data.

# Latitudes and Longitudes in NetCDF files
All NetCDF files (input files, mask files, and output files) must use the same spatial extent and resolution. Latitudes and longitudes should be stored as floats in the 'lat' and 'lon' variables. Lats range from 89.75 to -89.75. Lons range from -179.75 to 179.75. The grid spacing is every 0.5 degrees for a total of 360 lats and 720 lons.


