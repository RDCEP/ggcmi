#!/usr/bin/env python

# import modules
from os import listdir
from os.path import sep, exists
from optparse import OptionParser
from netCDF4 import Dataset as nc
from numpy.ma import isMaskedArray
from numpy import inf, float32, zeros, ones, where, diff, array_equal, logical_and

def climate_years(clim_names):
    clim_names_list = ['AgCFSR', 'AgMERRA', 'CFSR', 'ERAI', 'GRASP', \
                       'Princeton', 'WATCH', 'WFDEI.CRU', 'WFDEI.GPCC']
    year_range_list = [[1980, 2010], [1980, 2010], [1980, 2010], [1979, 2010], [1961, 2010], \
                       [1948, 2008], [1958, 2001], [1979, 2009], [1979, 2009]]
    yrs = [0] * len(clim_names)
    for i in range(len(clim_names)):
        if clim_names[i] in clim_names_list:
            yrs[i] = year_range_list[clim_names_list.index(clim_names[i])]
        else:
            raise Exception('Unknown climate')
    return yrs

def abr_crop_names(crop_names):
    full_names_list = ['maize', 'wheat', 'soy', 'rice', 'sorghum', 'millet', \
                       'managed_grass', 'sugarcane', 'barley', 'oat', 'rapeseed', \
                       'rye', 'sugar_beet', 'sunflower']
    abr_names_list = ['mai', 'whe', 'soy', 'ric', 'sor', 'mil', \
                      'mgr', 'sug', 'bar', 'oat', 'rap', \
                      'rye', 'sgb', 'sun']
    abr_names = [0] * len(crop_names)
    for i in range(len(crop_names)):
        if crop_names[i] in full_names_list:
            abr_names[i] = abr_names_list[full_names_list.index(crop_names[i])]
        else:
            raise Exception('Unknown crop')
    return abr_names

# parse inputs
parser = OptionParser()
parser.add_option("-d", "--dir", dest = "dir", default = "AgMIP.output", type = "string",
                  help = "Directory in which to perform verification")
parser.add_option("-m", "--mod", dest = "mod", default = "pDSSAT,pAPSIM", type = "string",
                  help = "Comma-separated list of crop models to verify (* = all models)")
parser.add_option("-w", "--weath", dest = "weath", default = "AgCFSR,AgMERRA,CFSR,ERAI,GRASP,Princeton,WATCH,WFDEI.CRU,WFDEI.GPCC", type = "string",
                  help = "Comma-separated list of weather datasets to verify (* = all weather datasets)")
parser.add_option("-c", "--crop", dest = "crop", default = "maize", type = "string",
                  help = "Comma-separated list of crops to verify (* = all crops, except 'others')")
parser.add_option("-s", "--sumdir", dest = "sumdir", default = "", type = "string",
                  help = "Where to save summary report file(s)")
options, args = parser.parse_args()

dir = options.dir
if dir[-1] == sep: dir = dir[: -1] # remove final separator

# crop models
if options.mod == '*':
    dirs = listdir(dir) # process all models
    if 'upload_stats' in dirs: dirs.remove('upload_stats')
    if 'aggregations' in dirs: dirs.remove('aggregations')
else:
    dirs = options.mod.split(',')

# climate data
if options.weath == '*':
    climmodels = listdir(dir + sep + dirs[0]) # process all climate data, based on first directory
else:
    climmodels = options.weath.split(',')
climyears = climate_years(climmodels)
numclims = len(climmodels)

# crops
if options.crop == '*':
    crops = listdir(dir + sep + dirs[0] + sep + climmodels[0])
    if 'others' in crops: crops.remove('others')
else:
    crops = options.crop.split(',')
cropabbr = abr_crop_names(crops)
numcrops = len(crops)

# scenarios
scens = ['fullharm_noirr', 'fullharm_firr', 'default_noirr', 'default_firr', 'harmnon_noirr', 'harmnon_firr']
dssat_scens = scens + ['fullharm_noirr_pt', 'fullharm_firr_pt'] # scenarios for pDSSAT and EPIC models
epic_scens = scens + \
             ['fullharm_noirr_pm', 'fullharm_firr_pm'] + \
             ['fullharm_noirr_pt', 'fullharm_firr_pt'] + \
             ['fullharm_noirr_hg', 'fullharm_firr_hg'] + \
             ['fullharm_noirr_br', 'fullharm_firr_br']

# variables
vlist = ['yield', 'pirrww', \
         'plant-day', 'maty-day', 'aet', 'initr', 'gsprcp', \
         'anth-day', 'biom', 'leach', 'gsrsds', 'sumt', 'sco2', 'sn2o']
vranges = [[0, 50], [0, 10000], \
           [1, 366], [1, 366], [0, 10000], [0, 1000], [0, 10000], \
           [1, 366], [0, 100], [0, 1000], [0, 100000], [-inf, inf], [-inf, inf], [-inf, inf]]
vunits = ['t ha-1 yr-1', 'mm yr-1', \
          'day of year', 'days from planting', 'mm yr-1', 'kg ha-1 yr-1', 'mm yr-1', \
          'days from planting', 't ha-1 yr-1', 'kg ha-1 yr-1', 'w m-2 yr-1', 'deg C-days yr-1', 'kg C ha-1', 'kg N2O-N ha-1']
numtotvars = len(vlist)
mandatoryvars = zeros((numtotvars,), dtype = bool)
mandatoryvars[: 2] = [True, True]

# indices for parsing filenames
fileindices = {'scen': [3, 4], 'var': 5, 'crop': 6, \
               'annual': 7, 'yr0': 8, 'yr1': 9}
fileindices2 = {'scen': [3, 4, 5], 'var': 6, 'crop': 7, \
                'annual': 8, 'yr0': 9, 'yr1': 10}

# iterate through all permutations
for d in dirs:
    msg_sims = '' # missing simulations
    msg_data = '' # missing data within a simulation
    changes = ''  # changes to make
    errors = ''   # errors in data
    info = ''     # basic file info
    fatal = ''    # fatal errors
    
    if 'pDSSAT' in d:
        scenarios = dssat_scens
    elif 'EPIC' in d:
        scenarios = epic_scens
    else:
        scenarios = scens

    # for marking missing data
    msg_datmat = ones((numclims, numcrops, numtotvars, len(scenarios)), dtype = bool)

    for cm in range(len(climmodels)):
        for cp in range(len(crops)):
            sim = '{:15s}{:15s}{:15s}'.format(d, climmodels[cm], crops[cp])
            subdir = dir + sep + d + sep + climmodels[cm] + sep + crops[cp]
            print 'Processing directory', subdir, '. . .'
            
            # check if simulation exists
            if not exists(subdir) or not len(listdir(subdir)):
                msg_sims += sim + '\n'
                continue
            
            # iterate over all files
            files = listdir(subdir)
            for f in files:
                fileparts = f.split('_')
                if len(fileparts) == 10:
                    fi = fileindices
                elif len(fileparts) == 11:
                    fi = fileindices2
                else:
                    continue # invalid file format
                scenidx = fi['scen']; varidx = fi['var']
                cropidx = fi['crop']; annidx = fi['annual']
                yr0idx = fi['yr0'];   yr1idx = fi['yr1']
                
                # for basic file info
                info_datmat = zeros((7,), dtype = '|S128')
                info_datmat[:] = 'N/A'
                info_datmat[0] = f
                
                # variable
                varname = fileparts[varidx].lower()
                if varname in vlist:
                    varfullidx = vlist.index(varname)
                else:
                    continue # cannot find variable
                
                # scenario
                scen = '_'.join([fileparts[i].lower() for i in scenidx]) # combine
                if scen in scenarios:
                    scenfullidx = scenarios.index(scen)
                else:
                    continue # cannot find scenario
                scen = scen.split('_') # resplit
                
                # mark as found
                msg_datmat[cm, cp, varfullidx, scenfullidx] = False
                sim2 = sim + '{:20s}{:15s}'.format(scenarios[scenfullidx], vlist[varfullidx])
                
                # change filename as necessary
                if fileparts[0] != d.lower():
                    fileparts[0] = d.lower() # change model name
                if fileparts[1] != climmodels[cm].lower():
                    fileparts[1] = climmodels[cm].lower() # change climate name
                if fileparts[2] != 'hist':
                    fileparts[2] = 'hist' # change to hist
                for i in range(len(scen)):
                    if fileparts[scenidx[i]] != scen[i]:
                        fileparts[scenidx[i]] = scen[i] # change scenario
                if fileparts[varidx] != varname:
                    fileparts[varidx] = varname # change var to lowercase
                if fileparts[cropidx] != cropabbr[cp]:
                    if fileparts[cropidx] != 'ri2' or cropabbr[cp] != 'ric': # special case of second season for rice called ri2
                        fileparts[cropidx] = cropabbr[cp] # change crop name
                if fileparts[annidx] != 'annual':
                    fileparts[annidx] = 'annual' # change to annual
                    
                yr0f = int(fileparts[yr0idx]); yr1f = int(fileparts[yr1idx].split('.')[0])
                info_datmat[5] = '[' + str(yr0f) + ', ' + str(yr1f) + ']' # save year range
                
                ncf = nc(subdir + sep + f) # load file
		fvars = ncf.variables.keys()
                
                lat = None; lon = None; time = None
                if not 'lat' in fvars: # check if latitude exists
                    if 'latitude' in fvars:
                        changes += sim2 + 'DIMENSION: Change latitude to lat\n'
                        lat = ncf.variables['latitude'][:]
                    else:
                        errors += sim2 + 'DIMENSION: No latitude dimension in file\n'
                else:
                    lat = ncf.variables['lat'][:]
                if not 'lon' in fvars: # check if longitude exists
                    if 'longitude' in fvars:
                        changes += sim2 + 'DIMENSION: Change longitude to lon\n'
                        lon = ncf.variables['longitude'][:]
                    else:
                        errors += sim2 + 'DIMENSION: No longitude dimension in file\n'
                else:
                    lon = ncf.variables['lon'][:]
                if not 'time' in fvars: # check if time exists
                    errors += sim2 + 'DIMENSION: No time dimension in file\n'
                else:
                    time = ncf.variables['time'][:]
                    info_datmat[6] = len(time)
                
                if not lat is None:
                    if any(diff(lat) > 0):
                        changes += sim2 + 'VALUES: Change latitudes to descending order\n'
                    if abs(lat[0] - 90) > 1. or abs(lat[-1] + 90) > 1.:
                        changes += sim2 + 'RANGE: Change latitude range to [90, -90]\n'
                if not lon is None:
                    if any(diff(lon) < 0):
                        changes += sim2 + 'VALUES: Change longitudes to ascending order\n'    
                    if abs(lon[0] + 180) > 1. or abs(lon[-1] - 180) > 1.:
                        changes += sim2 + 'RANGE: Change longitude range to [-180, 180]\n'
                if not time is None:
                    if not array_equal(time, range(1, len(time) + 1)):
                        changes += sim2 + 'VALUES: Change time values to increase uniformly from one\n'
                    timev = ncf.variables['time']
                    if not 'units' in timev.ncattrs():
                        changes += sim2 + 'UNITS: Add units to time\n'
                    else:
                        try:
                            tsplit = timev.units.split('growing seasons since ')[1].split(' ')
                            yr0, mth0, day0 = tsplit[0].split('-')[0 : 3]
                            hr0, min0, sec0 = tsplit[1].split(':')[0 : 3]
                            if mth0 != '01' or day0 != '01' or hr0 != '00' or min0 != '00' or sec0 != '00':
                                changes += sim2 + 'UNITS: Change time units date to "' + yr0 + '-01-01 00:00:00"\n'
                            yr0 = int(yr0); yr1 = yr0 + len(time) - 1 # check if filename is consistent
                            if yr0f != yr0:
                                fileparts[yr0idx] = str(yr0) # change start year
                            if yr1f != yr1:
                                fileparts[yr1idx] = str(yr1) + '.nc4' # change end year
                            cyr0 = climyears[cm][0]; cyr1 = climyears[cm][1] # check year range
                            if yr0 > cyr0:
                                errors += sim2 + 'YEARS: Start year > ' + str(cyr0) + '\n'
                            if yr1 < cyr1:
                                errors += sim2 + 'YEARS: End year < ' + str(cyr1) + '\n'
                        except:
                            changes += sim2 + 'UNITS: Change time units to "growing seasons since "<start_year>-01-01 00:00:00"\n'
                
                vfile = vlist[varfullidx] + '_' + fileparts[cropidx]
                if not vfile in fvars:
                    errors += sim2 + 'VARIABLE NAME: Variable name in file is inconsistent with filename\n'
                else:
                    nvars = 0; vfile = '' # try to find variable name
                    for n in fvars:
                        if not n in ['time', 'lat', 'latitude', 'lon', 'longitude']:
                            nvars += 1
                            vfile = n
                    if nvars > 1: vfile = '' # found more than one variable
                
                if vfile != '':
                    v = ncf.variables[vfile]
                    
                    # check dimensions
                    if not lat is None and not 'lat' in v.dimensions and not 'latitude' in v.dimensions:
                        errors += sim2 + 'DIMENSION: Variable must be function of lat\n'
                    if not lon is None and not 'lon' in v.dimensions and not 'longitude' in v.dimensions:
                        errors += sim2 + 'DIMENSION: Variable must be function of lon\n'
                    if not time is None and not 'time' in v.dimensions:
                        errors += sim2 + 'DIMENSION: Variable must be function of time\n'
                    
                    # check fill value
                    if not '_FillValue' in v.ncattrs():
                        changes += sim2 + 'FILL VALUE: Add fill value = 1e20 attribute to variable\n'
                    elif v.getncattr('_FillValue') != float32(1e20):
                        changes += sim2 + 'FILL VALUE: Change fill value to 1e20\n'
                    
                    # check units
                    if not 'units' in v.ncattrs():
                        changes += sim2 + 'UNITS: Add units = "' + vunits[varfullidx] + '" to variable\n'
                    elif v.units != vunits[varfullidx]:
                        changes += sim2 + 'UNITS: Change units to "' + vunits[varfullidx] + '"\n'
                    
                    v = v[:] # convert to numpy array
                    npts = v.size / 4 # divide by 4 to get rough approximation of land points
                    
		    # check if masked array
		    if not isMaskedArray(v):
		        fatal += sim2 + '\n'
		        continue
                    
		    # check percent unmasked
                    punmasked = 100. * (v.size - v.mask.sum()) / npts
                    if punmasked < 10.:
                        errors += sim2 + 'COVERAGE: Spatial coverage is less than 10%\n'
                    
                    # check negative values
                    if vlist[varfullidx] != 'sumt':
                        nneg = (v < 0.).sum()
                        if nneg:
                            pneg = 100. * nneg / npts
                            pneg = '{:.2f}'.format(pneg)
                            errors += sim2 + 'NEGATIVE: ' + pneg + '% of values are negative\n'
                    
                    # check range
                    plower = 100. * (v < vranges[varfullidx][0]).sum() / npts
                    if plower > 0.1:
                        plower = '{:.2f}'.format(plower)
                        errors += sim2 + 'RANGE: ' + plower + '% of values < ' + str(vranges[varfullidx][0]) + '\n'
                    phigher = 100. * logical_and(v > vranges[varfullidx][1], v != 1e20).sum() / npts
                    if phigher > 0.1:
                        phigher = '{:.2f}'.format(phigher)
                        errors += sim2 + 'RANGE: ' + phigher + '% of values > ' + str(vranges[varfullidx][1]) + '\n'
                        
                    # maximum and minimum values
                    vmax = v.max()
		    info_datmat[1] = '{:<9.3f}'.format(vmax)
                    if len(info_datmat[1]) > 9: info_datmat[1] = '{:<9.3e}'.format(vmax)
                    vmin = v.min()
		    info_datmat[2] = '{:<9.3f}'.format(vmin)
                    if len(info_datmat[2]) > 9: info_datmat[2] = '{:<9.3e}'.format(vmin)
                    
                    # number of points within range
                    info_datmat[3] = str(logical_and(v >= vranges[varfullidx][0], v <= vranges[varfullidx][1]).sum() / len(time))
                    
                    # variable name
                    info_datmat[4] = vfile
                    
                ncf.close()
                
                newname = '_'.join(fileparts)
                if newname != f:
                    changes += sim2 + 'FILENAME: Change to ' + newname + '\n'
                
                info += '{:80s}{:10s}{:10s}{:18s}{:15s}{:14s}{:8s}\n'.format(info_datmat[0], info_datmat[1], \
                    info_datmat[2], info_datmat[3], info_datmat[4], info_datmat[5], info_datmat[6])

    # mark missing data
    cm, cp, vr, sc = where(msg_datmat)
    for i in range(len(cm)):
        msg_data += '{:15s}{:15s}{:15s}{:20s}{:15s}\n'.format(d, climmodels[cm[i]], crops[cp[i]], scenarios[sc[i]], vlist[vr[i]])

    sumfilename = options.sumdir + sep + 'report_' + d + '.txt'
    sumfile = open(sumfilename, 'w')
    sumfile.write('!!! Summary produced by ggcmi.QC.py !!!\n')
    sumfile.write('!!! for model ' + d + ' !!!\n\n')

    header = '{:15s}{:15s}{:15s}{:20s}{:8s}'.format('Model', 'Climate', 'Crop', 'Scenario', 'Variable')
    sumfile.write('FATAL ERRORS\n\n')
    sumfile.write(header + '\n')
    sumfile.write('=' * len(header) + '\n')
    if fatal != '':
        sumfile.write(fatal + '\n')
    else:
        sumfile.write('None\n\n')

    header = '{:15s}{:15s}{:15s}{:20s}{:15s}{:s}'.format('Model', 'Climate', 'Crop', 'Scenario', 'Variable', 'Error')
    sumfile.write('DATA QUALITY WARNINGS\n\n')
    sumfile.write(header + '\n')
    sumfile.write('=' * len(header) + '\n')
    if errors != '':
        sumfile.write(errors + '\n')
    else:
        sumfile.write('None\n\n')

    header = '{:15s}{:15s}{:15s}{:20s}{:15s}{:s}'.format('Model', 'Climate', 'Crop', 'Scenario', 'Variable', 'Change')
    sumfile.write('METADATA ISSUES\n\n')
    sumfile.write(header + '\n')
    sumfile.write('=' * len(header) + '\n')
    if changes != '':
        sumfile.write(changes + '\n')
    else:
        sumfile.write('None\n\n')
    
    header = '{:15s}{:15s}{:s}'.format('Model', 'Climate', 'Crop')
    sumfile.write('MISSING SIMULATIONS\n\n')
    sumfile.write(header + '\n')
    sumfile.write('=' * len(header) + '\n')
    if msg_sims != '':
        sumfile.write(msg_sims + '\n')
    else:
        sumfile.write('None\n\n')
        
    header = '{:80s}{:10}{:10}{:18}{:15}{:14}{:8}'.format('File', 'MaxVal', 'MinVal', 'GridcellsInRange', 'VarName', 'YearRange', 'NumYears')
    sumfile.write('BASIC INFO\n\n')
    sumfile.write(header + '\n')
    sumfile.write('=' * len(header) + '\n')
    if info != '':
        sumfile.write(info + '\n')
    else:
        sumfile.write('None\n\n')
    
    header = '{:15s}{:15s}{:15s}{:20s}{:s}'.format('Model', 'Climate', 'Crop', 'Scenario', 'Variable')
    sumfile.write('MISSING FILES\n\n')
    sumfile.write(header + '\n')
    sumfile.write('=' * len(header) + '\n')
    if msg_data != '':
        sumfile.write(msg_data[: -1])
    else:
        sumfile.write('None')

    sumfile.close()
