#!/usr/bin/env python

# import modules
from os import listdir
from os.path import sep, exists
from optparse import OptionParser
from netCDF4 import Dataset as nc
from numpy import inf, float32, zeros, ones, where, diff, array_equal, logical_and

def climate_years(clim_names):
    yrs = [0] * len(clim_names)
    for i in range(len(clim_names)):
        if clim_names[i] == 'AgCFSR':
            yrs[i] = [1980, 2010]
        elif clim_names[i] == 'AgMERRA':
            yrs[i] = [1980, 2010]
        elif clim_names[i] == 'CFSR':
            yrs[i] = [1980, 2010]
        elif clim_names[i] == 'ERAI':
            yrs[i] = [1979, 2010]
        elif clim_names[i] == 'GRASP':
            yrs[i] = [1961, 2010]
        elif clim_names[i] == 'Princeton':
            yrs[i] = [1948, 2008]
        elif clim_names[i] == 'WATCH':
            yrs[i] = [1958, 2001]
        elif clim_names[i] == 'WFDEI.CRU':
            yrs[i] = [1979, 2012]
        elif clim_names[i] == 'WFDEI.GPCC':
            yrs[i] = [1979, 2009]
        else:
            raise Exception('Unknown climate')
    return yrs

def abr_crop_names(crop_names):
    abr_names = [0] * len(crop_names)
    for i in range(len(crop_names)):
        if crop_names[i] == 'maize':
            abr_names[i] = 'mai'
        elif crop_names[i] == 'wheat':
            abr_names[i] = 'whe'
        elif crop_names[i] == 'soy':
            abr_names[i] = 'soy'
        elif crop_names[i] == 'rice':
            abr_names[i] = 'ric'
        elif crop_names[i] == 'sorghum':
            abr_names[i] = 'sor'
        elif crop_names[i] == 'millet':
            abr_names[i] = 'mil'
        elif crop_names[i] == 'managed_grass':
            abr_names[i] = 'mgr'
        else:
            raise Exception('Unknown crop')
    return abr_names

# parse inputs
parser = OptionParser()
parser.add_option("-d", "--dir", dest = "dir", default = "AgMIP.output", type = "string",
                  help = "Directory in which to perform verification")
parser.add_option("-m", "--mod", dest = "mod", default = "pDSSAT.pm,pDSSAT.pt,pAPSIM", type = "string",
                  help = "Comma-separated list of crop models to verify")
parser.add_option("-w", "--weath", dest = "weath", default = "AgCFSR,AgMERRA,CFSR,ERAI,GRASP,Princeton,WATCH,WFDEI.CRU,WFDEI.GPCC", type = "string",
                  help = "Comma-separated list of weather datasets to verify")
parser.add_option("-c", "--crop", dest = "crop", default = "maize", type = "string",
                  help = "Comma-separated list of crops to verify")
parser.add_option("-s", "--sum", dest = "sum", default = "summary.txt", type = "string",
                  help = "Name of summary file", metavar = "FILE")
options, args = parser.parse_args()

dir = options.dir
if dir[-1] == sep: dir = dir[: -1] # remove final separator

# crop models
dirs = options.mod.split(',')

# climate data
climmodels = options.weath.split(',')
climyears = climate_years(climmodels)
numclims = len(climmodels)

# crops
crops = options.crop.split(',')
cropabbr = abr_crop_names(crops)
numcrops = len(crops)

# scenarios
dssat_pm_scens = ['fullharm_noirr', 'fullharm_firr', 'default_noirr', 'default_firr', 'harmnon_noirr', 'harmnon_firr']
dssat_pt_scens = ['fullharm_noirr', 'fullharm_firr']
apsim_scens = dssat_pm_scens

# variables
vlist = ['yield', 'pirrww', \
         'plant-day', 'maty-day', 'aet', 'initr', 'gsprcp', \
         'anth-day', 'biom', 'leach', 'gsrsds', 'sumt', 'sco2', 'sn2o']
vranges = [[0, 100], [0, 10000], \
           [1, 366], [1, 366], [0, 10000], [0, 1000], [0, 10000], \
           [1, 366], [0, 100], [0, 1000], [0, 10000], [-inf, inf], [-inf, inf], [-inf, inf]]
vunits = ['t ha-1 yr-1', 'mm yr-1', \
          'day of year', 'days from planting', 'mm yr-1', 'kg ha-1 yr-1', 'mm ha-1 yr-1', \
          'days from planting', 't ha-1 yr-1', 'kg ha-1 yr-1', 'w m-2 yr-1', 'deg C-days yr-1', 'kg C ha-1', 'kg N2O-N ha-1']
numtotvars = len(vlist)
mandatoryvars = zeros((numtotvars,), dtype = bool)
mandatoryvars[: 2] = [True, True]

msg_sims = '' # missing simulations
msg_data = '' # missing data within a simulation
changes = ''  # changes to make
errors = ''   # errors in data

msg_datmat_dssat_pm = ones((numclims, numcrops, numtotvars, len(dssat_pm_scens)), dtype = bool)
msg_datmat_dssat_pt = ones((numclims, numcrops, numtotvars, len(dssat_pt_scens)), dtype = bool)
msg_datmat_apsim = ones((numclims, numcrops, numtotvars, len(apsim_scens)), dtype = bool)

# iterate through all permutations
for d in dirs:
    if d == 'pDSSAT.pm':
        scenarios = dssat_pm_scens
    elif d == 'pDSSAT.pt':
        scenarios = dssat_pt_scens
    elif d == 'pAPSIM':
        scenarios = apsim_scens
    else:
        continue # skip directory

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
                if len(fileparts) != 10:
                    continue
                
                # variable
                varname = fileparts[5].lower()
                if varname in vlist:
                    varidx = vlist.index(varname)
                else:
                    continue # cannot find variable
                
                # scenario
                scen = fileparts[3].lower() + '_' + fileparts[4].lower() # combine
                if scen in scenarios:
                    scenidx = scenarios.index(scen)
                else:
                    continue # cannot find scenario
                scen = scen.split('_') # resplit
                
                # mark as found
                if d == 'pDSSAT.pm':
                    msg_datmat_dssat_pm[cm, cp, varidx, scenidx] = False
                elif d == 'pDSSAT.pt':
                    msg_datmat_dssat_pt[cm, cp, varidx, scenidx] = False
                elif d == 'pAPSIM':
                    msg_datmat_apsim[cm, cp, varidx, scenidx] = False
                sim2 = sim + '{:20s}{:15s}'.format(scenarios[scenidx], vlist[varidx])
                
                # change filename as necessary
                if fileparts[0] != d.lower():
                    fileparts[0] = d.lower() # change model name
                if fileparts[1] != climmodels[cm].lower():
                    fileparts[0] = climmodels[cm].lower() # change climate name
                if fileparts[2] != 'hist':
                    fileparts[2] = 'hist' # change to hist
                if fileparts[3] != scen[0]:
                    fileparts[3] = scen[0] # change first part of scenario
                if fileparts[4] != scen[1]:
                    fileparts[4] = scen[1] # change second part of scenario
                if fileparts[5] != varname:
                    fileparts[5] = varname # change var to lowercase
                if fileparts[6] != cropabbr[cp]:
                    fileparts[6] = cropabbr[cp] # change crop name
                if fileparts[7] != 'annual':
                    fileparts[7] = 'annual' # change to annual
                
                yr0 = int(fileparts[8]); yr1 = int(fileparts[9].split('.')[0])
                cyr0 = climyears[cm][0]; cyr1 = climyears[cm][1]
                if yr0 != cyr0:
                    fileparts[8] = str(cyr0) # change start year
                if yr1 != cyr1:
                    fileparts[9] = str(cyr1) + '.nc4' # change end year
                
                newname = '_'.join(fileparts)
                if newname != f:
                    changes += sim2 + 'FILENAME: Change to ' + newname + '\n'
                
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
                            if yr0 != str(cyr0) or mth0 != '01' or day0 != '01' or hr0 != '00' or min0 != '00' or sec0 != '00':
                                changes += sim2 + 'UNITS: Change time units date to "' + str(cyr0) + '-01-01 00:00:00"\n'
                        except:
                            changes += sim2 + 'UNITS: Change time units to "growing seasons since ' + str(cyr0) + '-01-01 00:00:00"\n'
                
                vfile = vlist[varidx] + '_' + cropabbr[cp]
                if not vfile in fvars:
                    errors += sim2 + 'VARIABLE NAME: Variable name in file is inconsistent with filename\n'
                else:
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
                        changes += sim2 + 'UNITS: Add units = "' + vunits[varidx] + '" to variable\n'
                    elif v.units != vunits[varidx]:
                        changes += sim2 + 'UNITS: Change units to "' + vunits[varidx] + '"\n'
                    
                    v = v[:] # convert to numpy array
                    npts = v.size
                    
                    # check percent unmasked
                    punmasked = 100. * (npts - v.mask.sum()) / npts
                    if punmasked < 10.:
                        errors += sim2 + 'COVERAGE: Spatial coverage is less than 10%\n'
                    
                    # check negative values
                    if vlist[varidx] != 'sumt':
                        nneg = (v < 0.).sum()
                        if nneg:
                            pneg = 100. * nneg / npts
                            pneg = '{:.2f}'.format(pneg)
                            errors += sim2 + 'NEGATIVE: ' + pneg + '% of values are negative\n'
                    
                    # check range
                    plower = 100. * (v < vranges[varidx][0]).sum() / npts
                    if plower > 0.1:
                        plower = '{:.2f}'.format(plower)
                        errors += sim2 + 'RANGE: ' + plower + '% of values < ' + str(vranges[varidx][0]) + '\n'
                    phigher = 100. * logical_and(v > vranges[varidx][1], v != 1e20).sum() / npts
                    if phigher > 0.1:
                        phigher = '{:.2f}'.format(phigher)
                        errors += sim2 + 'RANGE: ' + phigher + '% of values > ' + str(vranges[varidx][1]) + '\n'
                    
                ncf.close()

# mark missing data
if 'pDSSAT.pm' in dirs:
    cm, cp, vr, sc = where(msg_datmat_dssat_pm)
    for i in range(len(cm)):
        msg_data += '{:15s}{:15s}{:15s}{:20s}{:15s}\n'.format('pDSSAT.pm', climmodels[cm[i]], crops[cp[i]], dssat_pm_scens[sc[i]], vlist[vr[i]])
if 'pDSSAT.pt' in dirs:
    cm, cp, vr, sc = where(msg_datmat_dssat_pt)
    for i in range(len(cm)):
        msg_data += '{:15s}{:15s}{:15s}{:20s}{:15s}\n'.format('pDSSAT.pt', climmodels[cm[i]], crops[cp[i]], dssat_pt_scens[sc[i]], vlist[vr[i]])
if 'pAPSIM' in dirs:
    cm, cp, vr, sc = where(msg_datmat_apsim)
    for i in range(len(cm)):
        msg_data += '{:15s}{:15s}{:15s}{:20s}{:15s}\n'.format('pAPSIM', climmodels[cm[i]], crops[cp[i]], apsim_scens[sc[i]], vlist[vr[i]])

sumfile = open(options.sum, 'w')
sumfile.write('!!! Summary produced by GGCMIverifier.py !!!\n\n')

header = '{:15s}{:15s}{:s}'.format('Model', 'Climate', 'Crop')
sumfile.write('MISSING SIMULATIONS\n\n')
sumfile.write(header + '\n')
sumfile.write('=' * len(header) + '\n')
if msg_sims != '':
    sumfile.write(msg_sims + '\n')
else:
    sumfile.write('None\n\n')

header = '{:15s}{:15s}{:15s}{:20s}{:s}'.format('Model', 'Climate', 'Crop', 'Scenario', 'Variable')
sumfile.write('MISSING DATA\n\n')
sumfile.write(header + '\n')
sumfile.write('=' * len(header) + '\n')
if msg_data != '':
    sumfile.write(msg_data + '\n')
else:
    sumfile.write('None\n\n')

header = '{:15s}{:15s}{:15s}{:20s}{:15s}{:s}'.format('Model', 'Climate', 'Crop', 'Scenario', 'Variable', 'Change')
sumfile.write('REQUIRED CHANGES\n\n')
sumfile.write(header + '\n')
sumfile.write('=' * len(header) + '\n')
if changes != '':
    sumfile.write(changes + '\n')
else:
    sumfile.write('None\n\n')

header = '{:15s}{:15s}{:15s}{:20s}{:15s}{:s}'.format('Model', 'Climate', 'Crop', 'Scenario', 'Variable', 'Error')
sumfile.write('ERRORS\n\n')
sumfile.write(header + '\n')
sumfile.write('=' * len(header) + '\n')
if errors != '':
    sumfile.write(errors[: -1])
else:
    sumfile.write('None')

sumfile.close()