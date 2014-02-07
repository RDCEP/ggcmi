#!/usr/bin/env python

# import modules
import itertools, re, sys, numpy.ma as ma
from os import listdir
from netCDF4 import Dataset as nc
from optparse import OptionParser
from os.path import sep, isfile, isdir
from numpy import ones, where, isnan, cos, pi, resize, ceil, double, reshape, repeat

def tmap(X, M, g):
    """ X is X(g, t), M is M(lat, lon), g is g(g) """
    nlats, nlons = M.shape
    ng, nt = X.shape
    ret = ones((nt, nlats, nlons)) # default to one
    for i in range(ng):
        latidx, lonidx = where(M == g[i])
        for j in range(nt):
            if not X.mask[i, j]:
                ret[j, latidx, lonidx] = X[i, j]
    return ret
def filterfiles(listing):
    files = []
    for l in listing:
        if isfile(l):
            files.append(l)
    return files
def long2short(cropin):
    long_list = ['maize', 'wheat', 'soy', 'rice', 'sorghum', 'millet', \
                 'managed_grass', 'sugarcane', 'barley', 'rapeseed', \
                 'rye', 'sugar_beet']
    short_list = ['mai', 'whe', 'soy', 'ric', 'sor', 'mil', \
                  'mgr', 'sug', 'bar', 'rap', \
                  'rye', 'sgb']
    cropout = short_list[long_list.index(cropin)] if cropin in long_list else ''
    return cropout
def createnc(filename, time, lat, lon, scen, irr, tunits):
    with nc(filename, 'w', format = 'NETCDF4_CLASSIC') as f:
        f.createDimension('time', len(time)) # time
        timevar = f.createVariable('time', 'i4', ('time',))
        timevar[:] = time
        timevar.units = tunits
        timevar.long_name = 'time'
        f.createDimension('lat', len(lat)) # latitude
        latvar = f.createVariable('lat', 'f4', ('lat',))
        latvar[:] = lat
        latvar.units = 'degrees_north'
        latvar.long_name = 'latitude'
        f.createDimension('lon', len(lon)) # longitude
        lonvar = f.createVariable('lon', 'f4', ('lon',))
        lonvar[:] = lon
        lonvar.units = 'degrees_east'
        lonvar.long_name = 'longitude'    
        f.createDimension('scen', len(scen)) # scenario
        scenvar = f.createVariable('scen', 'i4', ('scen',))
        scenvar[:] = range(1, len(scen) + 1)
        scenvar.units = 'mapping'
        scenvar.long_name = ', '.join(scen)
        f.createDimension('irr', len(irr)) # irr
        irrvar = f.createVariable('irr', 'i4', ('irr',))
        irrvar[:] = range(1, len(irr) + 1)
        irrvar.units = 'mapping'
        irrvar.long_name = ', '.join(irr)

parser = OptionParser()
parser.add_option("-b", "--batch", dest = "batch", default = "1", type = "int",
                  help = "Batch to process")
parser.add_option("-n", "--numbatches", dest = "num_batches", default = "64", type = "int",
                  help = "Total number of batches")
parser.add_option("-d", "--dir", dest = "dir", default = "", type = "string",
                  help = "Directory in which to perform downscaling")
parser.add_option("-m", "--mod", dest = "mod", default = "pDSSAT,pAPSIM", type = "string",
                  help = "Comma-separated list of crop models (* = all models)")
parser.add_option("-w", "--weath", dest = "weath", default = "AgCFSR,AgMERRA,CFSR,ERAI,GRASP,Princeton,WATCH,WFDEI.CRU,WFDEI.GPCC", type = "string",
                  help = "Comma-separated list of weather datasets (* = all weather datasets)")
parser.add_option("-c", "--crop", dest = "crop", default = "maize", type = "string",
                  help = "Comma-separated list of crops (* = crops, excluding 'others')")
parser.add_option("-i", "--landuseir", dest = "landuseir", default = "", type = "string",
                  help = "Landuse (weight) mask file for irrigation", metavar = "FILE")
parser.add_option("-r", "--landuserf", dest = "landuserf", default = "", type = "string",
                  help = "Landuse (weight) mask file for rainfed", metavar = "FILE")
parser.add_option("-v", "--variables", dest = "variables", default = "", type = "string",
                  help = "Comma-separated list of variables")
parser.add_option("-o", "--outdir", dest = "outdir", default = "", type = "string",
                  help = "Output directory to save results")
parser.add_option("--corrdir", dest = "corrdir", default = "", type = "string",
                  help = "Directory of bias-corrected files")
parser.add_option("--mask", dest = "mask", default = "", type = "string",
                  help = "Mask file", metavar = "FILE")
parser.add_option("--corrmethod", dest = "corrmethod", default = "", type = "string",
                  help = "Correction method to apply")
options, args = parser.parse_args()

rootdir = options.dir # root directory
if rootdir[-1] == sep: rootdir = rootdir[: -1] # remove final separator

if options.mod == '*': # model
    models = listdir(rootdir) # process all models
    if 'upload_stats' in models: models.remove('upload_stats')
    if 'aggregations' in models: models.remove('aggregations')
else:
    models = options.mod.split(',') # model, climate, and crop names

corrfiles = [f for f in listdir(options.corrdir) if f.endswith('.nc4')] # bias-corrected files

totfiles = [] # total files to create
for i in models:
    if options.weath == '*':
        climates = listdir(rootdir + sep + i)
    else:
        climates = options.weath.split(',')
    for j in climates:
        if options.crop == '*':
            crops = listdir(sep.join([rootdir, i, j]))
	    if 'others' in crops: crops.remove('others')
        else:
            crops = options.crop.split(',')
        for k in crops:
            d = sep.join([i, j, k])
            dname = rootdir + sep + d
            if isdir(dname): # if directory
                listing = [dname + sep + l for l in listdir(dname)]
                files = filterfiles(listing)
                if len(files): # if contains files
                    freg = '_'.join([i.lower(), j.lower(), 'hist', long2short(k)])
                    bcf = [f for f in corrfiles if re.search(freg, f)]
                    if len(bcf): # found bias-corrected file
                        totfiles.append([d, bcf[0]])
nfiles = len(totfiles)

batch = options.batch # find out start and end indices for batch
numbatches = options.num_batches
bz = int(ceil(double(nfiles) / numbatches))
si = bz * (batch - 1)
ei = nfiles if batch == numbatches else min(si + bz, nfiles)

if si >= nfiles: # no work for processor to do
    print 'No jobs for processor to perform. Exiting . . .'
    sys.exit()

totfiles = totfiles[si : ei] # select files for batch
nfiles = len(totfiles)

uqcrops = list(set([f[0].split(sep)[2].capitalize() for f in totfiles])) # find unique crops (capitalize first letters)
ncrops = len(uqcrops)

landmasksir = [0] * ncrops # load weight masks
landuseir = nc(options.landuseir) # land use IR mask
for i in range(ncrops):
    irvars = landuseir.variables.keys()
    varidx = [uqcrops[i] in v for v in irvars].index(True)
    landmasksir[i] = landuseir.variables[irvars[varidx]][:]
landuseir.close()
landmasksrf = [0] * ncrops
landuserf = nc(options.landuserf) # land use RF mask
for i in range(ncrops):
    rfvars = landuserf.variables.keys()
    varidx = [uqcrops[i] in v for v in rfvars].index(True)
    landmasksrf[i] = landuserf.variables[rfvars[varidx]][:]
landuserf.close()

with nc(options.mask) as f: # load aggregation mask
    gadm_map = f.variables['gadm0'][:]
    lat = f.variables['lat'][:]
    lon = f.variables['lon'][:]

vars = options.variables.split(',') # get variables

irr = ['noirr', 'firr', 'sum'] # irr values
rf_idx = irr.index('noirr')
ir_idx = irr.index('firr')
sum_idx = irr.index('sum')
    
nlats, nlons = gadm_map.shape # get dimensions
nirr = len(irr)

area = 100 * (111.2 / 2) ** 2 * cos(pi * lat / 360) # get area(lat, lon)
area = resize(area, (nlons, nlats)).T

print 'START INDEX = ' + str(si) + ', END IDX = ' + str(ei) + ' . . .'
for i in range(nfiles): # iterate over subdirectories
    dir = totfiles[i][0]
    bcfile = totfiles[i][1]
    print 'PROCESSING INDEX = ' + str(si + i) + ': DIRECTORY ' + str(rootdir + sep + dir) + ' . . .'

    cidx = uqcrops.index(dir.split(sep)[2].capitalize()) # crop
    cropname = long2short(uqcrops[cidx].lower())

    dname = rootdir + sep + dir # get files
    files = listdir(dname)
    fileslist = filterfiles([dname + sep + l for l in files])
    if not len(fileslist):
        print 'No files found. Skipping directory . . .'
        continue

    with nc(options.corrdir + sep + bcfile) as f: # load bias-corrected file
        gadm = f.variables['gadm0_index'][:]
        time = f.variables['time'][:]
        tunits = f.variables['time'].units
        scen = f.variables['scen'].long_name.split(', ')
        detrend_methods = f.variables['detrend_methods'].long_name.split(', ')
        corr_methods = f.variables['corr_methods'].long_name.split(', ')
        yield_sim = f.variables['yield_sim'][:]
        yield_retrend = f.variables['yield_retrend'][:]

    nodt_idx = detrend_methods.index('none') # select "none" detrending
    yield_sim = yield_sim[:, :, :, nodt_idx]

    corr_idx = corr_methods.index(options.corrmethod) # select correction method
    yield_retrend = yield_retrend[:, :, :, corr_idx]

    time += int(re.findall(r'\d+', tunits)[0]) # get time in years

    yield_sim = ma.masked_where(isnan(yield_sim), yield_sim) # convert NaNs to masked
    yield_retrend = ma.masked_where(isnan(yield_retrend), yield_retrend)

    nt = len(time); nscen = len(scen) # get dimensions

    filename = options.outdir + sep + bcfile.split('.')[0] + '.rescaled.nc4' # create netCDF4 file
    createnc(filename, time - time[0], lat, lon, scen, irr, tunits)

    areair = resize(area * landmasksir[cidx], (nt, nlats, nlons)) # get areas
    areair = reshape(repeat(areair, nscen), (nt, nlats, nlons, nscen))
    arearf = resize(area * landmasksrf[cidx], (nt, nlats, nlons))
    arearf = reshape(repeat(arearf, nscen), (nt, nlats, nlons, nscen))

    areatot = areair + arearf # get total area
    areatot = ma.masked_where(areatot == 0, areatot)
    
    var = ma.zeros((nt, nlats, nlons, nscen, nirr)) # masked array
    for v in vars: # iterate over variables
        var[:] = 0. # reset
        var.mask = ones(var.shape) # mask all
        
        for s in itertools.product(scen, irr[: 2]): # iterate over scenarios
            sidx = scen.index(s[0])
            iidx = irr.index(s[1])
            scenname = s[0].split('_')
            if len(scenname) == 1:
                scen_irr = scenname[0] + '_' + s[1]
            else:
                scen_irr = '_'.join([scenname[0], s[1], scenname[1]])

            varfile = [f for f in files if re.search(scen_irr + '_' + v, f)] # load file
            if not len(varfile): continue
            varfile = varfile[0]
            print 'Processing', varfile, '. . .'
            
            with nc(sep.join([dname, varfile])) as f:
                fvars = f.variables.keys() # load data
                vidx = [v in s for s in fvars].index(True)
                fvar = f.variables[fvars[vidx]]

                funits = fvar.units if 'units' in fvar.ncattrs() else '' # load units and long name
                flongname = fvar.long_name if 'long_name' in fvar.ncattrs() else ''

                fvar = fvar[:] # change NaNs to zeros
                fvar[isnan(fvar)] = 0.

                ftime = f.variables['time'][:] # load time
                ftime_units = f.variables['time'].units

            ftime += int(re.findall(r'\d+', ftime_units)[0]) - 1 # get time in years

            for t in range(nt): # cut to times in bias-corrected file
                tidx = where(ftime == time[t])[0][0]
                var[t, :, :, sidx, iidx] = fvar[tidx, :, :]

            yr = yield_retrend[:, :, sidx] # yields for rescaling
            ys = yield_sim[:, :, sidx]

            num = tmap(yr, gadm_map, gadm) # scale
            denom = tmap(ys, gadm_map, gadm)
            var[:, :, :, sidx, iidx] *= num / denom
        
        # get sum
        var[:, :, :, :, sum_idx] = (arearf * var[:, :, :, :, rf_idx] + \
                                    areair * var[:, :, :, :, ir_idx]) / \
                                    areatot

        with nc(filename, 'a') as f: # append variables
            vvar = f.createVariable(v + '_' + cropname, 'f4', ('time', 'lat', 'lon', 'scen', 'irr',), zlib = True, complevel = 9, fill_value = 1e20)
            vvar[:] = var
            vvar.units = funits
            vvar.long_name = flongname