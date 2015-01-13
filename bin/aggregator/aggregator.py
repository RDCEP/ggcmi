#!/usr/bin/env python

# import modules
import sys, random, time as tm
from os import listdir
from datetime import datetime
from optparse import OptionParser
from netCDF4 import Dataset as nc
from numpy.ma import masked_array, unique, masked_where
from os.path import sep, getmtime, isfile, isdir, basename
from numpy import pi, zeros, ones, cos, resize, mod, where, ceil, double, isnan, logical_not, logical_and

# HELPER FUNCTIONS
def createnc(filename, time, tunits, scens, rdata, rnames, runits, rlongnames):
    f = nc(filename, 'w', format = 'NETCDF4_CLASSIC') # create file
    f.createDimension('time', len(time)) # time
    timevar = f.createVariable('time', 'i4', ('time',))
    timevar[:] = time
    timevar.units = tunits
    timevar.long_name = 'time'
    f.createDimension('scen', len(scens)) # scenario
    scenvar = f.createVariable('scen', 'i4', ('scen',))
    scenvar[:] = range(1, len(scens) + 1)
    scenvar.units = 'mapping'
    scenvar.long_name = ', '.join(scens)
    f.createDimension('irr', 3) # irr
    irrvar = f.createVariable('irr', 'i4', ('irr',))
    irrvar[:] = range(1, 4)
    irrvar.units = 'mapping'
    irrvar.long_name = 'ir, rf, sum'
    for i in range(len(rnames)): # index variables
        f.createDimension(rnames[i], len(rdata[i]))
        rvar = f.createVariable(rnames[i], 'i4', rnames[i])
        rvar[:] = rdata[i]
        rvar.units = runits[i]
        rvar.long_name = rlongnames[i]
    f.close() # close file
def filterfiles(listing):
    files = []
    for l in listing:
        if isfile(l) and not l.endswith('.old'): # skip old files
            files.append(l)
    return files
def findfile(files, scen_irr, var):
    for f in files:
        if '_%s_%s_' % (scen_irr, var) in f:
            return f
    return []
def getcropabbr(crop):
    cmap = {'maize': 'mai', 'wheat': 'whe', 'soy': 'soy', 'rice': 'ric', \
            'sorghum': 'sor', 'millet': 'mil', 'managed_grass': 'mgr'}
    if crop in cmap:
        return cmap[crop]
    else:
        raise Exception('Unrecognized crop name')
def getgrowing(dir, crop, ir):
    cmap = {'maize': 'Maize', 'wheat': 'Wheat', 'soy': 'Soybeans', 'rice': 'Rice', \
            'sorghum': 'Sorghum', 'millet': 'Millet'}
    irmap = {'firr': 'ir', 'noirr': 'rf'}
    if crop in cmap and ir in irmap:
        filestart = cmap[crop] + '_' + irmap[ir]
        flist = listdir(dir)
        file = [f for f in flist if f.startswith(filestart) and f.endswith('nc4')]
        if file:
            return dir + sep + file[0] # first return
        else:
            raise Exception('Cannot map to growing season file: %s, %s, %s' % (dir, crop, ir))
    else:
        raise Exception('Cannot map to growing season file: %s, %s, %s' % (dir, crop, ir))
def shiftdata(data, pdate, hdate, yearthr = 1):
    latidx, lonidx = where(logical_and(pdate >= hdate, hdate >= yearthr))
    datac = data.copy()
    shiftd = datac[:, latidx, lonidx]
    shiftd[1 :] = shiftd[: -1] # shift data down one
    shiftd[0].mask = True
    datac[:, latidx, lonidx] = shiftd
    return datac

class AggMask(object):
    def __init__(self, filename):
        f = nc(filename)
        varnames = f.variables.keys()
        varnames = [v for v in varnames if not v in ['lat', 'lon']] # remove lat, lon
        self.dat = {'names': [], 'units': [], 'longnames': [], 'data': []}
        for v in varnames:
            var = f.variables[v]
            self.dat['names'].append(v)
            self.dat['units'].append(var.units if 'units' in var.ncattrs() else '')
            self.dat['longnames'].append(var.long_name if 'long_name' in var.ncattrs() else '')
            self.dat['data'].append(var[:])
        self.lat = f.variables['lat'][:] # load latitude and longitude
        self.lon = f.variables['lon'][:]        
        f.close()
    def names(self):     return self.dat['names']
    def units(self):     return self.dat['units']
    def longnames(self): return self.dat['longnames']
    def data(self):      return self.dat['data']
    def udata(self): # unique data values
        ud = []
        for i in range(len(self.dat['names'])):
            udi = unique(self.dat['data'][i])
            udi = udi[~udi.mask] # remove fill value from list
            ud.append(udi)
        return ud

# parse inputs
parser = OptionParser()
parser.add_option("-b", "--batch", dest = "batch", default = "1", type = "int",
                  help = "Batch to process")
parser.add_option("-n", "--numbatches", dest = "num_batches", default = "64", type = "int",
                  help = "Total number of batches")
parser.add_option("-d", "--dir", dest = "dir", default = "", type = "string",
                  help = "Directory in which to perform aggregration")
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
parser.add_option("-a", "--agg", dest = "agg", default = "", type = "string",
                  help = "Comma-separated list of aggregation mask files")
parser.add_option("-t", "--tsfile", dest = "tsfile", default = "timestamps.txt", type = "string",
                  help = "File of timestamps", metavar = "FILE")
parser.add_option("-g", "--growdir", dest = "growdir", default = "", type = "string",
                  help = "Directory of growing season files")
parser.add_option("-o", "--outdir", dest = "outdir", default = "", type = "string",
                  help = "Comma-separated list of output directories to save results for different aggregation masks")
options, args = parser.parse_args()

print 'RUNNING AGGREGATOR'
print '=================='

fillv, yieldthr1, yieldthr2, yearthr = 1e20, 0.1, 30, 1 # some constants

rootdir = options.dir # root directory
if rootdir[-1] == sep: rootdir = rootdir[: -1] # remove final separator

if options.mod == '*': # model
    models = listdir(rootdir) # process all models
    if 'upload_stats' in models: models.remove('upload_stats')
    if 'aggregations' in models: models.remove('aggregations')
else:
    models = options.mod.split(',') # model, climate, and crop names

aggmasks = options.agg.split(',') # aggregration masks and output directories
aggmaskdirs = options.outdir.split(',')

totfiles = [] # total files to create
for i in models:
    if options.weath == '*':
        climates = listdir(rootdir + sep + i)
        if 'output.not.corrected' in climates: climates.remove('output.not.corrected')
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
                    for m in range(len(aggmasks)):
                        totfiles.append([d, aggmasks[m], aggmaskdirs[m]])
nfiles = len(totfiles)

random.seed(0) # shuffle directories to distribute large jobs
random.shuffle(totfiles)

batch = options.batch # find out start and end indices for batch
numbatches = options.num_batches
bz = int(ceil(double(nfiles) / numbatches))
si = bz * (batch - 1)
ei = nfiles if batch == numbatches else min(si + bz, nfiles)

if si >= nfiles: # no work for processor to do
    print 'No jobs for processor to perform. Exiting . . .'
    sys.exit()

totfiles = totfiles[si : ei] # select files for batch
nfiles   = len(totfiles)

uqcrops = list(set([f[0].split(sep)[2].capitalize() for f in totfiles])) # find unique crops (capitalize first letters)
ncrops  = len(uqcrops)

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

uqaggmasks = list(set([f[1] for f in totfiles])) # find unique aggregration masks
naggmasks  = len(uqaggmasks)

aggmaskobjs = [] # load aggregration masks
for i in range(naggmasks): aggmaskobjs.append(AggMask(uqaggmasks[i]))

nlats, nlons = len(aggmaskobjs[0].lat), len(aggmaskobjs[0].lon) # number of latitudes and longitudes

lat  = aggmaskobjs[0].lat # compute area as function of latitude
area = 100 * (111.2 / 2) ** 2 * cos(pi * lat / 180)
area = resize(area, (nlons, nlats)).T

tslines = open(options.tsfile).readlines() # load timestamps file
tsdic = {}
for i in range(len(tslines)):
    tsparts = tslines[i].strip('\n').split(' ')
    if len(tsparts) == 3: # contains directory, date, time
        tsdic[tsparts[0]] = ' '.join(tsparts[1 : 3])

tsfile = open('timestamps_' + str(batch) + '.txt', 'w') # create new timestamp file

print 'START INDEX = ' + str(si) + ', END IDX = ' + str(ei) + ' . . .'
for i in range(nfiles): # iterate over subdirectories
    dir = totfiles[i][0]
    aggmask = totfiles[i][1]
    outdir = totfiles[i][2]
    print 'PROCESSING INDEX = ' + str(si + i) + ': DIRECTORY ' + str(rootdir + sep + dir) + ' . . .'

    ts = datetime.fromtimestamp(getmtime(rootdir + sep + dir)).strftime('%Y-%m-%d %H:%M')
    tsfile.write(dir + ' ' + ts + '\n')
    if dir in tsdic: # check if file has been modified
        if ts == tsdic[dir]: # same timestamp means already processed
            print 'Timestamp has not changed. Skipping directory . . .'
            continue

    cropname = dir.split(sep)[2] # crop
    cropabbr = getcropabbr(cropname)
    cidx = uqcrops.index(cropname.capitalize())

    aidx = uqaggmasks.index(aggmask) # aggregration mask
    amask = aggmaskobjs[aidx]
    anames = amask.names()
    aunits = amask.units()
    alongnames = amask.longnames()
    adata = amask.data()
    audata = amask.udata()
    nmasks = len(anames)

    dname = rootdir + sep + dir
    files = listdir(dname)
    fileslist = filterfiles([dname + sep + l for l in files]) # get files
    if not len(fileslist):
        print 'No files found. Skipping directory . . .'
        continue
    files = [basename(l) for l in fileslist]
    files = [f for f in files if '_' + cropabbr + '_' in f]

    vars = []; scens = []; scens_full = [] # get variables and scenarios
    for j in range(len(files)):
        fs = files[j].split('_')
        if 'noirr' in fs:
            fs.remove('noirr')
            ir = 'noirr'
        else:
            fs.remove('firr')
            ir = 'firr'
    	if len(fs) == 9:
            vars.append(fs[4])
            scens.append(fs[3])
            scens_full.append(fs[3] + '_' + ir)
        elif len(fs) == 10: # pt files, etc.
            vars.append(fs[5])
            scens.append(fs[3] + '_' + fs[4])
 	    scens_full.append('_'.join([fs[3], ir, fs[4]]))
	else:
	    continue # skip irregular file
    vars = list(set(vars)); vars.sort()
    scens = list(set(scens)); scens.sort()
    scens_full = list(set(scens_full)); scens_full.sort()

    with nc(sep.join([rootdir, dir, files[0]])) as ncf: # use first file to get time and units
        time = ncf.variables['time'][:]
        tunits = ncf.variables['time'].units
        if tunits == 'year as %Y.%f': # patch for GEPIC
            tunits = 'growing seasons since %d-01-01 00:00:00' % int(time[0])
            time = time - time[0] + 1

    nv, nt, ns = len(vars), len(time), len(scens) # number of variables, times, scenarios

    filename = files[0].split('_') # use first file to get filename
    if len(filename) == 11: filename.remove(filename[3]) # remove the extra label for pt
    filename.remove(filename[3])
    filename.remove(filename[3])
    filename.remove(filename[3])
    filename = outdir + sep + '_'.join(filename) # save in output directory
    createnc(filename, time, tunits, scens, audata, anames, aunits, alongnames) # create nc file

    print 'Preallocating . . .'
    t0 = tm.time()
    averages = [0] * nmasks; areas = [0] * nmasks # final averages and areas
    aselect  = [0] * nmasks; vartmp = [0] * nmasks
    for j in range(nmasks):
        sz = audata[j].size
        averages[j] = masked_array(zeros((nv, sz, nt, ns, 3))) # preallocate
        averages[j].mask = ones(averages[j].shape) # all points are masked
        areas[j] = zeros((sz, nt, ns, 3))
        vartmp[j] = zeros((sz, nlats, nlons)) # for storing temporary variable
        aselect[j] = zeros((sz, nlats, nlons), dtype = bool)
        for k in range(sz):
            aselect[j][k] = adata[j] == audata[j][k]
    print '  Elapsed time =', tm.time() - t0, 'seconds . . .'

    vunits = [''] * nv
    for j in range(len(scens_full)): # iterate over scenarios
        scen_irr = scens_full[j]

        scen_irr_split = scen_irr.split('_') # find scenario and irr
        if len(scen_irr_split) == 2:
            sidx = scens.index(scen_irr_split[0])
        else:
            sidx = scens.index('_'.join(scen_irr_split[0 :: 2]))
	iidx = int(scen_irr_split[1] != 'firr')

        weight = landmasksir[cidx] if not iidx else landmasksrf[cidx] # weights

        print 'Loading planting and maturity dates . . .'
        t0 = tm.time()
        plantingfile = findfile(files, scen_irr, 'plant-day')
        if plantingfile:
            with nc(sep.join([rootdir, dir, plantingfile])) as pf:
               pdate = pf.variables['plant-day_' + cropabbr][0] # first time
        else:
            plantingfile = getgrowing(options.growdir, cropname, scen_irr_split[1])
            with nc(plantingfile) as pf:
                lats = pf.variables['lat'][:]
                lsi = [idxl[0] for idxl in sorted(enumerate(lats), reverse = True, key = lambda x : x[1])]
                pdate = pf.variables['planting day'][lsi] # flip latitudes
        harvestfile = findfile(files, scen_irr, 'maty-day')
        if harvestfile:
            with nc(sep.join([rootdir, dir, harvestfile])) as hf:
                hdate = hf.variables['maty-day_' + cropabbr][0] # first time
        else:
            harvestfile = getgrowing(options.growdir, cropname, scen_irr_split[1])
            with nc(harvestfile) as hf:
                lats = hf.variables['lat'][:]
                lsi = [idxl[0] for idxl in sorted(enumerate(lats), reverse = True, key = lambda x : x[1])]
                hdate = hf.variables['growing season length'][lsi] # flip latitudes
        hdate = mod(pdate + hdate, 366) # convert to Julian day
        hdate[hdate == 0] = 1
        print '  Elapsed time =', tm.time() - t0, 'seconds  . . .'

        print 'Loading yield mask . . .'
        t0 = tm.time()
        yieldfile = findfile(files, scen_irr, 'yield') # pull yield mask
        if yieldfile == []:
            continue
        else:
            yf = nc(sep.join([rootdir, dir, yieldfile]))
	yieldvar = yf.variables['yield_' + cropabbr][:]
        yieldvar = shiftdata(yieldvar, pdate, hdate, yearthr) # shift data
        yieldvar = masked_where(yieldvar < yieldthr1, yieldvar) # mask yields based on thresholds
        yieldvar = masked_where(yieldvar > yieldthr2, yieldvar)
        yieldmask = logical_not(yieldvar.mask)
        yf.close()
        print '  Elapsed time =', tm.time() - t0, 'seconds . . .'

        print 'Computing areas . . .'
        t0 = tm.time()
        for k in range(nmasks):
            for m in range(audata[k].size):
	        areas[k][m, :, sidx, iidx] = (weight * area * aselect[k][m] * yieldmask).sum(axis = 2).sum(axis = 1)
	        areas[k] = masked_where(areas[k] == 0, areas[k])
        print '  Elapsed time =', tm.time() - t0, 'seconds . . .'

        for k in range(nv): # iterate over variables
            t0 = tm.time()

            varfile = findfile(files, scen_irr, vars[k])
            if not len(varfile): continue
            print 'Processing', varfile, '. . .'

            with nc(sep.join([rootdir, dir, varfile])) as vf: # load file
                var = vf.variables[vars[k] + '_' + cropabbr]
                if not j: vunits[k] = var.units if 'units' in var.ncattrs() else ''
                var = var[:]
            var[isnan(var)] = 0. # change NaN values to zero
            var = shiftdata(var, pdate, hdate, yearthr) # shift data

            print '  Processing aggregration masks . . .'
            for m in range(nmasks):
                print '    ' + anames[m] + ' . . .'
                ridx, latidx, lonidx = where(aselect[m])
                vartmp[m][:] = 0. # reset
		for t in range(nt):
                    vartmp[m][ridx, latidx, lonidx] = var[t, latidx, lonidx] * weight[latidx, lonidx] * area[latidx, lonidx] * yieldmask[t, latidx, lonidx] * aselect[m][ridx, latidx, lonidx]
                    vsum = vartmp[m].sum(axis = 2).sum(axis = 1)
                    averages[m][k, :, t, sidx, iidx] = vsum / areas[m][:, t, sidx, iidx]
            print '  Elapsed time =', tm.time() - t0, 'seconds . . .'

    print 'Writing file . . .'
    with nc(filename, 'a', format = 'NETCDF4_CLASSIC') as f: # append variables
        for j in range(nmasks):
            name = anames[j]
            dims = (name, 'time', 'scen', 'irr')
            area1, area2 = areas[j][:, :, :, 0], areas[j][:, :, :, 1]
	    area1[logical_and(area1.mask, ~area2.mask)] = 0.
            area2[logical_and(area2.mask, ~area1.mask)] = 0. 
            areas[j][:, :, :, 2] = area1 + area2 # sum area
            areav = f.createVariable('area_' + name, 'f4', dims, zlib = True, complevel = 9)
            areav[:] = areas[j].squeeze()
            areav.units = 'hectares'
            areav.long_name = name + ' harvested area'
            for k in range(nv):
                av1 = areas[j][:, :, :, 0] * averages[j][k, :, :, :, 0]
		av2 = areas[j][:, :, :, 1] * averages[j][k, :, :, :, 1]
		av1[logical_and(av1.mask, ~av2.mask)] = 0.
		av2[logical_and(av2.mask, ~av1.mask)] = 0.
                averages[j][k, :, :, :, 2] = (av1 + av2) / areas[j][:, :, :, 2] # area average
                avev = f.createVariable(vars[k] + '_' + name, 'f4', dims, fill_value = fillv, zlib = True, complevel = 9)
                avev[:] = averages[j][k].squeeze()            
                avev.units = vunits[k]
                avev.long_name = 'average ' + name + ' ' + vars[k]

tsfile.close() # close timestamps file

print 'DONE!'