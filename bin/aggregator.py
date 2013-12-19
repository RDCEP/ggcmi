#!/usr/bin/env python

# import modules
import itertools, sys, re, time as tm
from os import listdir
from os.path import sep, exists
from optparse import OptionParser
from netCDF4 import Dataset as nc
from numpy.ma import masked_array, unique, masked_where
from numpy import pi, zeros, ones, cos, resize, where, ceil, double

# HELPER FUNCTION
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
    f.createDimension('irr', 2) # irr
    irrvar = f.createVariable('irr', 'i4', ('irr',))
    irrvar[:] = range(1, 3)
    irrvar.units = 'mapping'
    irrvar.long_name = 'ir, rf'
    for i in range(len(rnames)): # index variables
        rname = rnames[i] + '_index'
        f.createDimension(rname, len(rdata[i]))
        rvar = f.createVariable(rname, 'i4', (rname,))
        rvar[:] = rdata[i]
        rvar.units = runits[i]
        rvar.long_name = rlongnames[i]
    f.close() # close file

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
        nlats = self.lat.size
        nlons = self.lon.size
        self.dat['names'].append('global') # add global mask as a trick
        self.dat['units'].append('')
        self.dat['longnames'].append('')
        self.dat['data'].append(masked_array(zeros((nlats, nlons)), mask = zeros((nlats, nlons))))
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
parser.add_option("-m", "--mod", dest = "mod", default = "pDSSAT.pm,pDSSAT.pt,pAPSIM", type = "string",
                  help = "Comma-separated list of crop models")
parser.add_option("-w", "--weath", dest = "weath", default = "AgCFSR,AgMERRA,CFSR,ERAI,GRASP,Princeton,WATCH,WFDEI.CRU,WFDEI.GPCC", type = "string",
                  help = "Comma-separated list of weather datasets")
parser.add_option("-c", "--crop", dest = "crop", default = "maize", type = "string",
                  help = "Comma-separated list of crops")
parser.add_option("-i", "--landuseir", dest = "landuseir", default = "", type = "string",
                  help = "Landuse (weight) mask file for irrigation", metavar = "FILE")
parser.add_option("-r", "--landuserf", dest = "landuserf", default = "", type = "string",
                  help = "Landuse (weight) mask file for rainfed", metavar = "FILE")
parser.add_option("-a", "--agg", dest = "agg", default = "", type = "string",
                  help = "Comma-separated list of aggregation mask files")
parser.add_option("-o", "--outdir", dest = "outdir", default = "", type = "string",
                  help = "Comma-separated list of output directories to save results for different aggregation masks")
options, args = parser.parse_args()

# CONSTANTS
# =========
fillv = 1e20
yieldthr1 = 0.1 # t/ha
yieldthr2 = 50

dssat_pm_scens = ['fullharm', 'default', 'harmnon'] # scenarios
dssat_pt_scens = ['fullharm']
apsim_scens = dssat_pm_scens
dssat_pm_scens_full = ['fullharm_noirr', 'fullharm_firr', 'default_noirr', 'default_firr', 'harmnon_noirr', 'harmnon_firr']
dssat_pt_scens_full = ['fullharm_noirr', 'fullharm_firr']
apsim_scens_full = dssat_pm_scens_full

dssatvars = ['yield', 'pirrww', 'plant-day', 'maty-day', 'aet', 'gsprcp', 'anth-day', 'biom', 'gsrsds', 'sumt'] # variables
apsimvars = ['yield', 'pirrww', 'plant-day', 'maty-day', 'aet', 'gsprcp', 'anth-day', 'biom', 'gsrsds', 'sumt', 'initr', 'leach', 'sco2', 'sn2o']
# =========

rootdir = options.dir # root directory
if rootdir[-1] == sep: rootdir = rootdir[: -1] # remove final separator

models = options.mod.split(',') # model, climate, and crop names
climates = options.weath.split(',')
crops = options.crop.split(',')
aggmasks = options.agg.split(',') # aggregration masks and output directories
aggmaskdirs = options.outdir.split(',')

totfiles = [] # total files to create
for i in itertools.product(models, climates, crops):
    d = sep.join(i)
    if exists(rootdir + sep + d) and len(listdir(rootdir + sep + d)):
        for j in range(len(aggmasks)):
            totfiles.append([d, aggmasks[j], aggmaskdirs[j]])
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

uqcrops = list(set([f[0].split(sep)[2].title() for f in totfiles])) # find unique crops (capitalize first letters)
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

uqaggmasks = list(set([f[1] for f in totfiles])) # find unique aggregration masks
naggmasks = len(uqaggmasks)

aggmaskobjs = [] # load aggregration masks
for i in range(naggmasks): aggmaskobjs.append(AggMask(uqaggmasks[i]))

nlats = len(aggmaskobjs[0].lat) # number of latitudes and longitudes
nlons = len(aggmaskobjs[0].lon)

lat = aggmaskobjs[0].lat # compute area as function of latitude
area = 100 * (111.2 / 2)**2 * cos(pi * lat / 360)
area = resize(area, (nlons, nlats)).T

print 'START INDEX = ' + str(si) + ', END IDX = ' + str(ei) + ' . . .'
for i in range(nfiles): # iterate over subdirectories
    dir = totfiles[i][0]
    aggmask = totfiles[i][1]
    outdir = totfiles[i][2]
    print 'PROCESSING INDEX = ' + str(si + i) + ': DIRECTORY ' + str(rootdir + sep + dir) + ' . . .'
    
    cidx = uqcrops.index(dir.split(sep)[2].title()) # crop
    
    aidx = uqaggmasks.index(aggmask) # aggregration mask
    amask = aggmaskobjs[aidx]
    anames = amask.names()
    aunits = amask.units()
    alongnames = amask.longnames()
    adata = amask.data()
    audata = amask.udata()
    nmasks = len(anames)

    mod = dir.split(sep)[0] # variables and scenarios
    if mod == 'pDSSAT.pm':
        vars = dssatvars
        scens = dssat_pm_scens
        scens_full = dssat_pm_scens_full
    elif mod == 'pDSSAT.pt':
        vars = dssatvars
        scens = dssat_pt_scens
        scens_full = dssat_pt_scens_full
    elif mod == 'pAPSIM':
        vars = apsimvars
        scens = apsim_scens
        scens_full = apsim_scens_full

    files = listdir(rootdir + sep + dir) # get files

    f = sep.join([rootdir, dir, files[0]]) # use first file to get time and units
    ncf = nc(f)
    time = ncf.variables['time'][:]
    tunits = ncf.variables['time'].units
    ncf.close()

    nv = len(vars)  # number of variables
    nt = len(time)  # number of times
    ns = len(scens) # number of scenarios    

    filename = files[0].split('_') # use first file to get filename
    filename.remove(filename[3])
    filename.remove(filename[3])
    filename.remove(filename[3])
    filename = outdir + sep + '_'.join(filename) # save in output directory

    createnc(filename, time, tunits, scens, audata[: -1], anames[: -1], aunits[: -1], alongnames[: -1]) # create nc file

    averages = [0] * nmasks; areas = [0] * nmasks # final averages and areas
    adataresize = [0] * nmasks
    for j in range(nmasks):
        sz = audata[j].size
        adataresize[j] = resize(adata[j], (nt, nlats, nlons)) # resize maps
        averages[j] = fillv * ones((nv, sz, nt, ns, 2)) # preallocate
        areas[j] = zeros((sz, nt, ns, 2))

    vunits = [''] * nv
    for j in range(len(scens_full)): # iterate over scenarios
        scen_irr = scens_full[j]
        
        scen_irr_split = scen_irr.split('_') # find scenario and irr
        sidx = scens.index(scen_irr_split[0])
        iidx = int(scen_irr_split[1] != 'firr')

        weight = landmasksir[cidx] if not iidx else landmasksrf[cidx] # weights

        yieldfile = [f for f in files if re.search(scen_irr + '_yield', f)][0] # pull yield mask
        yf = nc(sep.join([rootdir, dir, yieldfile]))
        yvars = yf.variables.keys()
        yidx = ['yield' in v for v in yvars].index(True)
        yieldvar = yf.variables[yvars[yidx]][:]
        yieldvar = yieldvar.transpose((2, 0, 1))
        yieldvar = masked_where(yieldvar < yieldthr1, yieldvar) # mask yields based on thresholds
        yieldvar = masked_where(yieldvar > yieldthr2, yieldvar)
        yieldmask = yieldvar.mask
        yf.close()

        aselect = [0] * nmasks; vartmp = [0] * nmasks
        for k in range(nmasks):
            sz = audata[k].size
            acommon = masked_where(yieldmask, adataresize[k]) # mask
            vartmp[k] = zeros((sz, nlats, nlons)) # for storing temporary variable
            aselect[k] = zeros((sz, nt, nlats, nlons), dtype = bool)
            for m in range(sz):
                aselect[k][m] = acommon == audata[k][m]
                areas[k][m, :, sidx, iidx] = (weight * area * aselect[k][m]).sum(axis = 2).sum(axis = 1)

        for k in range(nv): # iterate over variables
            t0 = tm.time()
            
            varfile = [f for f in files if re.search(scen_irr + '_' + vars[k], f)]
            if not len(varfile): continue
            varfile = varfile[0]
            print 'Processing', varfile, '. . .'

            vf = nc(sep.join([rootdir, dir, varfile])) # load file

            fvars = vf.variables.keys() # get variable and units
            vidx = [vars[k] in v for v in fvars].index(True)
            var = vf.variables[fvars[vidx]]
            if not j: vunits[k] = var.units if 'units' in var.ncattrs() else ''
            var = var[:].transpose((2, 0, 1))
            vf.close()
                       
            print '  Processing aggregration masks . . .'
            for m in range(nmasks):
                print '    ' + anames[m] + ' . . .'
                for t in range(nt):
                    ridx = areas[m][:, t, sidx, iidx] > 0.
                    if ridx.sum():
                        vartmp[m][:] = 0.
                        idx1, idx2, idx3 = where(aselect[m][:, t, :, :])
                        vartmp[m][idx1, idx2, idx3] = var[t, idx2, idx3] * weight[idx2, idx3] * area[idx2, idx3] * aselect[m][idx1, t, idx2, idx3]        
                        vsum = vartmp[m].sum(axis = 2).sum(axis = 1)
                        averages[m][k, ridx, t, sidx, iidx] = vsum[ridx] / areas[m][ridx, t, sidx, iidx]

            print '  Time to process =', tm.time() - t0, 'seconds . . .'

    f = nc(filename, 'a', format = 'NETCDF4_CLASSIC') # append variables
    for j in range(nmasks):
        name = anames[j]
        dims = ('time', 'scen', 'irr')
        if name != 'global': dims = (name + '_index',) + dims
        areav = f.createVariable('area_' + name, 'f4', dims, zlib = True, complevel = 9)
        areav[:] = areas[j].squeeze()
        areav.units = 'hectares'
        areav.long_name = name + ' harvested area'
        for k in range(nv):
            avev = f.createVariable(vars[k] + '_' + name, 'f4', dims, fill_value = fillv, zlib = True, complevel = 9)
            avev[:] = averages[j][k].squeeze()
            avev.units = vunits[k]
            avev.long_name = 'average ' + name + ' ' + vars[k]
    f.close()

print 'DONE!'