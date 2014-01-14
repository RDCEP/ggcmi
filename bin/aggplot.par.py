#!/usr/bin/env python

# import modules
import csv, sys
from os import sep, listdir
from os.path import basename
import matplotlib.pyplot as plt
from optparse import OptionParser
from netCDF4 import Dataset as nc
from numpy import nan, zeros, ceil, double

class Plotter(object):
    tickfs = 10 # fontsizes
    titlefs = 18
    colors = ['b', 'g', 'r', 'k', 'y']
    def __init__(self, filename, outdir, fmt):
        self.filename = filename
        self.outdir = outdir
        self.fmt = fmt
        with nc(filename) as f:
            self.varkeys = f.variables.keys() # get available variables
            self.scennames = [s.strip(' ') for s in f.variables['scen'].long_name.split(',')] # get scenario names
            self.irr = f.variables['irr'][:] # get irrigation values and names
            self.irrnames = [i.strip(' ') for i in f.variables['irr'].long_name.split(',')]
            self.time = f.variables['time'][:] # get time and units
            self.tunits = f.variables['time'].units
    def splot(self, dat, var, varunits, labels, title, filename, lines = None, ncol = 3):
        # first dimension of dat is time, second dimension is ensemble
        for i in range(len(labels)): # time series
            if lines is None:
                plt.plot(self.time, dat[:, i], lw = 2, label = labels[i])
            else:
                plt.plot(self.time, dat[:, i], lines[i], lw = 2, label = labels[i])
        plt.ylim(plt.ylim()[0], plt.ylim()[1] * 4. / 3) # extend y upper limit by third
        plt.legend(ncol = ncol, prop = {'size': self.tickfs})
        plt.grid()
        plt.xticks(fontsize = self.tickfs)
        plt.yticks(fontsize = self.tickfs)
        plt.xlim(self.time[0], self.time[-1])
        plt.xlabel(self.tunits, fontsize = self.titlefs)
        plt.ylabel(var + ' (' + varunits + ')', fontsize = self.titlefs)
        plt.title(title, fontsize = self.titlefs)
        plt.tight_layout()
        plt.savefig(filename)
        plt.close()
    def bplot(self, dat, var, varunits, labels, title, filename, xlabel = 'Region', rotation = 0):
        plt.boxplot(dat) # box-and-whisker
        plt.grid()
        plt.xticks(range(1, len(labels) + 1), labels, fontsize = self.tickfs, rotation = rotation)
        plt.yticks(fontsize = self.tickfs)
        plt.xlabel(xlabel, fontsize = self.titlefs)
        plt.ylabel(var + ' (' + varunits + ')', fontsize = self.titlefs)
        plt.title(title, fontsize = self.titlefs)
        plt.tight_layout()
        plt.savefig(filename)
        plt.close()

class MaskPlotter(Plotter):
    def __init__(self, filename, meta, mask, outdir, fmt):
        super(MaskPlotter, self).__init__(filename, outdir, fmt)
        with nc(self.filename) as f:
            self.area = f.variables['area_' + mask][:]
            self.areaglob = f.variables['area_global'][:]
            self.indices = f.variables[mask + '_index'][:]
        metadata = [] # load meta data
        with open(meta, 'rU') as f:
            for row in csv.reader(f, delimiter = '\t'):
                metadata.append(row)
        self.metadic = {}
        for i in range(len(metadata)):
            md = metadata[i][0].split(',')
            self.metadic[md[0]] = md[1]
        filesplit = basename(filename).split('_') # get model, climate, crop
        self.model = filesplit[0]
        self.climate = filesplit[1]
        self.crop = filesplit[3]
        self.yr0 = filesplit[5]
	self.yr1 = filesplit[6].split('.')[0]
        totarea = self.area[:, 0, 0, :].sum(axis = 1)
        self.sidx = [a[0] for a in sorted(enumerate(totarea), key = lambda x: x[1], reverse = True)]
        self.mask = mask
    	self.prefix = '_'.join([self.model, self.climate, self.yr0, self.yr1]) + '_'
    def plot(self, var, n):
        with nc(self.filename) as f:
            dat = f.variables[var + '_' + self.mask][:]
            glob = f.variables[var + '_global'][:]
            varunits = f.variables[var + '_' + self.mask].units
        dat[dat.mask] = nan # replace fill values with nan
        
        nt = len(self.time)
        ns = len(self.scennames)
        nir = len(self.irrnames)
        pdat = zeros((nt, ns, nir + 1, n + 1))

        for i in range(len(self.scennames)): # get data
            scen = self.scennames[i]
            for j in range(len(self.irr)):
                ir = self.irrnames[j]
                for k in range(n):
                    si = self.sidx[k]
                    pdat[:, i, j, k] = dat[si, :, i, j]
                    pdat[:, i, nir, k] += self.area[si, :, i, j] * dat[si, :, i, j]
                pdat[:, i, j, n] = glob[:, i, j] # add global separately
                pdat[:, i, nir, n] += self.areaglob[:, i, j] * glob[:, i, j]
            for j in range(n):
                si = self.sidx[j]
                pdat[:, i, nir, j] /= self.area[si, :, i, :].sum(axis = 1)
            pdat[:, i, nir, n] /= self.areaglob[:, i, :].sum(axis = 1)

        regionlabels = [0] * (n + 1) # region labels
        for i in range(n):
            si = self.sidx[i]
            regionlabels[i] = self.metadic[str(self.indices[si])]
        regionlabels[n] = 'global'
        
        seriesscenlabels = [0] * (ns * nir)
        boxscenlabels = [0] * (ns * (nir + 1))
        scenlines = [0] * (ns * nir)
        cnt1 = 0; cnt2 = 0
        for i in range(len(self.scennames)): # scenario labels and colors
            scen = self.scennames[i]
            col = self.colors[i]
            for j in range(len(self.irr)):
                ir = self.irrnames[j]
                seriesscenlabels[cnt1] = scen + ' / ' + ir
                boxscenlabels[cnt2] = scen + '\n' + ir
                scenlines[cnt1] = col + '-' if not j else col + '--'
                cnt1 += 1; cnt2 += 1
            boxscenlabels[cnt2] = scen + '\nsum'
            cnt2 += 1

        for i in range(len(self.scennames)): # time-region plots
            scen = self.scennames[i]
            for j in range(len(self.irr)):
                ir = self.irrnames[j]
                title = '-'.join([self.model, self.climate, self.crop]) + '\n' + ' / '.join([scen, ir])
                filehead = sep.join([self.outdir, self.crop, var, self.mask, 'reg', self.prefix])
		filename1 = filehead + '_'.join(['tseries', scen, ir]) + '.' + self.fmt
		filename2 = filehead + '_'.join(['bplot', scen, ir]) + '.' + self.fmt
                self.splot(pdat[:, i, j, :], var, varunits, regionlabels, title, filename1, ncol = 4) 
                self.bplot(pdat[:, i, j, :], var, varunits, regionlabels, title, filename2, rotation = 30)        

        seriesdat = zeros((nt, ns * nir, n + 1)) # reorder data
        boxdat = zeros((nt, ns * (nir + 1), n + 1))
        cnt1 = 0; cnt2 = 0
        for i in range(len(self.scennames)):
            for j in range(len(self.irr)):
                seriesdat[:, cnt1, :] = pdat[:, i, j, :]
                boxdat[:, cnt2, :] = pdat[:, i, j, :]
                cnt1 += 1; cnt2 += 1
            boxdat[:, cnt2, :] = pdat[:, i, nir, :]
            cnt2 += 1
        
        for i in range(n + 1): # time-scenario plots
            si = self.indices[self.sidx[i]]
            title = '-'.join([self.model, self.climate, self.crop]) + '\n' + regionlabels[i]
            app = str(si) if i != n else 'global'
	    filehead = sep.join([self.outdir, self.crop, var, self.mask, 'scen', self.prefix])
            filename1 = filehead + '_'.join(['tseries', 'r' + app]) + '.' + self.fmt
            filename2 = filehead + '_'.join(['bplot', 'r' + app]) + '.' + self.fmt
            self.splot(seriesdat[:, :, i], var, varunits, seriesscenlabels, title, filename1, scenlines, ncol = ns)    
            self.bplot(boxdat[:, :, i], var, varunits, boxscenlabels, title, filename2, 'Scenario', rotation = 30) 

# parse inputs
parser = OptionParser()
parser.add_option("-b", "--batch", dest = "batch", default = "1", type = "int",
                  help = "Batch to process")
parser.add_option("-n", "--numbatches", dest = "num_batches", default = "64", type = "int",
                  help = "Total number of batches")
parser.add_option("-d", "--inputdir", dest = "inputdir", default = "", type = "string",
                  help = "Directory of aggregated files to be plotted")
parser.add_option("-v", "--vars", dest = "vars", default = "yield", type = "string",
                  help = "Comma-separated list of variables to plot")
parser.add_option("-a", "--aggmasks", dest = "aggmasks", default = "fpu", type = "string",
                  help = "Comma-separated list of aggregation masks to plot (e.g., fpu, basin, region)")
parser.add_option("-m", "--meta", dest = "meta", default = "", type = "string",
                  help = "Comma-separated list of meta files corresponding to masks")
parser.add_option("-p", "--nplots", dest = "nplots", default = 5, type = "int",
                  help = "Number of plots to make (created in order of descending land area)")
parser.add_option("-f", "--fmt", dest = "fmt", default = "png", type = "string",
                  help = "File format in which to save plot(s) (e.g., png, pdf)")
parser.add_option("-o", "--outdir", dest = "outdir", default = "", type = "string",
                  help = "Output directory to save plots")
options, args = parser.parse_args()

vars = options.vars.split(',')      # variable names
masks = options.aggmasks.split(',') # masks
meta = options.meta.split(',')      # meta files
nplots = options.nplots             # number of plots

files = listdir(options.inputdir) # files to plot
nfiles = len(files)

batch = options.batch # find out start and end indices for batch
numbatches = options.num_batches
bz = int(ceil(double(nfiles) / numbatches))
si = bz * (batch - 1)
ei = nfiles if batch == numbatches else min(si + bz, nfiles)

if si >= nfiles: # no work for processor to do
    print 'No jobs for processor to perform. Exiting . . .'
    sys.exit()

for f in files[si : ei]: # iterate over files
    print 'Plotting file', f, '. . .'
    for m in range(len(masks)): # iterate over masks
        mp = MaskPlotter(options.inputdir + sep + f, meta[m], masks[m], options.outdir, options.fmt)
        for v in vars: # iterate over variables
            mp.plot(v, nplots)
