#!/usr/bin/env python

# import modules
import csv, re
import numpy.ma as ma
from os import listdir, sep
from optparse import OptionParser
from netCDF4 import Dataset as nc
from numpy import zeros, ones, arange, polyfit, diff

# utility functions
def winave(d, nstart, nend):
    # perform n-point window average of data d
    nd = len(d)
    dave = d.copy()
    numramp = (nend - nstart) / 2
    n = nstart - 2
    for i in range(nd):
        if i <= numramp:
            n += 2
        elif i >= nd - numramp:
            n -= 2
        idx0 = max(i - n / 2, 0)
        idx1 = min(i + n / 2 + 1, nd)
        dave[i] = d[idx0 : idx1].mean()
    dave = dave[n / 2 : -(n / 2)] # remove left and right edges
    return dave

class FAOStatFile(object):
    def __init__(self, csvfile):
        self.data = []
        with open(csvfile) as f: # load data
            for row in csv.reader(f):
                self.data.append(row)
        self.gadmidx = self.data[0].index('gadm code')
        self.yr0idx = self.data[0].index('element') + 1
        self.years = [int(y) for y in self.data[0][self.yr0idx :]]
        self.ndata = len(self.data) - 1

        self.gadm = [] # gadm
        for i in range(self.ndata):
            g = self.data[i + 1][self.gadmidx]
            if g != 'NA': self.gadm.append(int(g))
        self.gadm = sorted(self.gadm) # sort
        
        m = re.match(r'.* \((.*)\)', self.data[1][self.yr0idx - 1]) # units
        self.units = m.group(1) if m else ''

    def getData(self):
        data = zeros((len(self.gadm), len(self.years)))
        for i in range(self.ndata):
            g = self.data[i + 1][self.gadmidx]
            if g != 'NA':
                gidx = self.gadm.index(int(g))
                yld = self.data[i + 1][self.yr0idx :]
                data[gidx] = [0 if y == '' else float(y) for y in yld]
        data = ma.masked_where(data == 0, data) # mask 0 or blank values
        return data

# parse inputs
parser = OptionParser()
parser.add_option("-i", "--input", dest = "input", default = "", type = "string",
                  help = "Input directory of csv files of yields and areas")
options, args = parser.parse_args()

files = listdir(options.input)
filedic = {}
for f in files:
    if f.endswith('.csv'):
        crop = f.split('.')[1]
        if f.startswith('yield.'):
            if not crop in filedic: filedic[crop] = {}
            filedic[crop]['yield'] = options.input + sep + f
        elif f.startswith('area.'):
            if not crop in filedic: filedic[crop] = {}
            filedic[crop]['area'] = options.input + sep + f

firstfile = filedic.values()[0].values()[0] # pull years and gadm from first file
faotest = FAOStatFile(firstfile)
years = faotest.years
gadm = faotest.gadm

filename = 'faostat.{:d}-{:d}.nc4'.format(years[0], years[-1]) # filename

detrend = ['none', 'linear', 'quadratic', 'moving average', 'fractional first difference', 'trend-removed fractional first difference']

f = nc(filename, 'w', format = 'NETCDF4_CLASSIC') # create file
f.createDimension('gadm0_index', len(gadm)) # create gadm index
gadmvar = f.createVariable('gadm0_index', 'i4', ('gadm0_index',), zlib = True, shuffle = False, complevel = 9)
gadmvar[:] = gadm
gadmvar.units = 'GADM L0 index'
gadmvar.long_name = '253 countries'
f.createDimension('time', len(years)) # create time
timevar = f.createVariable('time', 'i4', ('time',), zlib = True, shuffle = False, complevel = 9)
timevar[:] = [y - years[0] for y in years]
timevar.units = 'years since {:d}-01-01'.format(years[0])
timevar.long_name = 'time'
f.createDimension('detrend', len(detrend)) # create detrend
detrendvar = f.createVariable('detrend', 'i4', ('detrend',), zlib = True, shuffle = False, complevel = 9)
detrendvar[:] = range(1, len(detrend) + 1)
detrendvar.units = 'mapping'
detrendvar.long_name = ', '.join(detrend)

for crop, dic in filedic.iteritems():
    if not 'yield' in dic or not 'area' in dic: continue

    yieldfile = dic['yield']
    faoyield = FAOStatFile(yieldfile) # get yield
    yieldarr = faoyield.getData() / 10000. # divide by 10,000 to get tonnes
    if faoyield.years != years:
        print 'Warning: Years not the same in file', yieldfile
    if faoyield.gadm != gadm:
        print 'Warning: gadm not the same in file', yieldfile
    
    areafile = dic['area']
    faoarea = FAOStatFile(dic['area']) # get area
    areaarr = faoarea.getData()
    if faoarea.years != years:
        print 'Warning: Years not the same in file', areafile
    if faoarea.gadm != gadm:
        print 'Warning: gadm not the same in file', areafile
    
    yielddt = ma.zeros((yieldarr.shape + (len(detrend),)))
    yielddt.mask = ones((yielddt.shape))
    for i in range(len(yielddt)): # detrend
        y = yieldarr[i]
        x = arange(0, len(y))
        for j in range(len(detrend)):
            if detrend[j] == 'none':
                yielddt[i, :, j] = y
            elif detrend[j] == 'linear':
                fm = polyfit(x, y, 1)
                yielddt[i, :, j] = y - fm[0] * x - fm[1]
            elif detrend[j] == 'quadratic':
                fm = polyfit(x, y, 2)
                yielddt[i, :, j] = y - fm[0] * x ** 2 - fm[1] * x - fm[2]
            elif detrend[j] == 'moving average':
                yielddt[i, 2 : -2, j] = y[2 : -2] - winave(y, 5, 7)
            elif detrend[j] == 'fractional first difference':
                yielddt[i, 1 :, j] = diff(y) / y[: -1]
            elif detrend[j] == 'trend-removed fractional first difference':
                ffd = diff(y) / y[: -1]
                fm = polyfit(x[1 :], ffd, 1)
                yielddt[i, 1 :, j] = ffd - fm[0] * x[1 :] - fm[1]
            else:
                raise Exception('Unrecognized detrend method')

    yieldvar = f.createVariable('yield' + '_' + crop, 'f4', ('gadm0_index', 'time', 'detrend',), zlib = True, shuffle = False, complevel = 9, fill_value = 1e20) # create yield
    yieldvar[:] = yielddt
    yieldvar.units = 't ha-1 yr-1'
    yieldvar.long_name = 'average gadm0 yield'
    areavar = f.createVariable('area' + '_' + crop, 'f4', ('gadm0_index', 'time',), zlib = True, shuffle = False, complevel = 9, fill_value = 1e20) # create area
    areavar[:] = areaarr
    areavar.units = faoarea.units
    areavar.long_name = 'area harvested'

f.close() # close file