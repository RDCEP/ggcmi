from numpy import zeros, ones
from netCDF4 import Dataset as nc
from numpy.ma import masked_array, unique

class AggMaskLoader(object):
    def __init__(self, filename, varnames = None):
        f = nc(filename) # open file
        
        if varnames is None: # no variables specified
            varnames = f.variables.keys()
            varnames = [v for v in varnames if not v in ['lat', 'lon']] # remove lat, lon
            varnames += ['global']
        else:
            if not isinstance(varnames, list): # make list
                varnames = [varnames]
        
        self.dat = {'names': [], 'units': [], 'longnames': [], 'data': []}

        for v in varnames:
            if v != 'global':
                var = f.variables[v]
                self.dat['names'].append(v)
                self.dat['units'].append(var.units if 'units' in var.ncattrs() else '')
                self.dat['longnames'].append(var.long_name if 'long_name' in var.ncattrs() else '')
                self.dat['data'].append(var[:])
            else:
                nlats = f.variables['lat'].size
                nlons = f.variables['lon'].size

                self.dat['names'].append('global') # global mask
                self.dat['units'].append('')
                self.dat['longnames'].append('')
                self.dat['data'].append(masked_array(ones((nlats, nlons)), mask = zeros((nlats, nlons))))

        f.close()

    def names(self):     return self.dat['names']
    
    def units(self):     return self.dat['units']
    
    def longnames(self): return self.dat['longnames']
    
    def data(self):      return self.dat['data']
    
    def udata(self):
        ud = []
        for i in range(len(self.dat['names'])):
            udi = unique(self.dat['data'][i])
            udi = udi[~udi.mask] # remove fill value from list
            ud.append(udi)
        return ud