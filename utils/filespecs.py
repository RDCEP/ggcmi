from netCDF4 import Dataset as nc

class FileSpec(object):
    def __init__(self, filename): self.filename = filename

    def append(self, varname, var, dims, units, longname):
        with nc(self.filename, 'a') as f: # append to file
            vvar = f.createVariable(varname, 'f4', dims, zlib = True, shuffle = False, complevel = 9, fill_value = 1e20)
            vvar[:] = var
            vvar.units = units
            vvar.long_name = longname

class BiasCorrectFile(FileSpec):
    def __init__(self, filename, gadm, time, scen, dt, mp, cr):
        super(BiasCorrectFile, self).__init__(filename)

        with nc(filename, 'w', format = 'NETCDF4_CLASSIC') as f:
            f.createDimension('gadm0', len(gadm)) # create gadm
            gadmvar = f.createVariable('gadm0', 'i4', 'gadm0')
            gadmvar[:] = gadm
            gadmvar.units = 'GADM L0 index'
            gadmvar.long_name = '253 countries'

            f.createDimension('time', len(time)) # create time
            timevar = f.createVariable('time', 'i4', 'time')
            timevar[:] = time - time[0]
            timevar.units = 'years since {:d}-01-01'.format(int(time[0]))
            timevar.long_name = 'time'

            f.createDimension('scen', len(scen)) # create scen
            scenvar = f.createVariable('scen', 'i4', 'scen')
            scenvar[:] = range(1, len(scen) + 1)
            scenvar.units = 'mapping'
            scenvar.long_name = ', '.join(scen)

            f.createDimension('dt', len(dt)) # create dt
            dtvar = f.createVariable('dt', 'i4', 'dt')
            dtvar[:] = range(1, len(dt) + 1)
            dtvar.units = 'mapping'
            dtvar.long_name = ', '.join(dt)

            f.createDimension('mp', len(mp)) # create mp
            mpvar = f.createVariable('mp', 'i4', 'mp')
            mpvar[:] = range(1, len(mp) + 1)
            mpvar.units = 'mapping'
            mpvar.long_name = ', '.join(mp)

            f.createDimension('cr', len(cr)) # create cr
            crvar = f.createVariable('cr', 'i4', 'cr')
            crvar[:] = range(1, len(cr) + 1)
            crvar.units = 'mapping'
            crvar.long_name = ', '.join(cr)

class MultimetricsFile(FileSpec):
    def __init__(self, filename, gadm, scen, times, climate, crop, dt, mp, cr):
        super(MultimetricsFile, self).__init__(filename)

        with nc(filename, 'w', format = 'NETCDF4_CLASSIC') as f:
            gadmvar = f.createVariable('gadm0', 'i4', 'gadm0')
            gadmvar[:] = gadm
            gadmvar.units = 'GADM L0 index'
            gadmvar.long_name = '253 countries'

            f.createDimension('scen', len(scen))
            scenvar = f.createVariable('scen', 'i4', 'scen')
            scenvar[:] = range(1, len(scen) + 1)
            scenvar.units = 'mapping'
            scenvar.long_name = ', '.join(scen)

            f.createDimension('time_range', len(times))
            timesvar = f.createVariable('time_range', 'i4', 'time_range')
            timesvar[:] = range(1, len(times) + 1)
            timesvar.units = 'mapping'
            timesvar.long_name = ', '.join(times)

            f.createDimension('climate', len(climate))
            climatevar = f.createVariable('climate', 'i4', 'climate')
            climatevar[:] = range(1, len(climate) + 1)
            climatevar.units = 'mapping'
            climatevar.long_name = ', '.join(climate)

            f.createDimension('crop', len(crop))
            cropvar = f.createVariable('crop', 'i4', 'crop')
            cropvar[:] = range(1, len(crop) + 1)
            cropvar.units = 'mapping'
            cropvar.long_name = ', '.join(crop)

            f.createDimension('dt', len(dt))
            dtvar = f.createVariable('dt', 'i4', 'dt')
            dtvar[:] = range(1, len(dt) + 1)
            dtvar.units = 'mapping'
            dtvar.long_name = ', '.join(dt)

            f.createDimension('mp', len(mp))
            mpvar = f.createVariable('mp', 'i4', 'mp')
            mpvar[:] = range(1, len(mp) + 1)
            mpvar.units = 'mapping'
            mpvar.long_name = ', '.join(mp)

            f.createDimension('cr', len(cr))
            crvar = f.createVariable('cr', 'i4', 'cr')
            crvar[:] = range(1, len(cr) + 1)
            crvar.units = 'mapping'
            crvar.long_name = ', '.join(cr)