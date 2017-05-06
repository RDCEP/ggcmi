from netCDF4 import Dataset


class FileSpec(object):
    def __init__(self, filename):
        self.filename = filename

    def append(self, varname, var, dims, units, longname):
        with Dataset(self.filename, 'a') as f:
            vvar = f.createVariable(varname, 'f4', dims, zlib=True, shuffle=False, complevel=9, fill_value=1e20)
            vvar[:] = var
            vvar.units = units
            vvar.long_name = longname


class AggregationFile(FileSpec):
    def __init__(self, filename, time, tunits, aggs, aggname, aggunits, agglongname):
        super(AggregationFile, self).__init__(filename)

        with Dataset(filename, 'w', format='NETCDF4_CLASSIC') as f:

            f.createDimension(aggname, len(aggs))
            aggsvar = f.createVariable(aggname, 'i4', aggname)
            aggsvar[:] = aggs
            aggsvar.units = aggunits
            aggsvar.long_name = agglongname

            f.createDimension('time', len(time))
            timevar = f.createVariable('time', 'i4', 'time')
            timevar[:] = time
            timevar.units = tunits
            timevar.long_name = 'time'


class BiasCorrectFile(FileSpec):
    def __init__(self, filename, aggs, aggname, aggunits, agglongname, time, scen, dt, mp, cr):
        super(BiasCorrectFile, self).__init__(filename)

        with Dataset(filename, 'w', format='NETCDF4_CLASSIC') as f:

            f.createDimension(aggname, len(aggs))
            aggsvar = f.createVariable(aggname, 'i4', aggname)
            aggsvar[:] = aggs
            aggsvar.units = aggunits
            aggsvar.long_name = agglongname

            f.createDimension('time', len(time))
            timevar = f.createVariable('time', 'i4', 'time')
            timevar[:] = time - time[0]
            timevar.units = 'years since {:d}-01-01'.format(int(time[0]))
            timevar.long_name = 'time'

            f.createDimension('scen', len(scen))
            scenvar = f.createVariable('scen', 'i4', 'scen')
            scenvar[:] = range(1, len(scen) + 1)
            scenvar.units = 'mapping'
            scenvar.long_name = ', '.join(scen)

            f.createDimension('dt', len(dt))
            dtvar = f.createVariable('dt', 'i4', 'dt')
            dtvar[:] = range(1, len(dt) + 1)
            dtvar.units = 'mapping'
            dtvar.long_name = ', '.join(dt)
            dtvar.note = 'detrend method'

            f.createDimension('mp', len(mp))
            mpvar = f.createVariable('mp', 'i4', 'mp')
            mpvar[:] = range(1, len(mp) + 1)
            mpvar.units = 'mapping'
            mpvar.long_name = ', '.join(mp)
            mpvar.note = 'mean-preserving method'

            f.createDimension('cr', len(cr))
            crvar = f.createVariable('cr', 'i4', 'cr')
            crvar[:] = range(1, len(cr) + 1)
            crvar.units = 'mapping'
            crvar.long_name = ', '.join(cr)
            crvar.note = 'correction method'


class MultimetricsFile(FileSpec):
    def __init__(self, filename, aggs, aggname, aggunits, agglongname, scen, times, dt, mp, cr):
        super(MultimetricsFile, self).__init__(filename)

        with Dataset(filename, 'w', format='NETCDF4_CLASSIC') as f:

            f.createDimension(aggname, len(aggs))
            aggsvar = f.createVariable(aggname, 'i4', aggname)
            aggsvar[:] = aggs
            aggsvar.units = aggunits
            aggsvar.long_name = agglongname

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

            f.createDimension('dt', len(dt))
            dtvar = f.createVariable('dt', 'i4', 'dt')
            dtvar[:] = range(1, len(dt) + 1)
            dtvar.units = 'mapping'
            dtvar.long_name = ', '.join(dt)
            dtvar.note = 'detrend method'

            f.createDimension('mp', len(mp))
            mpvar = f.createVariable('mp', 'i4', 'mp')
            mpvar[:] = range(1, len(mp) + 1)
            mpvar.units = 'mapping'
            mpvar.long_name = ', '.join(mp)
            mpvar.note = 'mean-preserving method'

            f.createDimension('cr', len(cr))
            crvar = f.createVariable('cr', 'i4', 'cr')
            crvar[:] = range(1, len(cr) + 1)
            crvar.units = 'mapping'
            crvar.long_name = ', '.join(cr)
            crvar.note = 'correction method'


class ModelEnsembleFile(FileSpec):
    def __init__(self, filename, metric, aggs, aggname, aggunits, agglongname, time, dt, mp, cr, nm):
        super(ModelEnsembleFile, self).__init__(filename)

        with Dataset(filename, 'w', format='NETCDF4_CLASSIC') as f:

            f.createDimension(aggname, len(aggs))
            aggsvar = f.createVariable(aggname, 'i4', aggname)
            aggsvar[:] = aggs
            aggsvar.units = aggunits
            aggsvar.long_name = agglongname

            f.createDimension('time', len(time))
            timevar = f.createVariable('time', 'i4', 'time')
            timevar[:] = time - time[0]
            timevar.units = 'years since {:d}-01-01'.format(int(time[0]))
            timevar.long_name = 'time'

            f.createDimension('dt', len(dt))
            dtvar = f.createVariable('dt', 'i4', 'dt')
            dtvar[:] = range(1, len(dt) + 1)
            dtvar.units = 'mapping'
            dtvar.long_name = ', '.join(dt)
            dtvar.note = 'detrend method'

            f.createDimension('mp', len(mp))
            mpvar = f.createVariable('mp', 'i4', 'mp')
            mpvar[:] = range(1, len(mp) + 1)
            mpvar.units = 'mapping'
            mpvar.long_name = ', '.join(mp)
            mpvar.note = 'mean-preserving method'

            f.createDimension('cr', len(cr))
            crvar = f.createVariable('cr', 'i4', 'cr')
            crvar[:] = range(1, len(cr) + 1)
            crvar.units = 'mapping'
            crvar.long_name = ', '.join(cr)
            crvar.note = 'correction method'

            f.createDimension('nm', nm)
            nmvar = f.createVariable('nm', 'i4', 'nm')
            nmvar[:] = range(1, nm + 1)
            nmvar.long_name = 'top models used'

            f.createDimension('wt', 2)
            weightedvar = f.createVariable('wt', 'i4', 'wt')
            weightedvar[:] = [1, 2]
            weightedvar.units = 'mapping'
            weightedvar.long_name = 'unweighted, %s-weighted' % metric


class MultimetricsEnsembleFile(FileSpec):
    def __init__(self, filename, aggs, aggname, aggunits, agglongname, times, dt, mp, cr, nm, wt):
        super(MultimetricsEnsembleFile, self).__init__(filename)

        with Dataset(filename, 'w', format='NETCDF4_CLASSIC') as f:

            f.createDimension(aggname, len(aggs))
            aggsvar = f.createVariable(aggname, 'i4', aggname)
            aggsvar[:] = aggs
            aggsvar.units = aggunits
            aggsvar.long_name = agglongname

            f.createDimension('time_range', len(times))
            timesvar = f.createVariable('time_range', 'i4', 'time_range')
            timesvar[:] = range(1, len(times) + 1)
            timesvar.units = 'mapping'
            timesvar.long_name = ', '.join(times)

            f.createDimension('dt', len(dt))
            dtvar = f.createVariable('dt', 'i4', 'dt')
            dtvar[:] = range(1, len(dt) + 1)
            dtvar.units = 'mapping'
            dtvar.long_name = ', '.join(dt)
            dtvar.note = 'detrend method'

            f.createDimension('mp', len(mp))
            mpvar = f.createVariable('mp', 'i4', 'mp')
            mpvar[:] = range(1, len(mp) + 1)
            mpvar.units = 'mapping'
            mpvar.long_name = ', '.join(mp)
            mpvar.note = 'mean-preserving method'

            f.createDimension('cr', len(cr))
            crvar = f.createVariable('cr', 'i4', 'cr')
            crvar[:] = range(1, len(cr) + 1)
            crvar.units = 'mapping'
            crvar.long_name = ', '.join(cr)
            crvar.note = 'correction method'

            f.createDimension('nm', nm)
            nmvar = f.createVariable('nm', 'i4', 'nm')
            nmvar[:] = range(1, nm + 1)
            nmvar.long_name = 'top models used'

            f.createDimension('wt', len(wt))
            weightedvar = f.createVariable('wt', 'i4', 'wt')
            weightedvar[:] = range(1, len(wt) + 1)
            weightedvar.units = 'mapping'
            weightedvar.long_name = ', '.join(wt)


class RescaledFile(FileSpec):
    def __init__(self, filename, time, lat, lon, irr):
        super(RescaledFile, self).__init__(filename)

        with Dataset(filename, 'w', format='NETCDF4_CLASSIC') as f:
            f.createDimension('time', len(time))
            timevar = f.createVariable('time', 'i4', 'time')
            timevar[:] = time - time[0]
            timevar.units = 'years since {:d}-01-01'.format(int(time[0]))
            timevar.long_name = 'time'

            f.createDimension('lat', len(lat))
            latvar = f.createVariable('lat', 'f4', 'lat')
            latvar[:] = lat
            latvar.units = 'degrees_north'
            latvar.long_name = 'latitude'

            f.createDimension('lon', len(lon))
            lonvar = f.createVariable('lon', 'f4', 'lon')
            lonvar[:] = lon
            lonvar.units = 'degrees_east'
            lonvar.long_name = 'longitude'

            f.createDimension('irr', len(irr))
            irrvar = f.createVariable('irr', 'i4', 'irr')
            irrvar[:] = range(1, len(irr) + 1)
            irrvar.units = 'mapping'
            irrvar.long_name = ', '.join(irr)
