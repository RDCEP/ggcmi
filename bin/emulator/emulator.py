#!/usr/bin/env python

import argparse
import numpy as np
import os
import shutil
import sys
from netCDF4 import Dataset
from subprocess import Popen


def aggregate(inputfile, outputfile, variable, level, region, weight_var, tempdir, out_stream):
    """
    Perform aggregation
    Args:
        inputfile (string): Input filename
        outputfile (string): Output filename
        variable (string): Variable to aggregate
        level (string): Aggregation level
        region (int): Aggregation region
        weight_var (string): Name of weight variable (typically 'weight' or 'area')
        tempdir (string): Directory to store temporary intermediate files
        out_stream (file): File to redirect stdout of nco commands
    """
    with Dataset(inputfile, 'a') as f:
        vidx = np.where(f.variables['lev_' + level][:] == region)[0][0]
        latlower, latupper = f.variables['latbnds_' + level][vidx]
        lonlower, lonupper = f.variables['lonbnds_' + level][vidx]
        var = f.variables[variable]
        vararr = var[:]
        if variable == "HWAM":
            vararr[vararr == -99] = 0
        else:
            vararr = np.ma.masked_equal(vararr, -99)
        var[:] = vararr
    regselect = '-d lat,%f,%f -d lon,%f,%f' % (latlower, latupper, lonlower, lonupper)

    # average rainfed, irrigated
    options = '-h -B \'%s == %d\' -w %s -a lat,lon %s' % (level, region, weight_var, regselect)
    run_nco('ncwa', inputfile, outputfile, options, out_stream)
    with Dataset(outputfile) as f:
        udim = str(f.variables[variable].dimensions[0])
    run_nco('ncks', outputfile, outputfile, '-O -h --mk_rec_dim irr', out_stream)
    run_nco('ncpdq', outputfile, outputfile, '-O -h -a irr,%s' % udim, out_stream)
    run_nco('ncks', outputfile, outputfile, '-O -h -v %s' % variable, out_stream)

    # replace NaNs with zeros
    with Dataset(outputfile, 'a') as f:
        var = f.variables[variable]
        vararr = var[:]
        vararr = np.ma.masked_where(np.isnan(vararr), vararr)
        var[:] = vararr

    # total area
    tmpfile1 = os.path.join(tempdir, "tmpfile1.%s.%s" % (region, level))
    run_nco('ncks', inputfile, tmpfile1, '-h -v %s,%s' % (level, weight_var), out_stream)
    run_nco('ncwa', tmpfile1, tmpfile1, '-O -h -B \'%s == %d\' -N -a lat,lon %s' % (level, region, regselect),
            out_stream)
    run_nco('ncks', tmpfile1, outputfile, '-3 -h -A', out_stream)
    os.remove(tmpfile1)

    # sum
    tmpfile2 = os.path.join(tempdir, "tmpfile2.%s.%s" % (region, level))
    run_nco('ncwa', outputfile, tmpfile2, '-h -a irr -w %s' % weight_var, out_stream)
    run_nco('ncecat', tmpfile2, tmpfile2, '-O -h -u irr', out_stream)
    run_nco('ncap2', tmpfile2, tmpfile2, '-O -h -s "irr[irr]=3"', out_stream)

    # concatenate sum
    run_nco('ncks', outputfile, outputfile, '-O -h -v %s' % variable, out_stream)
    run_nco('ncrcat', '%s %s' % (outputfile, tmpfile2), outputfile, '-O -h', out_stream)
    os.remove(tmpfile2)

    # create level dimension
    run_nco('ncecat', outputfile, outputfile, '-O -h -u %s' % level, out_stream)
    run_nco('ncap2', outputfile, outputfile, '-O -h -s "%s[%s]=%d"' % (level, level, region), out_stream)

    # Create summed variable
    tmpfile3 = os.path.join(tempdir, "tmpfile3.%s.%s" % (region, level))
    run_nco('ncks', inputfile, tmpfile3, '-h -v %s,%s %s' % (weight_var, variable, regselect), out_stream)
    run_nco('ncpdq', tmpfile3, tmpfile3, '-h -O -a time,scen,irr,lat,lon', out_stream)

    with Dataset(tmpfile3, 'a') as nc:
        fill = getattr(nc.variables[variable], "_FillValue")
        irrlist = nc.variables['irr'].long_name.replace(' ', '').split(',')
        nirr = len(nc.dimensions['irr'])
        area = nc.variables['area'][:]
        prodvar = nc.createVariable("%s_prod" % variable, "f8", ("irr", "lat", "lon"),
                                    fill_value=fill)
        var = nc.variables[variable][:]
        sumval = np.ma.array(np.zeros((1, nirr+1), dtype=np.float64))

        for irridx, irr in enumerate(irrlist):
            prodvar[irridx] = area[irridx].astype(np.float64) * var[irridx].astype(np.float64)
            sumval[0, irridx] = np.nansum(prodvar[irridx])
        sumval[0, -1] = np.nansum(prodvar[:])

    with Dataset(outputfile, 'a') as nc:
        sumvar = nc.createVariable("%s_%s_sum" % (variable, level), "f8", (level, 'irr'), fill_value=fill)
        sumvar[:] = sumval[:]
    os.remove(tmpfile3)


def evaluate(val, equation, c, t, w, n):
    """
    Evaluation function for emulator
    Args:
        val (ndarray): ndarray with the shape (funcs, nlats, nlons)
        equation (list): List of terms to be evaluated
        c (float): Carbon level
        t (float): Temperature level
        w (float): Water level
        n (float): Nitrogen level
    Returns:
        ndarray representing
    """
    mods = []
    for idx, term in enumerate(equation):
        # constant: '1', '22'
        if is_float(term):
            mods.append(float(term))
        # variable: 'c', 't'
        elif len(term) == 1:
            mods.append(locals()[term])
        # exponent: c2, c23
        elif len(term) > 1 and term[0].isalpha() and is_float(term[1:]):
            mods.append(locals()[term[0]] ** float(term[1:]))
        # cross terms: cn, tw, nt
        else:
            prod = 1.0
            for termitem in term:
                prod *= locals()[termitem]
            mods.append(prod)

    for midx, m in enumerate(mods):
        val[midx] *= m
    return np.nansum(val, axis=0)


def float_or_netcdf(val, desc, minimum, maximum, parser):
    """
    Arguments can be either floats or NetCDF files.
    If an argument should be a float, cast it and check its range
    If it's a NetCDF file, open the variable, set min/max, and return a numpy array with its data
    Arguments:
        val (string): Argument to examine
        desc (string): A textual description of val
        minimum (float): Minimum acceptable value for val
        maximum (float): Maximum acceptable value for val
        parser (ArgumentParser): Parser object to throw errors to
    Returns:
        val, as either a float or a numpy array
    """
    if is_float(val):
        val = float(val)
        if minimum <= val <= maximum:
            return val
        else:
            parser.error("%s value must be between %s and %s" % (desc, minimum, maximum))
    elif is_netcdf(val):
        nc = Dataset(val)
        if desc not in nc.variables:
            parser.error("File %s does not have a variable named %s" % (val, desc))
        data = nc.variables[desc][:].astype(np.float64)
        maxval = np.nanmax(data)
        minval = np.nanmin(data)
        if maxval > maximum or minval < minimum:
            print "Limiting %s values to between %s and %s" % (val, minimum, maximum)
            data[data > maximum] = maximum
            data[data < minimum] = minimum
        nc.close()
        return data
    else:
        parser.error("%s is neither a float or a NetCDF file" % val)


def is_float(val):
    try:
        float(val)
        return True
    except (ValueError, TypeError):
        return False


def is_netcdf(f):
    try:
        Dataset(f)
        return True
    except IOError:
        return False


def process_var(var, desc, nc, nlats, nlons, c, t, w, n):
    """
    Given a variable, netcdf, and c,t,w,n levels, create a numpy array of outputs
    Args:
        var (string): Fit file variable name
        desc (string): Fit file variable description
        nc (Dataset): NetCDF Dataset
        nlats (int): Number of latitudes
        nlons (int): Number of longitudes
        c (int): Carbon level
        t (int): Temperature level
        w (int): Water level
        n (int): Nitrogen level
    Returns:
        If variable exists in nc, return a masked ndarray containing the evaluated results
        Otherwise, return an empty masked ndarray with all masked values
    """
    if var in nc.variables:
        print "Processing %s '%s'" % (desc, var)
        data = nc.variables[var][:]
        funcdim = nc.variables[var].dimensions[0]
        equation = [str(x) for x in nc.variables[funcdim].long_name.split(',')]
        data.unshare_mask()
        data = evaluate(data, equation, c, t, w, n)
    else:
        print "Unable to find %s '%s', skipping" % (desc, var)
        data = np.ma.masked_array(np.ndarray((nlats, nlons)), mask=True)
    return data


def run_nco(operator, input_file, output_file, nco_options, out_stream):
    """
    Execute an nco operator
    Args:
        operator (string): NCO operator to use
        input_file (string): Input filename
        output_file (string): Ouptut filename
        nco_options (string): NCO command line options
        out_stream (file): File to redirect stdout
    """
    command = "%s %s %s %s" % (operator, nco_options, "--output=%s" % output_file, input_file)
    p = Popen(command, shell=True, stdout=out_stream)
    exitcode = p.wait()
    if exitcode != 0:
        sys.exit("Command failed: %s" % command)


def main():
    rfvar = "fitparms"
    irrvar = "fitparmsirr"
    lats = np.arange(89.75, -90, -0.5)
    lons = np.arange(-179.75, 180, 0.5)
    irr = np.arange(2)
    mask_path = os.path.join(os.path.dirname(__file__), "masks", "aggmasks")
    np.seterr(invalid='ignore')

    parser = argparse.ArgumentParser()
    parser.add_argument("-agg", "--agg", help="Aggregation level (eg. gadm0)")
    parser.add_argument("-aggmask", "--aggmask", default="%s/gadm.nc4" % mask_path, help="Aggregation mask")
    parser.add_argument("-aggweight", "--aggweight", help="Aggregation weight")
    parser.add_argument("-aggoutput", "--aggoutput", default="output.agg.nc", help="Aggregation output filename")
    parser.add_argument("-c", "--c", required=True, help="Atmospheric co2 concentration (ppm)")
    parser.add_argument("-t", "--t", required=True,
                        help="Temperature anomalies NetCDF, or a fixed value from -1 to 6 C")
    parser.add_argument("-w", "--w", required=True,
                        help="Precipitation anomolies NetCDF, or a fixed percentage from -0.5 to 0.3")
    parser.add_argument("-n", "--n", required=True,
                        help="Fertilizer NetCDF, or fixed amount from 10 to 200 kg/ha")
    parser.add_argument("-fit", "--fit", required=True, help="Fit file")
    parser.add_argument("-o", "--o", default="output.nc", help="Output file")
    args = parser.parse_args()

    if args.agg or args.aggweight:
        if not args.agg or not args.aggweight:
            parser.error("Aggregation requires setting --agg and --aggweight")

    c = float_or_netcdf(args.c, "co2", 360, 510, parser)
    t = float_or_netcdf(args.t, "temperature", -1.0, 6.0, parser)
    w = float_or_netcdf(args.w, "precipitation", -0.5, 0.3, parser)
    n = float_or_netcdf(args.n, "fertilizer", 10.0, 200.0, parser)

    fit = args.fit
    if not is_netcdf(fit):
        parser.error("Fit file %s is not a NetCDF file" % fit)

    fitnc = Dataset(fit)
    out_var = fitnc.variable if 'variable' in fitnc.ncattrs() else "result"
    long_name = fitnc.long_name if 'long_name' in fitnc.ncattrs() else None
    units = fitnc.units if 'units' in fitnc.ncattrs() else None

    rfdata = process_var(rfvar, "rainfed variable", fitnc, lats.size, lons.size, c, t, w, n)
    irrdata = process_var(irrvar, "irrigated variable", fitnc, lats.size, lons.size, c, t, w, n)
    result = np.ma.array(np.ndarray((irr.size, lats.size, lons.size)))
    result[0] = rfdata[:]
    result[1] = irrdata[:]
    result = np.ma.masked_where(result < 0, result)

    with Dataset(args.o, "w") as nc:
        nc.createDimension("lat", lats.size)
        nc.createDimension("lon", lons.size)
        nc.createDimension("irr", irr.size)
        nc.createVariable("irr", "i4", "irr")[:] = np.arange(irr.size)
        nc.variables['irr'].long_name = "rf, irr"
        latvar = nc.createVariable("lat", "f4", "lat")
        latvar[:] = lats[:]
        latvar.axis = "Y"
        latvar.long_name = "latitude"
        latvar.standard_name = "latitude"
        latvar.units = "degrees_north"
        lonvar = nc.createVariable("lon", "f4", "lon")
        lonvar[:] = lons[:]
        lonvar.axis = "X"
        lonvar.long_name = "longitude"
        lonvar.standard_name = "longitude"
        lonvar.units = "degrees_east"
        var = nc.createVariable(out_var, "f4", ("irr", "lat", "lon"), fill_value=1e20, zlib=True, complevel=9)
        var[:] = result[:]
        if units:
            var.units = units
        if long_name:
            var.long_name = long_name

    if args.agg:
        tempdir = os.path.join(os.getcwd(), 'agg.temp')
        if os.path.exists(tempdir):
            shutil.rmtree(tempdir)
        os.makedirs(tempdir)
        merged = os.path.join(tempdir, "merged.nc")
        shutil.copyfile(args.o, merged)
        mergedtmp = "%s.tmp" % merged
        devnull = open(os.devnull, 'w')
        for f in [args.aggweight, args.aggmask]:
            options = '-h'
            run_nco('ncks', f, "%s.tmp" % merged, options, devnull)
            run_nco('ncks', "%s.tmp" % merged, merged, '-h -A', devnull)
            os.remove(mergedtmp)

        regions = {}
        levels = args.agg.split(',')
        with Dataset(merged) as f:
            for lev in levels:
                regions[lev] = f.variables['lev_' + lev][:]

        # aggregate each level
        for lev in levels:
            reglist = list(set(regions[lev].tolist()))
            for regidx, reg in enumerate(reglist):
                sys.stdout.write('\rAggregating %s region %d/%d' % (lev, regidx+1, len(reglist)))
                sys.stdout.flush()
                tmpfile = os.path.join(tempdir, 'tempfile.%s.%09d' % (lev, reg))
                aggregate(merged, tmpfile, out_var, lev, reg, 'area', tempdir, devnull)
            print

        # append results
        if os.path.isfile(args.aggoutput):
            os.remove(args.aggoutput)
        for i in range(len(levels)):
            tmpfile = os.path.join(tempdir, 'tempfile.%s' % levels[i])
            run_nco('ncrcat', os.path.join(tempdir, 'tempfile.%s.*' % levels[i]), tmpfile, '-h', devnull)
            run_nco('ncpdq', tmpfile, tmpfile, '-O -h -a %s' % levels[i], devnull)
            run_nco('ncrename', tmpfile, tmpfile, '-O -h -v %s,%s_%s' % (out_var, out_var, levels[i]),
                    devnull)
            run_nco('ncks', tmpfile, args.aggoutput, '-h -A', devnull)
        shutil.rmtree(tempdir)

    print "Done"


if __name__ == "__main__":
    main()
