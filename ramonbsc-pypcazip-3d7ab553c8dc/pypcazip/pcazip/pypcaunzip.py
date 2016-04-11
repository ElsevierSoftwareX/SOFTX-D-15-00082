#!/usr/bin/python

import sys
import os.path as op
import logging as log

import numpy as np

import pcazip
from MDPlus import Universe, Writer
from MDPlus.analysis.pca import pczfile

#######################################################
# MAIN FUNCTION
#######################################################

def pcaunzip(args):

    if args.verbosity > 0:
        log.basicConfig(format="%(levelname)s: %(message)s", level=log.INFO)
        log.info("Verbose output.")
    else:
        log.basicConfig(format="%(levelname)s: %(message)s")

    if ((args.topology is None) or (args.compressed is None) or (args.preload is None)):
        log.error('')
        log.error(
            'All or any of the mandatory command line arguments is missing. The correct usage of PCAUNZIP should be:')
        log.error(
            'python ./pcaunzip.py -t|--topology <topology-file>, -c|--compressed <compressed-file>, -p|--preload <True|False>, [optional arguments]')
        log.error('')
        log.error('Type "python ./pcaunzip.py -h" or "python ./pcaunzip.py --help" for further details.')
        log.error('')
        sys.exit(-1)

    # If no name is provided for the output file, the extension of the npz file is just changed into .dcd.
    if not args.output:
        dir = op.dirname(args.compressed)
        base = op.basename(args.compressed)
        name = op.splitext(base)[0]
        args.output = op.join(dir, name + ".dcd")

    try:
        import netCDF4
        nonetCDF4 = False
    except ImportError:
        nonetCDF4 = True

    if op.splitext(args.output)[1].lower() == '.ncdf' and nonetCDF4:
        log.error('netcdf4-python with the netCDF4 and HDF5 libraries must be installed to read AMBER .ncdf files.\nSee installation instructions at https://code.google.com/p/mdanalysis/wiki/netcdf')
        exit(1)

    # PCAunzip
    log.info("PCAunzipping")
    pfile = pczfile.Pczfile(args.compressed, args.preload)
    kwargs = {'dt':1.0}
    w = Writer(args.output, pfile.natoms, **kwargs)

    for ts_index in xrange(pfile.nframes):
        time_step = []
        time_step.append(pfile.frame(ts_index))
        kwargs = {'nparray':np.array(time_step), 'dt':1.0}
        universe = Universe(args.topology, 'whatever.Numpy', **kwargs)
        for ts in universe.trajectory:
            w.write(ts)
    w.close()

