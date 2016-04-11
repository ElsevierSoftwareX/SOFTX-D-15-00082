#!/usr/bin/env python -W ignore
'''
                 *** The command line interface for pyPczclust ***

                       Adapted to use the mapping module.
'''

import logging as log

import numpy as np
from scipy import ndimage

from MDPlus.analysis.pca import pczfile
from MDPlus.analysis import mapping

def pczclust(args): 
    '''
    Performs histogram/watershed based clustering on data from a .pcz
    file.
    '''

    if args.verbosity > 0:
        log.basicConfig(format="%(levelname)s: %(message)s", level=log.DEBUG)
        log.info("Verbose output")
    else:
        log.basicConfig(format="%(levelname)s: %(message)s")

    p = pczfile.Pczfile(args.pczfile)

    projs = np.zeros((p.nframes,args.dims))
    for i in range(args.dims):
        projs[:,i] = p.proj(i)

    m = mapping.Map(projs,resolution=args.bins, boundary=1)
    mapping.watershed(m)
    out = [m.cluster_id(id) for id in projs]

    np.savetxt(args.outfile,np.c_[projs,out], fmt=("%8.3f"*args.dims + "%5d"))
