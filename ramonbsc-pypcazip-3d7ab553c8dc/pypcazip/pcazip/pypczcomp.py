#!/usr/bin/env python

import sys
import logging as log

import numpy as np

import pcazip
from MDPlus.analysis.pca import pczfile
from MDPlus.core.cofasu import rmsd


def pczcomp(args):
    '''
    Pczcomp compares two .pcz files. It reports on the RMSD between the
    two average structures, the dot product matrix, the subspace overlap,
    and the average maximum dot product. In addition (unless the '--quick'
    option is given) it reports a number of metrics related to Mahalanobis
    distances. Firstly it reports the average Mahalanobis distance of the
    snapshots in each .pcz file from their respective averages, and from
    each other's average structures. Then it reports on the percentage of
    snapshots in each .pcz file that lie within the 90% envelope of the
    snapshots in the other.
    '''

    if args.verbosity > 0:
        log.basicConfig(format="%(levelname)s: %(message)s", level=log.DEBUG)
        log.info("Verbose output.")
    else:
        log.basicConfig(format="%(levelname)s: %(message)s")
     
    listLong = ['--input', '--quick', '--nvecs']
    listShort = ['-i']
    	
    for i in range(1, len(sys.argv)):
        if sys.argv[i] in listShort and listLong[listShort.index(sys.argv[i])] in sys.argv:
            log.error('''Please, use either the long or the short form of an option but never both! Try again!''')
            sys.exit(-1)
	
    if (args.input is None):
        log.error('')
        log.error('''All or any of the mandatory command line arguments is missing. The correct usage of pyPczcomp should be:''')
        log.error('pyPczcomp -i|--input <input-file> [optional arguments]')
        log.error('')
        log.error('''Type "pyPczcomp -h" or "pyPczcomp --help" for further details.''')
        log.error('')
        sys.exit(-1)

    px = pczfile.Pczfile(args.input[0])
    py = pczfile.Pczfile(args.input[1])

    if px.natoms != py.natoms:
        print 'Error: the number of atoms in the two files is different.'
        exit(1)

    if px.nvecs < args.nvecs:
        print 'Error: {0} only contains {1} vectors'.format(args.input[0],px.nvecs)
        exit(1)

    if py.nvecs < args.nvecs:
        print 'Error: {0} only contains {1} vectors'.format(args.input[1],py.nvecs)
        exit(1)

    natoms = px.natoms

    print 'Comparison of X: {0} and Y: {1}'.format(args.input[0],args.input[1])

    rms,R = rmsd(px.avg(),py.avg(),R=True)
    print 'Rmsd between <X> and <Y>: {0:6.2f}'.format(rms)

    if not args.quick:
        print 'Mahalanobis distances:'
# mxx, myy, mxy and myx will contain the Mahalanobis distances of each snapshot
# in X from <X>, of Y from <Y>, of X from <Y>, and of Y from <X>.
        mxx = np.zeros(px.nframes)
        myy = np.zeros(px.nframes)
        myx = np.zeros(py.nframes)
        mxy = np.zeros(px.nframes)
        evalx = px.evals()
        evaly = py.evals()
        for k in range(py.nframes):
            xydif = py.frame(k)
            syy = py.scores(k)
            syx = px.map(xydif)
            smyy = 0.0
            smyx = 0.0
            for j in range(args.nvecs):
                smyx += syx[j]*syx[j]/evalx[j]
                smyy += syy[j]*syy[j]/evaly[j]
            myx[k] = np.sqrt(smyx)
            myy[k] = np.sqrt(smyy)

        for k in range(px.nframes):
            xydif = px.frame(k)
            sxx = px.scores(k)
            sxy = py.map(xydif)
            smxx = 0.0
            smxy = 0.0
            for j in range(args.nvecs):
                smxy += sxy[j]*sxy[j]/evaly[j]
                smxx += sxx[j]*sxx[j]/evalx[j]
            mxy[k] = np.sqrt(smxy)
            mxx[k] = np.sqrt(smxx)

        print 'X from <X>: {0:6.2f} +/- {1:6.2f}'.format(mxx.mean(),mxx.std())
        print 'Y from <Y>: {0:6.2f} +/- {1:6.2f}'.format(myy.mean(),myy.std())
        print 'X from <Y>: {0:6.2f} +/- {1:6.2f}'.format(mxy.mean(),mxy.std())
        print 'Y from <X>: {0:6.2f} +/- {1:6.2f}'.format(myx.mean(),myx.std())

#    
# Calculation of the 90% enevlopes
#
        tx = np.sort(mxx)[int(0.9*px.nframes)]
        ty = np.sort(myy)[int(0.9*py.nframes)]

        i = 0
        for j in mxy:
            if j < tx:
                i += 1
        print '{0:3d}% of X lies within the 90% envelope of Y'.format(i*100/px.nframes)

        i = 0
        for j in myx:
            if j < ty:
                i += 1
        print '{0:3d}% of Y lies within the 90% envelope of X'.format(i*100/py.nframes)

#
# Dot product stuff. We rotate the eigenvectors in the second (Y) .pcz file
# using the rotation matrix found in least-squares fitting the average
# structures in the earlier RMSD calculation.
#
    evx = np.zeros((args.nvecs,3*natoms))
    evy = np.zeros((args.nvecs,3*natoms))

    for i in range(args.nvecs):
        evx[i] = px.evec(i)
        evy[i] = py.evec(i)
        evy[i] = np.dot(evy[i].reshape(natoms,3),R).flatten()

    dp = np.zeros((args.nvecs,args.nvecs))
    ss = 0.0
    for i in range(args.nvecs):
        for j in range(args.nvecs):
            dp[i,j] = np.absolute(np.dot(evx[i],evy[j]))
            ss += dp[i,j]*dp[i,j]

    ss = ss/args.nvecs

    print 'Dot product matrix:'
    print '    ',
    for i in range(args.nvecs):
        print '{0:4d}'.format(i),
    print ' '
    for i in range(args.nvecs):
        print '{0:4d}'.format(i),
        for j in range(args.nvecs):
            print '{0:4.2f}'.format(dp[i,j]),
        print ' '

    print 'Subspace overlap: {0:6.3f}'.format(np.sqrt(ss))
    amdp = (dp.max(axis=0).sum()+dp.max(axis=1).sum())/(2*args.nvecs)
    print 'Average maximum dot product: {0:6.3f}'.format(amdp)
