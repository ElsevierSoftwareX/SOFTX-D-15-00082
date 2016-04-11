#!/usr/bin/env python
#
# Compares two trajectory files, to see if they "agree".
#
# "Agreement" means:
#    a) They contain the same number of snapshots
#    b) All pairs of snapshots have an rmsd of less than tol
import warnings

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import MDAnalysis as mda
    from MDAnalysis.analysis.align import alignto
    import sys

    topfile = sys.argv[1]
    trajfile1 = sys.argv[2]
    trajfile2 = sys.argv[3]
    tol = float(sys.argv[4])

    u1 = mda.Universe(topfile, trajfile1)
    u2 = mda.Universe(topfile, trajfile2)

    n1 = 0
    for ts in u1.trajectory:
        n1 += 1

    n2 = 0
    for ts in u1.trajectory:
        n2 += 1

    if  n1 != n2:
        print "Error: only {0} snapshots in first file but {1} in second.".format(n1,n2)
        exit(1)

    u1.trajectory.rewind()
    u2.trajectory.rewind()
    for i in range(n1):
       if i > 0:
           u1.trajectory.next()
           u2.trajectory.next()

       before,after = alignto(u1,u2)
       if after > tol:
           print "Error: at snapshot {0} the two coordinate sets have an RMSD greater than {1} (value={2}).".format(i, tol, after)
           exit(1)
