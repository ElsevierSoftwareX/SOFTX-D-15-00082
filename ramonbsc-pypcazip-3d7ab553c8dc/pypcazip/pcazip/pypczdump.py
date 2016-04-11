#!/usr/bin/python
'''																	
 				Module including functionality that allows to elaborate					
 					and extract useful information from pcz files.						
'''

import sys
import logging as log

import numpy as np
import math

import pcazip
from MDPlus.analysis.pca import pczfile

def pczdump(args):

    listLong = ['--input', '--evec', '--fluc', '--evals', '--avg', '--coll', '--proj', '--info', '--verbosity',
                '--output', '--rms', '--anim', '--pdb']
    listShort = ['-i', '-e', '-f', '-l', '-a', '-c', '-p', '-n', '-v', '-o', '-r', '-m']

    numArgs = 0
    for i in range(1, len(sys.argv)):
        if sys.argv[i] in listShort and listLong[listShort.index(sys.argv[i])] in sys.argv:
            log.error('''Please, use either the long or the short form of an option but never both! Try again!''')
            sys.exit(-1)

    listUniqOpts = ['--evec', '--fluc', '--evals', '--avg', '--coll', '--proj', '--info', '--rms', '--anim', '-e', '-f',
                    '-l', '-a', '-c', '-p', '-n', '-r', '-m']
    for i in range(1, len(sys.argv)):
        if sys.argv[i] in listUniqOpts:
            numArgs += 1
            if numArgs == 2:
                log.error(
                    '''Please ask for an option at a time from pyPczdump. You have currently been asking for more than one option! Try again!''')
                sys.exit(-1)

    if args.verbosity > 0:
        log.basicConfig(format="%(levelname)s: %(message)s", level=log.DEBUG)
        log.info("Verbose output.")
    else:
        log.basicConfig(format="%(levelname)s: %(message)s")

    if (args.input is None):
        log.error('')
        log.error('''All or any of the mandatory command line arguments is missing. The
		correct usage of pyCoCo should be:''')
        log.error('python ../pczdump.py -i|--input <input-file> [optional arguments]')
        log.error('')
        log.error('''Type "pyPczdump -h" or "pyPczdump --help" for further details.''')
        log.error('')
        sys.exit(-1)

    if args.evec is not None:
        pcz_instance = pczfile.Pczfile(args.input)
        e = pcz_instance.evec(int(args.evec))
        if args.output is not None:
            np.savetxt(args.output, np.column_stack((e,)))
        else:
            log.info("Values of eigenvector %d follow:\n" % int(args.evec))
            for i in e:
                print i

    elif args.proj is not None:
        pcz_instance = pczfile.Pczfile(args.input)
        p = pcz_instance.proj(int(args.proj))
        if args.output is not None:
            np.savetxt(args.output, np.column_stack((p,)))
        else:
            log.info("Values of the projections of eigenvector %d follow:\n" % int(args.proj))
            for i in p:
                print i
    elif args.evals is True:
        pcz_instance = pczfile.Pczfile(args.input)
        evs = pcz_instance.evals()
        if args.output is not None:
            np.savetxt(args.output, np.column_stack((evs,)))
        else:
            log.info("The list of eigenvalues follows:\n")
            for i in evs:
                print i
    elif args.info is True:
        pcz_instance = pczfile.Pczfile(args.input)
        if pcz_instance.title is not None:
            log.info("Basic information on the compressed file:\n")
            print "%s\n" % pcz_instance.title
        print "The PCZ file format version is %s\n" % pcz_instance.version
        print "The number of atoms is %d\n" % pcz_instance.natoms
        print "The number of frames is %d.\n" % pcz_instance.nframes
        print "The number of eigenvectors is %d.\n" % pcz_instance.nvecs
        print "The variance captured inside this file is %.2f%%.\n" % (
        100 * pcz_instance.evals().sum() / pcz_instance.quality)
    elif args.fluc is not None:
        pcz_instance = pczfile.Pczfile(args.input)
        evec = pcz_instance.evec(int(args.fluc))
        fluc = np.empty((pcz_instance.natoms))
        for i in range(pcz_instance.natoms):
            j = 3 * (i + 1) - 3
            fluc[i] = evec[j] * evec[j] + evec[j + 1] * evec[j + 1] + evec[j + 2] * evec[j + 2]
            fluc[i] = np.sqrt(fluc[i])
        if args.output is not None:
            np.savetxt(args.output, np.column_stack((fluc,)))
        else:
            log.info("Values of the fluctuations associated with eigenvector %d follow:\n" % int(args.fluc))
            for i in fluc:
                print i
    elif args.avg is True:
        pcz_instance = pczfile.Pczfile(args.input)
        avg = pcz_instance.avg()

        atom_nr = 0

        if args.pdb is None:
            log.error(
                '''The "--avg" or "-a" option REQUIRES the "--pdb" option. The output average structure will be a pdb file either on disk or on screen! Try again including the "--pdb" option!''')
            sys.exit(-1)
        else:
            if args.output is not None:
                with open(args.output, 'w') as wf:
                    with open(args.pdb) as f:
                        for line in f:
                            if line[:4] == "ATOM":
                                wf.write(line[:30])
                                wf.write('%8.3f%8.3f%8.3f' % (avg[atom_nr][0], avg[atom_nr][1], avg[atom_nr][2]))
                                wf.write(line[54:])
                                atom_nr += 1
                            else:
                                wf.write(line)
                    f.close()
                wf.close()
            else:
                with open(args.pdb) as f:
                    for line in f:
                        if line[:4] == "ATOM":
                            print line[:30] + '{0:8.3f}{1:8.3f}{2:8.3f}'.format(avg[atom_nr][0], avg[atom_nr][1],
                                                                                avg[atom_nr][2]) + line[54:-1]
                            atom_nr += 1
                        else:
                            print line[:-1]
                f.close()

    elif args.coll is True:
        pcz_instance = pczfile.Pczfile(args.input)
        col = np.zeros((pcz_instance.nvecs))
        for i in range(pcz_instance.nvecs):
            e = pcz_instance.evec(i)
            for j in range(0, pcz_instance.natoms * 3, 3):
                r2 = e[j] * e[j] + e[j + 1] * e[j + 1] + e[j + 2] * e[j + 2]
                col[i] = col[i] - r2 * math.log(r2)
        col = np.exp(col) / pcz_instance.natoms
        if args.output is None:
            # np.savetxt(args.input+"_collMetric.out",np.column_stack((col,)))
            log.info("Collectivity metric values K for each eigenvector follow.:\n")
            log.info("Modes producing most of the collective motion in the system will have high K values.\n")
            for i in col:
                print i
        else:
            np.savetxt(args.output, np.column_stack((col,)))

    elif args.rms is not None:
        pcz_instance = pczfile.Pczfile(args.input)
        if ((int(args.rms) >= 0) and (int(args.rms) < pcz_instance.nframes)):
            s2 = pcz_instance.scores(int(args.rms))
        rmsd = np.zeros((pcz_instance.nframes))
        for i in range(pcz_instance.nframes):
            rmsd[i] = 0.0
            s1 = pcz_instance.scores(i)
            if ((int(args.rms) >= 0) and (int(args.rms) < pcz_instance.nframes)):
                s1 = s1 - s2
            rmsd[i] = np.sum(s1 * s1) / pcz_instance.natoms

        if args.output is not None:
            np.savetxt(args.output, np.column_stack((rmsd,)))
        else:
            for i in rmsd:
                print np.sqrt(i)

    elif args.anim is not None:
        pcz_instance = pczfile.Pczfile(args.input)
        avg = pcz_instance.avg()
        evec = pcz_instance.evec(int(args.anim))
        proj = pcz_instance.proj(int(args.anim))

        rmin = np.min(proj)
        rmax = np.max(proj)
        rinc = (rmax - rmin) * 0.1
        # re-create proj...
        proj = np.zeros(20)
        for i in range(1, 6, 1):
            proj[i] = proj[i - 1] + rinc
        for i in range(6, 16, 1):
            proj[i] = proj[i - 1] - rinc
        for i in range(16, 20, 1):
            proj[i] = proj[i - 1] + rinc

        if args.pdb is None:
            log.error(
                '''The "--anim" or "-m" option REQUIRES the "--pdb" option. The output average structure will be a pdb file either on disk or on screen! Try again including the "--pdb" option!''')
            sys.exit(-1)
        else:
            if args.output is not None:
                with open(args.output, 'w') as wf:
                    wf.write('REMARK   Animation of eigenvector %4d\n' % int(args.anim))
                    with open(args.pdb) as f:
                        line = f.readline()
                        prec_curs = 0
                        while (line[:4] != "ATOM"):
                            if line[:5] != "MODEL":
                                wf.write(line)
                            prec_curs = f.tell()
                            line = f.readline()

                        for j in range(20):
                            f.seek(prec_curs, 0)
                            x = avg.reshape((avg.shape[0] * avg.shape[1])) + evec * proj[j]
                            atom_nr = 0
                            wf.write('MODEL  %5d\n' % (j + 1))
                            for line in f:
                                if line[:4] == "ATOM":
                                    wf.write(line[:30])
                                    for m in range(3):
                                        wf.write('%8.3f' % (x[atom_nr + m]))
                                    wf.write(line[54:])
                                    atom_nr += 3
                                elif (line[:3] != "END"):
                                    wf.write(line)
                            wf.write('ENDMDL\n')
                    f.close()
                wf.close()
            else:
                print 'REMARK   Animation of eigenvector {0:4d}'.format(int(args.anim))
                with open(args.pdb) as f:
                    line = f.readline()
                    prec_curs = 0
                    while (line[:4] != "ATOM"):
                        if line[:5] != "MODEL":
                            print line[:-1]
                        prec_curs = f.tell()
                        line = f.readline()

                    for j in range(20):
                        f.seek(prec_curs, 0)
                        x = avg.reshape((avg.shape[0] * avg.shape[1])) + evec * proj[j]
                        atom_nr = 0
                        print 'MODEL  {0:5d}'.format(j + 1)
                        for line in f:
                            if line[:4] == "ATOM":
                                print line[:30] + '{0:8.3f}{1:8.3f}{2:8.3f}'.format(x[atom_nr], x[atom_nr + 1],
                                                                                    x[atom_nr + 2]) + line[54:-1]
                                atom_nr += 3
                            elif (line[:3] != "END"):
                                print line[:-1]
                        print 'ENDMDL'
                f.close()
