#!/usr/bin/python
"""
Routines to parse pcz format files. All routines return numpy
arrays, if appropriate.
"""
import struct
import MDAnalysis.lib.qcprot as qcp
import numpy as np
import sys
try:
    import h5py
    h5py_available = True
except ImportError:
    h5py_available = False

class Pczfile:
    def __init__(self, filename, preload=False):
        """
        Initialises a new pcz object with the data from the given
        file. With the optional keyword argument 'preload' set to the
        default value of False, data is extracted directly from the pcz
        file as and when required. If memory is not an issue, setting
        'preload=True' will read all data into private arrays and certain
        subsequent methods calls will be much faster.
        """
        self.filename = filename
        try:
            self.filehandle = open(self.filename, 'rb')
        except IOError as e:
            print "Problems while tried to open a pcz-format file."
            print "I/O error({0}): {1}".format(e.errno, e.strerror)
            return
        self.data = self.filehandle.read(100)
        if self.data[0:4] != 'PCZ4' and self.data[0:4] != 'PCZ6':
                if not h5py_available:
                    print "Error: h5py module not available. Please install the PCZ7 extra requirements to be able to use it. The command to do this requirements is: pip install pyPcazip[PCZ7]"
                    sys.exit(0)
                try:
                    self.filehandle = h5py.File(self.filename, "r")
                    self.data = 'PCZ7'
                except IOError:
                    print 'Error - unrecognised file format: ', self.data[0:4], 'while trying to open a pcz file!'
                    self.filehandle.close()
                    return

        if self.data[0:4] == 'PCZ4' or self.data[0:4] == 'PCZ6':
            self.keydata = struct.unpack('4s80s3if', self.data)
            self.version = self.keydata[0]
            self.title = self.keydata[1]
            self.natoms = self.keydata[2]
            self.nframes = self.keydata[3]
            self.nvecs = self.keydata[4]
            self.quality = self.keydata[5]

            self.data = self.filehandle.read(16)
            self.keydata = struct.unpack('4i', self.data)
            self.extra = self.keydata[3]

            if self.extra == 0:
                self.headerblocksize = 116
            else:
                self.headerblocksize = 116 + 16 * self.natoms

            self.avgblocksize = 3 * self.natoms * 4
            if self.version == 'PCZ4':
                self.evblocksize = 3 * self.natoms * 4 + 4 + self.nframes * 4
            else:
                self.evblocksize = 3 * self.natoms * 4 + 4 + 8 + self.nframes * 2

            self.preloaded = preload
            if preload:
                self.preloaded = False
                self._avg = self.avg()
                self._evals = self.evals()
                self._evecs = self.evecs()
                self._projs = np.zeros((self.nvecs, self.nframes))
                for i in range(self.nvecs):
                    self._projs[i, :] = self.proj(i)
                self.preloaded = True

        elif self.data[0:4] == 'PCZ7':
            self.version = self.filehandle.attrs['version']
            self.title = self.filehandle.attrs['title']
            self.natoms = self.filehandle.attrs['natoms']
            self.nframes = self.filehandle.attrs['nframes']
            self.numframes = self.nframes
            self.n_evecs = self.filehandle.attrs['nvecs']
            self.nvecs = self.filehandle.attrs['nvecs']
            self.quality = self.filehandle.attrs['quality']
            self._evecs = np.array(self.filehandle['evec_dataset'])
            self._evals = np.array(self.filehandle['eval_dataset'])
            self._avg = np.array(self.filehandle['ref_coord_dataset'])

            self.preloaded = preload
            if preload:
                self._projs = np.array(self.filehandle['proj_dataset'])
                self._p0 = np.array(self.filehandle['p0_dataset'])
                self._pinc = np.array(self.filehandle['pinc_dataset'])
        else:
            print 'Error - unrecognised file format: ', self.data[0:4], 'while trying to open a pcz file!'
            return

    def avg(self):
        """
        Returns the average structure contained in the pcz file
        as an (natoms,3) numpy array.
        """
        if self.preloaded:
            return self._avg
        elif self.version == 'PCZ7':
            return self._avg.reshape((self.natoms, 3))
        else:
            self.filehandle.seek(self.headerblocksize)
            self.data = self.filehandle.read(self.avgblocksize)
            return np.array(struct.unpack(str(self.natoms * 3) + 'f',
                                          self.data)).reshape((self.natoms, 3))

    def eval(self, ival):
        """
        Returns an eigenvalue from the file.
        """
        if ival >= self.nvecs:
            print 'Error - only ', self.nvecs, ' eigenvectors present'
            return 0.0
        else:
            if self.preloaded:
                return self._evals[ival]
            elif self.version == 'PCZ7':
                return self._evals[ival]
            else:
                self.filehandle.seek(self.headerblocksize + self.avgblocksize +
                                     ival * self.evblocksize + self.avgblocksize)
                return struct.unpack('f', self.filehandle.read(4))[0]

    def evals(self):
        """
        Returns an array of all eigenvalues in the file.
        """
        if self.preloaded:
            return self._evals
        else:
            evs = np.zeros(self.nvecs)
            for i in range(self.nvecs):
                evs[i] = self.eval(i)
            return evs

    def evec(self, ivec):
        """
        Returns a chosen eigenvector from the file in the
        form of a (3*natoms) numpy array.
        """
        if ivec >= self.nvecs:
            print 'Error - only ', self.nvecs, 'eigenvectors present'
            return None
        else:
            if self.preloaded:
                return self._evecs[ivec, :]
            elif self.version == 'PCZ7':
                return self._evecs[ivec, :]
            else:
                self.filehandle.seek(self.headerblocksize + self.avgblocksize +
                                     ivec * self.evblocksize)
                return np.array(struct.unpack(str(3 * self.natoms) + 'f',
                                              self.filehandle.read(self.avgblocksize)))

    def evecs(self):
        """
        Returns all eigenvectors in the file in the form of a
        (nvecs,3*natoms) numpy array.
        """
        if self.preloaded:
            return self._evecs
        else:
            evs = np.zeros((self.nvecs, self.natoms * 3))
            for i in range(self.nvecs):
                evs[i, :] = self.evec(i)
            return evs

    def proj(self, iproj):
        """
        Returns an array of the projections along a given eigenvector. There
        will be one value per snapshot.
        """
        if iproj >= self.nvecs:
            print 'Error - only ', self.nvecs, 'eigenvectors present'
            return None
        else:
            if self.preloaded:
                if self.version == 'PCZ7':
                    proj_list = []
                    for frame_index in xrange(self.nframes):
                        proj_list.append(self._p0[frame_index] + self._pinc[frame_index] * self._projs[frame_index][iproj])
                return self._projs[iproj, :]

            elif self.version == 'PCZ7':
                proj_list = []
                for frame_index in xrange(self.nframes):
                    p0 = np.array(self.filehandle['p0_dataset'][frame_index])
                    pinc = np.array(self.filehandle['pinc_dataset'][frame_index])
                    i = np.array(self.filehandle['proj_dataset'][frame_index][iproj])
                    proj_list.append(p0 + pinc * i)
                return np.array(proj_list)
            else:
                self.filehandle.seek(self.headerblocksize + self.avgblocksize +
                                     iproj * self.evblocksize + self.avgblocksize + 4)
                if self.version == 'PCZ4':
                    return np.array(struct.unpack(str(self.nframes) + 'f',
                                                  self.filehandle.read(4 * self.nframes)))
                else:
                    ip = list(0 for i in range(self.nframes))
                    p0, pinc = struct.unpack('2f', self.filehandle.read(8))
                    self.filehandle.seek(self.headerblocksize + self.avgblocksize +
                                         iproj * self.evblocksize + self.avgblocksize + 12)
                    ip = struct.unpack(str(self.nframes) + 'h',
                                       self.filehandle.read(self.nframes * 2))
                    return np.array(list(p0 + pinc * i for i in ip))

    def scores(self, framenumber):
        """
        Method that returns the scores (projections) corresponding to
        a chosen snapshot (zero based).
        """
        if framenumber >= self.nframes:
            return None
        elif self.version == 'PCZ7':
            p0 = np.array(self.filehandle['p0_dataset'][framenumber])
            pinc = np.array(self.filehandle['pinc_dataset'][framenumber])
            projs_comp = np.array(self.filehandle['proj_dataset'][framenumber])
            projs = [p0 + pinc * i for i in projs_comp]
            return np.array(projs)
        else:
            x = np.zeros(self.nvecs)
            for i in range(self.nvecs):
                x[i] = self.proj(i)[framenumber]
            return x

    def frame(self, framenumber):
        """
        Method to return the coordinates of the given frame - i.e.
        to decompress a snapshot. The data is returned as a (natoms,3)
        numpy array.
        """
        if (framenumber >= self.nframes):
            return None
        else:
            scores = self.scores(framenumber)
            return self.unmap(scores)

    def unmap(self, scores):
        """
        Method to return the coordinates corresponding to a given
        set of scores. If the scores vector has less than nvec elements,
        Pczfile.expand is called to fill in the missing values.
        """
        if self.version == 'PCZ7':
            projs = np.array(scores).T
            return (np.array((self._avg + np.dot(projs, self._evecs)))).reshape(self.natoms, 3)
        else:
            x = self.avg()
            if len(scores) < self.nvecs:
                scores = self.expand(scores)
            for i in range(self.nvecs):
                x = x + (self.evec(i) * scores[i]).reshape((self.natoms, 3))
            return x

    def expand(self, scores):
        """
        Method to complete a truncated list of scores with values for
        missing eigenvectors. Basically the pcz file is scanned to
        find the frame with the best match to the scores given, the
        missing values are then copied directly from this frame.
        """
        nvecs = self.nvecs
        ns = len(scores)
        if ns >= nvecs:
            return scores
        else:
            closestframe = self.closest(scores)
            newscores = self.scores(closestframe)
            newscores[0:ns] = scores
            return newscores

    def closest(self, scores):
        """
        Method to find the index of the frame with scores closest to the
        target values.
        """
        ns = len(scores)
        if self.preloaded:
            temp = self._projs
        else:
            temp = np.zeros((self.nvecs, self.nframes))
            for vec in range(self.nvecs):
                temp[vec, :] = self.proj(vec)

        best = 0
        err = ((temp[0:ns, 0] - scores) * (temp[0:ns, 0] - scores)).sum()
        for frame in range(self.nframes):
            newerr = ((temp[0:ns, frame] - scores) * (temp[0:ns, frame] - scores)).sum()
            if newerr < err:
                err = newerr
                best = frame
        return best

    def map(self, crds):
        """
        Method to map an arbitrary coordinate set onto the PC model. The
        coordinate set should be a (natom,3) array-like object that matches
        (for size) what's in the pczfile. An array of projections will be
        returned, one value for each eignevector in the pcz file.
        """
        crdset = np.array(crds)
        if np.shape(crdset) == (self.natoms, 3):
            avg = np.reshape(self.avg(), (self.natoms, 3))
            rot = np.zeros(9, dtype=np.float64)
            crdset_cog = crdset.sum(axis=0) / crdset.shape[0]
            avg_cog = avg.sum(axis=0) / avg.shape[0]
            c = crdset.copy()
            c -= crdset_cog
            rmsd = qcp.CalcRMSDRotationalMatrix((avg - avg_cog).T.astype(np.float64),
                                                c.T.astype(np.float64),
                                                crdset.shape[0], rot, None)
            R = np.matrix(rot.reshape(3, 3))
            c = c * R
            c += avg_cog
            c -= avg
            prj = np.zeros(self.nvecs)
            for i in range(self.nvecs):
                prj[i] = (np.dot(c.flatten(), self.evec(i)))
            return prj
        else:
            return None

