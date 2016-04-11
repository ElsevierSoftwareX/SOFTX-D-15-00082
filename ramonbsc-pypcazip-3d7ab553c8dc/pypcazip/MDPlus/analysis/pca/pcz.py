#!/usr/bin/env python
'''
The pcz class.
'''
from MDAnalysis.analysis.align import *
import numpy as np
from scipy.linalg import *
import struct
import logging as log
import sys
from time import time
import warnings

warnings.simplefilter("ignore")

try:
    import h5py
    h5py_available = True
except ImportError:
    h5py_available = False

class Pcz:
    def __init__(self, cofasu, version='PCZ6',target=None, covar=None, quality=90.0,
                 req_evecs=None, rank = 0, preload=False, fastmethod=False):
        '''
        Initialises a new pcz object with the data from the given
        cofasu object:
        >>> topfile = "../../../../test/2ozq.pdb"
        >>> trjfile = "../../../../test/2ozq.dcd"
        >>> from MDPlus.core import cofasu
        >>> c = cofasu.Cofasu(cofasu.Fasu(topfile, trjfile, filter='name CA and resid 1-10'))
        >>> p = Pcz(c)

        
        target can be a precalculated global average
        structure:

        >>> target = c.fitted_average()-10.0
        >>> p = Pcz(c, target=target)

        covar can be a precalculated covariance matrix, in which case the
        corresponding target structure must also be given:
        >>> cm = c.cov()
        >>> p = Pcz(c, target=target, covar=cm)

        >>> p = Pcz(c, covar=cm)
        Traceback (most recent call last):
           ...
        ValueError: A defined covariance matrix requires a defined target.


        The quality setting defaults to 90%:
        >>> p = Pcz(c, quality=95)
        >>> p = Pcz(c, quality=120)
        Traceback (most recent call last):
           ...
        ValueError: quality must lie in the range 0:100.

        If preload==True, then an attempt will be made to load all the
        trajectory data into memory. This means slower initialisation but
        faster access to data later.

        >>> c = cofasu.Cofasu(cofasu.Fasu(topfile, trjfile, filter='name CA and resid 1-10'))
        >>> p = Pcz(c)
        >>> ev1 = p.evals()
        >>> p2 = Pcz(c, preload=True)
        >>> ev2 = p2.evals()
        >>> print(np.allclose(ev1, ev2))
        True
 
        If fastmethod=True then a fast approximate diagonalisation
        method is used.

        >>> f = cofasu.Fasu(topfile, trjfile, filter='resid 1-30')
        >>> c = cofasu.Cofasu(f)
        >>> p1 = Pcz(c)
        >>> ev1 = p1.evals()
        >>> p2 = Pcz(c, fastmethod=True)
        >>> ev2 = p2.evals()
        >>> print(np.allclose(ev1, ev2))
        True

        '''

        self.preloaded = preload
        self.version = version
        self.cofasu = cofasu
        self.quality = quality
        self.rank = rank
        self.natoms = self.cofasu.natoms
        self.nframes = self.cofasu.numframes()

        if quality < 0 or quality > 100:
            raise ValueError('quality must lie in the range 0:100.')

        if rank == 0:
            log.info('Pcz: {0} atoms and {1} snapshots'.format(self.natoms, self.nframes))
        if covar is None:
            if rank == 0:
                log.info('Pcz: least-squares fitting snapshots')
            if target is None:
                time_avg_0 = time()
                self._avg = cofasu.fitted_average()
                time_avg_1 = time()
                if rank == 0:
                    log.info(
                'Pcz: Time for trajectory fitting: {0:.2f} s\n'.format(time_avg_1 - time_avg_0))

            else:
                self._avg = target

            if rank == 0:
                log.info('Pcz: calculating covariance matrix')
            if fastmethod:
                # adapted from Ian Dryden's R code. If you have
                # n atoms and p snapshots, then the conventional
                # way to do the pca is to calculate the [3n,3n]
                # covariance matrix and then diagonalise that.
                # However if p < 3n, then the last 3n-p eigenvectors
                # and values are meaningless anyway, and instead
                # you can calculate the [p,p] covariance matrix,
                # diagonalise that, and then do some
                # data massaging to recover the full eigenvectors
                # and eigenvalues from that. Here we extend the
                # approach to situations where 3n is just too big
                # for diagonalisation in a reasonable amount of
                # time, by just taking a selection of snapshots
                # from the full set (<< 3n) and applying this
                # approach. Obviously this is an approximate
                # method, but it may be good enough.
                if rank == 0:
                    log.info("Using fast approximate diagonalisation method")
                nsamples = min(100,self.nframes)
                tmptrj = np.zeros((3*self.natoms,nsamples))
                stepsize = self.nframes/nsamples
                j = 0
                for i in range(nsamples):
                    x2 = self.cofasu.coords(j,target=self._avg) - self._avg
                    tmptrj[:,i] = x2.flatten()
                    j += stepsize
                cv = np.dot(tmptrj.T,tmptrj.conj())/nsamples
            else:
                cv = self.cofasu.cov(self._avg, preload=preload)
        else:
            # Covariance matrix supplied. This requires that the corresponding
            # target structure is given as well.
            if target is None:
                raise ValueError('A defined covariance matrix requires a defined target.')
            else:
                self._avg = target
                cv = covar

        if rank == 0:
            log.info('Pcz: diagonalizing covariance matrix')
        time_diag_cov_0 = time()
        w, v = eigh(cv)
        if fastmethod:
            vv = np.zeros(nsamples)
            z = np.dot(tmptrj,v)
            for i in range(nsamples):
                vv[i] = np.sqrt((z[:,i]*z[:,i]).sum())
                z[:,i] = z[:,i]/vv[i]

            w2 = np.sqrt(abs(w/nsamples))*vv
            w = w2[w2.argsort()]
            v = z[:,w2.argsort()]

        cs = np.cumsum(w[::-1])
        self.totvar = cs[-1]
        tval = cs[-1] * self.quality / 100
        i = 0
        while cs[i] < tval:
            i += 1

        i += 1
        self.nvecs = i
        # override this with req_evecs, if given:
        if req_evecs is not None:
            if req_evecs > len(w):
                if rank == 0:
                    log.error('Pcz: you asked for {0} eigenvectors but there are only {1} available.'.format(req_evecs, len(w)))
            else:
                self.nvecs = req_evecs
                i = req_evecs

        self._evals = w[-1:-(i + 1):-1]
        self._evecs = v[:, -1:-(i + 1):-1].T
        time_diag_cov_1 = time()
        if rank == 0:
            log.info(
                'Pcz: Time for diagonalizing covariance matrix: {0:.2f} s\n'.format(time_diag_cov_1 - time_diag_cov_0))

        if rank == 0:
            log.info('Pcz: calculating projections')
        time_proj_calc_0 = time()
        if self.preloaded:
            trj = np.zeros((self.nframes, self.natoms * 3))
            for i in range(self.nframes):
                mob = self.cofasu.coords(i)
                mob -= mob.mean(axis=0)
                R, rms = rotation_matrix(self._avg, mob)
                trj[i, :] = np.dot(mob, R).flatten() - self._avg.flatten()

            self._projs = np.dot(trj, self._evecs.T).T
        else:
            if rank == 0:
                log.info("Calculating projections big files method")
            max_ts = 1000
            num_blocks = self.nframes/max_ts
            ts_extra = self.nframes%max_ts
            self._projs = np.zeros((self.nvecs, self.nframes))
            trj = np.zeros((max_ts, self.natoms * 3))
            for i in range(num_blocks):
                for j in range(max_ts):
                    index = i*max_ts + j
                    mob = self.cofasu.coords(index)
                    mob -= mob.mean(axis=0)
                    R, rms = rotation_matrix(self._avg, mob)
                    trj[j, :] = np.dot(mob, R).flatten() - self._avg.flatten()

                self._projs[:,i*max_ts:(i+1)*max_ts] = np.dot(trj, self._evecs.T).T
            trj = np.zeros((ts_extra, self.natoms * 3))
            for j in range(ts_extra):
                index = num_blocks*max_ts + j
                mob = self.cofasu.coords(index)
                mob -= mob.mean(axis=0)
                R, rms = rotation_matrix(self._avg, mob)
                trj[j, :] = np.dot(mob, R).flatten() - self._avg.flatten()

            self._projs[:,num_blocks*max_ts:] = np.dot(trj, self._evecs.T).T
        time_proj_calc_1 = time()
        if rank == 0:
            log.info('Pcz: Time for calculating projections: {0:.2f} s\n'.format(time_proj_calc_1 - time_proj_calc_0))

    def numframes(self):
        """
        Method to match the cofasu equivalent
        """
        return self.nframes

    def avg(self):
        """
        Returns the average structure contained in the pcz file
        as an (natoms,3) numpy array.
        >>> topfile = "../../../../test/2ozq.pdb"
        >>> trjfile = "../../../../test/2ozq.dcd"
        >>> from MDPlus.core import cofasu
        >>> f = cofasu.Fasu(topfile, trjfile, filter='name CA')
        >>> c = cofasu.Cofasu(f)
        >>> p = Pcz(c)
        >>> print(np.allclose(p.avg()[0], np.array([ 31.323149, 61.575380, 40.136298]), atol=0.001, rtol=0.001))
        True

        """
        return self._avg

    def eval(self, ival):
        """
        Returns an eigenvalue from the file.
        >>> topfile = "../../../../test/2ozq.pdb"
        >>> trjfile = "../../../../test/2ozq.dcd"
        >>> from MDPlus.core import cofasu
        >>> f = cofasu.Fasu(topfile, trjfile, filter='name CA')
        >>> c = cofasu.Cofasu(f)
        >>> p = Pcz(c)
        >>> print(np.allclose(p.eval(4), np.array([2.6183940]), rtol=0.001, atol=0.001))
        True

        """
        if ival >= self.nvecs:
            print 'Error - only ', self.nvecs, ' eigenvectors present'
            return 0.0
        else:
            return self._evals[ival]

    def evals(self):
        """
        Returns an array of all eigenvalues in the file.
        >>> topfile = "../../../../test/2ozq.pdb"
        >>> trjfile = "../../../../test/2ozq.dcd"
        >>> from MDPlus.core import cofasu
        >>> f = cofasu.Fasu(topfile, trjfile, filter='name CA')
        >>> c = cofasu.Cofasu(f)
        >>> p = Pcz(c)
        >>> print(np.allclose(p.evals()[4], np.array([2.6183940]), rtol=0.001, atol=0.001))
        True

        """
        return self._evals

    def evec(self, ivec):
        """
        Returns a chosen eigenvector from the file in the
        form of a (3*natoms) numpy array.
        >>> topfile = "../../../../test/2ozq.pdb"
        >>> trjfile = "../../../../test/2ozq.dcd"
        >>> from MDPlus.core import cofasu
        >>> f = cofasu.Fasu(topfile, trjfile, filter='name CA')
        >>> c = cofasu.Cofasu(f)
        >>> p = Pcz(c)
        >>> print(np.allclose(abs(p.evec(1)[12]), np.array([0.00865751377]), rtol=0.001, atol=0.001))
        True

        """
        if ivec >= self.nvecs:
            print 'Error - only ', self.nvecs, 'eigenvectors present'
            return None
        else:
            return self._evecs[ivec, :]

    def evecs(self):
        """
        Returns all eigenvectors in the file in the form of a
        (nvecs,3*natoms) numpy array.
        >>> topfile = "../../../../test/2ozq.pdb"
        >>> trjfile = "../../../../test/2ozq.dcd"
        >>> from MDPlus.core import cofasu
        >>> f = cofasu.Fasu(topfile, trjfile, filter='name CA')
        >>> c = cofasu.Cofasu(f)
        >>> p = Pcz(c)
        >>> e = p.evecs()
        >>> print(e.shape)
        (18, 471)

        >>> element = abs(e[1,12])
        >>> print(np.allclose(element, np.array([0.0086575138]), rtol=0.001, atol=0.001))
        True

        """
        return self._evecs

    def proj(self, iproj):
        """
        Returns an array of the projections along a given eigenvector. There
        will be one value per snapshot.
        >>> topfile = "../../../../test/2ozq.pdb"
        >>> trjfile = "../../../../test/2ozq.dcd"
        >>> from MDPlus.core import cofasu
        >>> f = cofasu.Fasu(topfile, trjfile, filter='name CA')
        >>> c = cofasu.Cofasu(f)
        >>> p = Pcz(c)
        >>> prj = abs(p.proj(3))
        >>> print(np.allclose(prj[21], 0.33430696, rtol=0.001, atol=0.001))
        True

        """
        if iproj >= self.nvecs:
            print 'Error - only ', self.nvecs, 'eigenvectors present'
            return None
        else:
            return self._projs[iproj, :]

    def scores(self, framenumber):
        """
        Method that returns the scores (projections) corresponding to
        a chosen snapshot (zero-based).
        >>> topfile = "../../../../test/2ozq.pdb"
        >>> trjfile = "../../../../test/2ozq.dcd"
        >>> from MDPlus.core import cofasu
        >>> f = cofasu.Fasu(topfile, trjfile, filter='name CA')
        >>> c = cofasu.Cofasu(f)
        >>> p = Pcz(c)
        >>> s = abs(p.scores(12))
        >>> print(np.allclose(s[3], 0.58309314, rtol=0.001, atol=0.001))
        True

        """
        if( framenumber >= self.nframes):
             return None
        else:
             x = np.zeros(self.nvecs)
             for i in range(self.nvecs):
                 x[i] = self.proj(i)[framenumber]
             return x

    def coords(self, framenumber):
        """
        Synonym for frame() method, to match cofasu.
        """
        return self.frame(framenumber)

    def frame(self,framenumber):
        """
        Method to return the coordinates of the given frame - i.e.
        to decompress a snapshot. The data is returned as a (natoms,3) 
        numpy array.
        >>> topfile = "../../../../test/2ozq.pdb"
        >>> trjfile = "../../../../test/2ozq.dcd"
        >>> from MDPlus.core import cofasu
        >>> f = cofasu.Fasu(topfile, trjfile, filter='name CA')
        >>> c = cofasu.Cofasu(f)
        >>> p = Pcz(c, quality=95)
        >>> ref = c.coords(5)
        >>> x = p.frame(5)
        >>> cofasu.rmsd(ref, x) < 0.19
        True

        """
        if(framenumber >= self.nframes):
            return None
        else:
            scores = self.scores(framenumber)
            return self.unmap(scores)


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

    def unmap(self,scores):
        """
        Method to return the coordinates corresponding to a given
        set of scores. If the scores vector has less than nvec elements,
        Pczfile.expand is called to fill in the missing values.
        >>> topfile = "../../../../test/2ozq.pdb"
        >>> trjfile = "../../../../test/2ozq.dcd"
        >>> from MDPlus.core import cofasu
        >>> f = cofasu.Fasu(topfile, trjfile, filter='name CA')
        >>> c = cofasu.Cofasu(f)
        >>> p = Pcz(c)
        >>> a = p.avg()
        >>> a2 = p.unmap(np.zeros(p.nvecs))
        >>> print(cofasu.rmsd(a, a2) < 0.001)
        True

        """
        x = self.avg()
        if len(scores) < self.nvecs:
            scores = self.expand(scores)
        for i in range(self.nvecs):
            x = x + (self.evec(i)*scores[i]).reshape((self.natoms,3))
        return x

    def expand(self, scores):
        """
        Method to complete a truncated list of scores with values for
        missing eigenvectors. Basically the pcz file is scanned to
        find the frame with the best match to the scores given, the
        missing values are then copied directly from this frame.
        >>> topfile = "../../../../test/2ozq.pdb"
        >>> trjfile = "../../../../test/2ozq.dcd"
        >>> from MDPlus.core import cofasu
        >>> f = cofasu.Fasu(topfile, trjfile, filter='name CA')
        >>> c = cofasu.Cofasu(f)
        >>> p = Pcz(c)
        >>> e = abs(p.expand([1.0, 2.0, 3.0]))
        >>> print(np.allclose(e[3:6], [2.5556779, 0.9558872, 2.2131005], rtol=0.001, atol=0.001))
        True

        """
        nvecs = self.nvecs
        ns = len(scores)
        if ns >= nvecs:
            return scores
        else:
            newscores = np.zeros(nvecs)
            newscores[0:ns] = scores
            temp = self._projs
            best = 0
            err = ((temp[0:ns,0]-scores)*(temp[0:ns,0]-scores)).sum()
            newscores[ns:nvecs] = temp[ns:nvecs,0]
            for frame in range(self.nframes):
                newerr = ((temp[0:ns,frame]-scores)*(temp[0:ns,frame]-scores)).sum()
                if newerr < err:
                    err = newerr
                    newscores[ns:nvecs] = temp[ns:nvecs,frame]
            return newscores

    def map(self,crds):
        """
        Method to map an arbitrary coordinate set onto the PC model. The
        coordinate set should be a (natom,3) array-like object that matches
        (for size) what's in the pczfile. An array of projections will be 
        returned, one value for each eignevector in the pcz file.
        >>> topfile = "../../../../test/2ozq.pdb"
        >>> trjfile = "../../../../test/2ozq.dcd"
        >>> from MDPlus.core import cofasu
        >>> f = cofasu.Fasu(topfile, trjfile, filter='name CA')
        >>> c = cofasu.Cofasu(f)
        >>> p = Pcz(c)
        >>> crds = (c.coords(10) + c.coords(20))* 0.5
        >>> print(np.allclose(abs(p.map(crds)[:3]),[0.3844803, 1.6537842, 2.29433996], rtol=0.001, atol=0.001))
        True

        """
        crdset=np.array(crds)
        if np.shape(crdset) == (self.natoms, 3):
          avg = np.reshape(self.avg(),(self.natoms,3))
          rot = np.zeros(9, dtype=np.float64)
          crdset_cog = crdset.sum(axis=0)/crdset.shape[0]
          avg_cog = avg.sum(axis=0)/avg.shape[0]
          c = crdset.copy()
          c -= crdset_cog
          rmsd = qcp.CalcRMSDRotationalMatrix((avg-avg_cog).T.astype(np.float64),
                                    c.T.astype(np.float64),
                                    crdset.shape[0], rot, None)
          R = np.matrix(rot.reshape(3,3))
          c = c * R
          c += avg_cog
          c -= avg
          prj = np.zeros(self.nvecs)
          for i in range(self.nvecs):
              prj[i]=(np.dot(c.flatten(),self.evec(i)))
          return prj
        else:
          return None
  

    def write(self, filename,  title='Created by pcz.write()'):
        """
        Write out the PCZ file. At the moment only the PCZ4 and PCZ6 formats
        are  implemented.
        """
        if self.version == 'PCZ7' and not h5py_available:
            log.info("WARNING: The PCZ6 format will be used because the h5py module required for PCZ7 is not available. Please install the PCZ7 extra requirements to be able to use it. The command to do this is: pip install pyPcazip[PCZ7]")
            self.version = 'PCZ6'

        if self.version == 'UNKN':
            if h5py_available:
                self.version = 'PCZ7'
            else:
                self.version = 'PCZ6'

        if self.version != 'PCZ4' and self.version != 'PCZ6' and self.version != 'PCZ7':
            raise TypeError('Only PCZ4/6/7 formats supported')

        if self.rank == 0:
            log.info("Using "+self.version+" format")

        if self.version == 'PCZ4' or self.version == 'PCZ6':
            f = open(filename, 'wb')
            f.write(struct.pack('4s80s3if', self.version, title, self.natoms, self.nframes, self.nvecs, self.totvar))
            f.write(struct.pack('4i', 0, 0, 0, 0))
            for v in self.avg().flatten():
                f.write(struct.pack('f', v))
            for i in range(self.nvecs):
                for v in self.evec(i):
                    f.write(struct.pack('f', v))
                f.write(struct.pack('f', self.eval(i)))

                # All in memory
                if self.preloaded:
                    projection = self.proj(i)

                # loop file (low memory demand)
                else:
                    projection = np.zeros(self.nframes)
                    for t in range(self.nframes):

                        mob = self.cofasu.coords(t)
                        mob -= mob.mean(axis=0)
                        R, rms = rotation_matrix(self._avg, mob)

                        time_step_coords = np.dot(mob, R).flatten() - self._avg.flatten()
                        projection[t] = np.dot(time_step_coords, self._evecs[i].T).T

                if self.version == 'PCZ4':
                    for v in projection:
                        f.write(struct.pack('f', v))
                elif self.version == 'PCZ6':
                    pinc = (projection.max() - projection.min()) / 65534
                    p0 = (projection.max() + projection.min()) / 2
                    f.write(struct.pack('2f', p0, pinc))
                    for v in projection:
                        f.write(struct.pack('h', np.int16((v - p0) / pinc)))

                else:
                    print 'Error - only PCZ4 and PCZ6 formats supported'

            f.close()
            return
        elif self.version == 'PCZ7':
            if not h5py_available:
                print "Error: h5py module not available. Please install the PCZ7 extra requirements to be able to use it. The command to do this is: pip install pyPcazip[PCZ7]"
                sys.exit(0)
            f = h5py.File(filename, "w")
            # Write evecs and evalues
            evec_array = []
            eval_array = []
            for evec_index in xrange(self.nvecs):
                evec_array.append([])
                for v in self.evec(evec_index):
                    evec_array[evec_index].append(v)
                eval_array.append(self.eval(evec_index))
            f.create_dataset("evec_dataset", (self.nvecs, (3 * self.natoms)), dtype='f', data=np.array(evec_array))
            f.create_dataset("eval_dataset", (self.nvecs,), dtype='f', data=np.array(eval_array))
            # Write reference coordinates
            f.create_dataset("ref_coord_dataset", (len(self.avg().flatten()),), dtype='f',
                             data=np.array(self.avg().flatten()))
            # Write properties
            f.attrs['version'] = self.version
            f.attrs['title'] = title
            f.attrs['natoms'] = self.natoms
            f.attrs['nframes'] = self.nframes
            f.attrs['nvecs'] = self.nvecs
            f.attrs['quality'] = self.totvar

            # Loop on every ts
            proj_dataset = f.create_dataset("proj_dataset", (self.nframes,
            self.nvecs), dtype='int16')
            p0_dataset = f.create_dataset("p0_dataset", (self.nframes,), dtype='f')
            pinc_dataset = f.create_dataset("pinc_dataset", (self.nframes,), dtype='f')



            for ts_index in xrange(self.nframes):
                projection_values = []
                for evec_index in xrange(self.nvecs):
                    # Prepare coords of ts
                    if self.preloaded:
                        projection_values.append(self.proj(evec_index)[ts_index])
                    else:
                        mob = self.cofasu.coords(ts_index)
                        mob -= mob.mean(axis=0)
                        R, rms = rotation_matrix(self._avg, mob)
                        time_step_coords = np.dot(mob, R).flatten() - self._avg.flatten()
                        projection_values.append((np.dot(time_step_coords, self._evecs[evec_index].T).T).item(0))

                # Write coords of ts
                projection_values = np.array(projection_values)
                pinc = (projection_values.max() - projection_values.min()) / 65534
                if pinc == 0:
                    pinc == numpy.nextafter(0, 1)
                p0 = (projection_values.min() + projection_values.max()) / 2
                p0_dataset[ts_index] = p0
                pinc_dataset[ts_index] = pinc
                projection_values = projection_values - p0
                proj_dataset[ts_index] = (projection_values / pinc).astype(np.int16)


        else:
            raise TypeError('Only PCZ4/6/7 formats supported')

if __name__ == "__main__":
    import __builtin__
    __builtin__.parallelbuiltin = False
    import doctest
    doctest.testmod()
