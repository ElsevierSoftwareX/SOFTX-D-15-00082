# The mapping module. Creates and does various useful things with
# multidimensional histograms.

import numpy as np
import scipy.ndimage as nd

class Map():
    
    def __init__(self, arr, resolution=10, boundary=0, limits=None):
        '''
        Initialise a distribution map. This is essentially a multidimensional 
        histogram produced from a set of data, arr(N,D) where N is the
        number of points and D is the number of dimensions. "resolution" 
        specifies the number of bins in each dimension, and may be a single 
        number or a list of length D. The optional limits parameter sets the
        histogram boundaries, which are otherwise set automatically to include
        all the data.
        '''
        self.ndim = arr.shape[1]
        self.boundary = boundary
        # resolution may be a single value or a list of values, 1 per dimension
        self.resolution = np.zeros(self.ndim)
        self.resolution[:] = resolution
        if limits is None:
            self.limits = []
            min = arr.min(axis=0)
            max = arr.max(axis=0)
            # set a boundary
            for i in range(self.ndim):
                buff = (max[i]-min[i])/(self.resolution[i]-2*boundary)*boundary*1.01
                self.limits.append((min[i]-buff,max[i]+buff))
        else:
            self.limits=limits

        # Create the histogram
        self._H, self._edges = np.histogramdd(arr, bins=self.resolution, range=self.limits)
        self.shape = self._H.shape
        # find the bin dimensions (this should be buff, but don't assume)
        self.cellsize = []
        for i in range(self.ndim):
            self.cellsize.append(self._edges[i][1]-self._edges[i][0])
        # calculate the bin volume, and number of sampled bins
        self.cellvol = self.cellsize[0]
        for l in self.cellsize[1:]:
            self.cellvol = self.cellvol*l
        self.coverage = self._H[np.where(self._H > 0)].size
        # correct coverage if there is a boundary:
        v0 = 1
        v1 = 1
        for d in self.resolution:
           v0 = v0*d
           v1 = v1*(d-2*self.boundary)
        self.deadspace = v0 - v1
        self.volume = v0
        # now give preliminary values to the cluster-related data.
        self.ID = self._H.copy()
        indMax = np.unravel_index(self.ID.argmax(),self.ID.shape)
        self.ID = np.where(self._H > 0, 1, 0)
        self.ID[indMax] = -1
        self.sizes = [0, arr.shape[0]]

    def map(self, vec):
        '''
        Returns the index of the bin that corresponds to point vec.
        '''
        if len(vec) != self.ndim:
            raise ValueError('Error - vector has wrong length.')
        indx = []
        i = 0
        for c in vec:
            indx.append(np.digitize((c,c),self._edges[i])[0]-1)
            i += 1
        return indx

    def unmap(self, indx):
        '''
        Returns a vector corresponding to the mid point of the
        bin with the given index.
        '''
        if len(indx) != self.ndim:
            raise ValueError('Error - index has wrong number of dimensions.')
        vec = []
        i = 0
        for c in indx:
            vec.append(self.limits[i][0]+c*self.cellsize[i]+ 0.5*self.cellsize[i])
            i += 1
        return vec


    def cluster_id(self, vec):
        '''
        Returns the cluster ID of the point vec.
        '''
        if len(vec) != self.ndim:
            raise ValueError('Error - vector has wrong length.')
        indx = []
        i = 0
        for c in vec:
            indx.append(np.digitize((c,c),self._edges[i])[0]-1)
            i += 1
        return self.ID[tuple(indx)]

    def cluster_centre(self, id):
        '''
        Returns a vector corresponding to the mid point of the
        cluster with the given id. It is assumed that these are
        the bins with negative cluster ids.
        '''
        if id < 0:
            id = -id
        w = np.where(self.ID == -id)
        e = []
        for i in range(self.ndim):
            e.append(self._edges[i][w[i][0]]+0.5*self.cellsize[i])
        return e 

    def cluster_size(self, id):
        '''
        Returns the number of samples in cluster = id
        '''
        return self.sizes[id-1]

def coco(map, method='coco', npoints=1):
    '''
    The CoCo (Complementary Coordinates) methods. The input is an
    N-dimensional histogram defining the occupancy of the space.
    Various CoCo methods will identify 'interesting' regions to
    be sampled next.
    '''
    if method == 'coco':
        '''
        returns new points, generated using the COCO procedure,
        in the form of an (npoints,D) numpy array, where D is the number of
        dimensions in the map.
        '''
        cp = np.zeros((npoints,map.ndim))
        # make a temporary binary image, and invert
        tmpimg = np.where(map._H > 0, 0, 1)
        for i in range(npoints):
            dis = nd.morphology.distance_transform_edt(tmpimg)
            indMax = np.unravel_index(dis.argmax(),dis.shape)
            for j in range(map.ndim):
                cp[i,j]=map._edges[j][0]+indMax[j]*map.cellsize[j]
            
            tmpimg[indMax] = 0
        return cp

    elif method == 'hpoints':
        '''
        hpoints returns new points that form a halo of unsampled space
        just beyond the sampled region.
        '''
        # This is the halo filter:
        def f(arr):
            cval = arr[len(arr)/2]
            if cval == 0 and np.max(arr) > 0:
                return 1
            else:
                return 0

        halo = nd.filters.generic_filter(map._H,f,size=3,mode='constant')
        npoints = int(np.sum(halo))
        hp = np.zeros((npoints,map.ndim))
        for i in range(npoints):
            indMax = np.unravel_index(halo.argmax(),map.shape)
            for j in range(map.ndim):
                hp[i,j]=map.edges[j][0]+indMax[j]*map.cellsize[j]
            
            halo[indMax] = 0
        return hp

    elif method == 'fpoints':
        '''
        fpoints returns new points at the frontier of sampled space
        '''
        # This is the frontier filter:
        def f(arr):
            cval = arr[len(arr)/2]
            if cval > 0 and np.min(arr) == 0:
                return 1
            else:
                return 0

        front = nd.filters.generic_filter(map._H,f,size=3,mode='constant')
        npoints = int(np.sum(front))
        fp = np.zeros((npoints,map.ndim))
        for i in range(npoints):
            indMax = np.unravel_index(front.argmax(),map.shape)
            for j in range(map.ndim):
                fp[i,j]=map._edges[j][0]+indMax[j]*map.cellsize[j]
            
            front[indMax] = 0
        return fp

    elif method == 'bpoints':
        '''
        bpoints() returns new points not at the frontier of sampled space
        '''
        # This is the buried filter:
        def f(arr):
            cval = arr[len(arr)/2]
            if cval > 0 and np.min(arr) > 0:
                return 1
            else:
                return 0

        bur = nd.filters.generic_filter(map._H,f,size=3,mode='constant')
        npoints = int(np.sum(bur))
        bp = np.zeros((npoints,map.ndim))
        for i in range(npoints):
            indMax = np.unravel_index(bur.argmax(),map.shape)
            for j in range(map.ndim):
                bp[i,j]=map._edges[j][0]+indMax[j]*map.cellsize[j]
            
            bur[indMax] = 0
        return bp

    elif method == 'rpoints':
        '''
        rpoints() returns one point per bin of sampled space, and its weight
        '''

        tmpimg = map._H.copy()
        hsum = np.sum(map._H)
        npoints = tmpimg[np.where(tmpimg > 0)].size
        wt = np.zeros((npoints))
        rp = np.zeros((npoints,map.ndim))
        for i in range(npoints):
            indMax = np.unravel_index(tmpimg.argmax(),map.shape)
            for j in range(map.ndim):
                rp[i,j]=map._edges[j][0]+indMax[j]*map.cellsize[j]
            
            tmpimg[indMax] = 0
            wt[i] = map._H[indMax]/hsum
        return rp,wt

    else:
        raise ValueError('Unknown method: {}'.format(method))




def watershed(map):
    '''
    Watershed clustering. The input is a map object that defines the
    sampling of the space. The method sets the ID data in the map.
    '''
    # The method requires a clear 1-bin boundary around the distribution:
    if map.boundary < 1:
        raise ValueError('The input map needs to have a boundary.')

    # define the function that will locate local maxima in the distibution: 
    def f(arr):
        if np.argmax(arr) == len(arr)/2:
            return 1
        else:
            return 0
    # now find them (calling them 'labels' but at the moment they are
    # really 'maxima', will become labels properly next):
    labels = nd.filters.generic_filter(map._H,f,size=3,mode='constant')
    # now turn them into proper labels:
    nd.measurements.label(labels,output=labels)
    # The array used for doing the watershedding needs to simultaneously
    # encode both the number of counts in the bin, and the bin label. We
    # do this by raising the bin counts by enough powers of ten to leave
    # room in the low order bits for the label.
    maxval = np.max(labels)
    scale = int(10**(np.ceil(np.log10(maxval+1))))
    
    # Now put all this into the ID array:
    map.ID = map.ID * scale + labels
    # Define the function that assigns bins to labels. For each currently
    # unlabelled bin, we search around to find the neighbour bin with
    # the highest occupancy. If this bin is already labelled, the current
    # bin gets this label too.
    def f2(arr,scale):
        cval = arr[len(arr)/2]
        if cval > 0 and cval%scale == 0:
            maxval = np.max(arr)
            if maxval%scale > 0:
                cval = cval +maxval%scale
        return cval

    # now we can do the watershed-like clustering. The process is iterative,
    # labels "spread out" in a downhill direction from the maxima until all
    # bins are labelled.
    IDnew = nd.filters.generic_filter(map.ID,f2,size=3,mode='constant',extra_arguments=(scale,))
    while np.any(IDnew != map.ID):
        map.ID = IDnew
        IDnew = nd.filters.generic_filter(map.ID,f2,size=3,mode='constant',extra_arguments=(scale,))

    # make ID nice for indexing:
    map.ID = map.ID.astype(np.int)
    # remove the bincount data, leaving just the labels:
    map.ID = map.ID%scale
    # Flatten for now to make the next few steps a bit easier:
    map.ID = map.ID.flatten()
    # now we find out big each cluster is:
    sizes = np.zeros(1+int(np.max(map.ID)), dtype=np.int)
    for i in range(len(map.ID)):
        sizes[int(map.ID[i])] += map._H.flatten()[i]
    # now reassign labels so largest group has label=1, etc.

    newlab = np.argsort(np.argsort(sizes))
    newlab = abs(newlab-newlab.max()) + 1
    map.sizes = np.sort(sizes)[::-1]

    for i in range(len(map.ID)):
        map.ID[i] = newlab[map.ID[i]]

    # put ID back into the right shape:
    map.ID = map.ID.reshape(map._H.shape)
    # now mark "root" structures with negative indices:
    map.ID = np.where(labels>0, -map.ID, map.ID)
    # set empty bins to have label=0:
    map.ID = np.where(map.ID==maxval+1, 0, map.ID)
