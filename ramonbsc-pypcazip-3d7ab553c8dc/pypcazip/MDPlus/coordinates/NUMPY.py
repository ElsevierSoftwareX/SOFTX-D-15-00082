import MDAnalysis
from MDAnalysis.coordinates.base import Timestep
import numpy as np

class NUMPYReader(MDAnalysis.coordinates.base.Reader):
    """
    Reads coordinates from a (num_frames, num_atoms x 3) numpy array
    into MDAnalysis universes. Modified from the XYZ reader in the MDAnalysis distribution
    As it is a numpy array, all in memory.
    """
    format = "Numpy"
    units = {'time': 'ns', 'length': 'Angstrom'}

    def __init__(self, filename='inmemory.numpy', **kwargs):

        self.nparray = kwargs['nparray']
        self.nparray = self.nparray.reshape(self.nparray.shape[0], -1)
        
        self.fixed = 0
        self.skip = 1
        self.periodic = False
        self.skip_timestep = 0
        self.delta = 0
        self.n_atoms = self.nparray[0].shape[0] / 3
        self.n_frames = self.nparray.shape[0]
        self.ts = Timestep(self.n_atoms)
        self.ts.frame = -1
        # Read in the first timestep
        self._read_next_timestep()

    def __iter__(self):
        self.ts.frame = -1 
        while True:
            try:
                yield self._read_next_timestep()
            except EOFError:
                self.close()
                raise StopIteration

    def _read_next_timestep(self, ts=None):
        # check that the timestep object exists
        if ts is None:
            ts = self.ts

        if ts.frame < self.n_frames - 1:
            # stop when the cursor has reached the end of that block
            ts.frame += 1
            ts._unitcell = np.zeros(6, np.float32)
            ts._x[:] = self.nparray[ts.frame][0::3]
            ts._y[:] = self.nparray[ts.frame][1::3]
            ts._z[:] = self.nparray[ts.frame][2::3]
            return ts
        raise EOFError

    def _jump_to_frame(self,frame):
        ts = self.ts
        ts._unitcell = np.zeros(6, np.float32)
        ts._x[:] = self.nparray[frame][0::3]
        ts._y[:] = self.nparray[frame][1::3]
        ts._z[:] = self.nparray[frame][2::3]
        ts.frame = frame
        return ts

    
    def populate_ts(self, ts, frame):
        ts._unitcell = np.zeros(6, np.float32)
        ts._x[:] = self.nparray[frame][0::3]
        ts._y[:] = self.nparray[frame][1::3]
        ts._z[:] = self.nparray[frame][2::3]
        ts.frame = frame 

    def rewind(self):
        """reposition on first frame"""
        # reset ts
        ts = self.ts
        ts.status = 1
        ts.frame = -1
        ts.step = 0
        ts.time = 0
        self.skip_timestep = 0
        # the next method is inherited from the Reader Class and calls _read_next_timestep
        self.next()

    def Writer(self, filename, numatoms=None, **kwargs):
        #NB. Just copied from DCD.py
        #numatoms = kwargs.pop('numatoms', self.numatoms)
        kwargs.setdefault('start', 0)
        kwargs.setdefault('step', self.skip)
        kwargs.setdefault('delta', 10)
        kwargs.setdefault('remarks', 'All a bit too much hardcoded')
        from MDAnalysis.coordinates.DCD import DCDWriter
        return DCDWriter(filename, numatoms, **kwargs)

    def __getitem__(self, frame):
        if (np.dtype(type(frame)) != np.dtype(int)) and (type(frame) != slice):
            raise TypeError
        if np.dtype(type(frame)) == np.dtype(int):
            if frame < 0:
                # Interpret similar to a sequence
                frame += len(self)
            if (frame < 0) or (frame >= len(self)):
                raise IndexError            
            return self._jump_to_frame(frame)
        elif type(frame) == slice: # if frame is a slice object
            if not ((type(frame.start) == int or frame.start is None) and
                    (type(frame.stop) == int or frame.stop is None) and
                    (type(frame.step) == int or frame.step is None)):
                raise TypeError("Slice indices are not integers")
            def iterDCD(start=frame.start, stop=frame.stop, step=frame.step):
                start, stop, step = self._check_slice_indices(start, stop, step)
                for i in xrange(start, stop, step):
                    yield self[i]
            return iterDCD()
