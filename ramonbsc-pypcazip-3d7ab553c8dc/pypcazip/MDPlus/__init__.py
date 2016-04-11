import sys
import os.path as op
import traceback

# The next line is just to make sure Deprecation warnings stay suppressed:
import MDAnalysis.analysis.align as _maa_
from MDAnalysis import Universe as MDAUniverse
from MDAnalysis import Writer
import MDAnalysis

from MDPlus.coordinates import NUMPY

try:
    import mdtraj as md
    mdtraj_available = True
except ImportError:
    mdtraj_available = False

#Register the new readers: 
# NPZ: numpy matrix binary file format, encoding coordinates.
# PCZ: PCAzip & PCAunzip
# NUMPY: (object), for various uses, among them support MDtraj file formats (Binpos, etc)
# TRJ: Upgraded version of MDAnalysis TRJ Reader/Writer (fixes bug in Reader, 
#      adds Writer).

#MDAnalysis.coordinates._trajectory_readers['NPZ'] = NPZ.NPZReader
#MDAnalysis.coordinates._trajectory_readers_permissive['NPZ'] = NPZ.NPZReader
#MDAnalysis.coordinates._trajectory_writers['NPZ'] = NPZ.NPZWriter

#MDAnalysis.coordinates._trajectory_readers['PCZ'] = PCZ.PCZReader
#MDAnalysis.coordinates._trajectory_readers_permissive['PCZ'] = PCZ.PCZReader
#MDAnalysis.coordinates._trajectory_writers['PCZ'] = PCZ.PCZWriter

MDAnalysis.coordinates._trajectory_readers['NUMPY'] = NUMPY.NUMPYReader
#MDAnalysis.coordinates._trajectory_readers_permissive['NUMPY'] = NUMPY.NUMPYReader
#MDAnalysis.coordinates._trajectory_readers['MDCRD'] = TRJ.TRJReader
#MDAnalysis.coordinates._trajectory_readers_permissive['MDCRD'] = TRJ.TRJReader
#MDAnalysis.coordinates._trajectory_writers['MDCRD'] = TRJ.TRJWriter
#MDAnalysis.coordinates._trajectory_readers['TRJ'] = TRJ.TRJReader
#MDAnalysis.coordinates._trajectory_readers_permissive['TRJ'] = TRJ.TRJReader
#MDAnalysis.coordinates._trajectory_writers['TRJ'] = TRJ.TRJWriter

# We need to hack the Universe
class Universe(MDAUniverse):

    def __init__(self, topologyfile, coordinatefile=None, **kwargs):

        if not coordinatefile:
            super(Universe,self).__init__(topologyfile, coordinatefile=None, **kwargs)
        
        else:
            format = kwargs.pop('format', None)
            if 'binpos' in op.splitext(coordinatefile)[1].lower():
                format = 'binpos'
            elif 'netcdf' in op.splitext(coordinatefile)[1].lower():
            # Extension *.netcdf may not be automatically recognised by MDAnalysis as "ncdf".
                format = 'ncdf'
            elif format is None:
                format = op.splitext(coordinatefile)[1][1:].lower()
                                   
            try:
                super(Universe, self).__init__(topologyfile, coordinatefile=coordinatefile, format=format, **kwargs)
                
            except ValueError, err:
                if 'binpos' in format:
                    if not mdtraj_available:
                        print "Error: mdtraj module not available. Please, install the BINPOS extra requirements to be able to use this format. The command to do this requirements is: pip install pyPcazip[BINPOS]"
                        sys.exit(0)
                    # Obtain an mdtraj trajectory object from mdtraj.
                    t = md.load_binpos(coordinatefile,topologyfile)
                    # Get coordinates in shape=(n_frames, n_atoms, 3).
                    # MDtraj works in Nanometer Units!! Need to rescale to Angstrom (x10).
                    trajNp = t.xyz * 10 
                    # Re-shape to (n_frames, 3 x n_atoms) where every frame is (a1_x, a1_y, a1_z, a2_x, ... aN_z).
                    trajNp = trajNp.reshape(t.n_frames,-1)
         
                
                    kwargs = {'nparray':trajNp, 'permissive':True} 
                    super(Universe,self).__init__(topologyfile,'x.numpy', **kwargs)
    
                
                else:
                    raise TypeError("Universe.load_new() cannot find an appropriate coordinate reader "
                                "for file %r.\n%r" % (coordinatefile, err))                

