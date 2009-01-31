# $Id$
"""The formats submodule contains code to read coordinates, either
single frames (e.g. the PDB module) or trajectories (such as the DCD
reader). All readers are supposed to expose a Reader object that
represents a common API to the other code.

The Universe contains the API entry point attribute

  Universe.trajectory

that points to the actual Reader object; all Readers are supposed to
be accessible through this entry point in the same manner ('If it
quacks and walks like a duck it probably is a duck as far as I am
concerned.')
"""
# TODO: define the API here (for the moment, look at DCD.DCDReader)

__all__ = ['DCD', 'PDB', 'XTC']

import PDB,DCD

# frame writers: export to single frame formats such as PDB, gro, crd
# Signature:
#   W = FrameWriter(filename)
#   W.write(AtomGroup)
_frame_writers = {'pdb': PDB.PrimitivePDBWriter,
                 }

# trajectory writers: export frames, typically only saving coordinates
# (although PDB movies are the exception); maybe look at OpenBabel as
# not to reinvent the wheel.
# Signature:
#   W = TrajectoryWriter(filename,numatoms,**kwargs)
#   W.write_next_timestep(TimeStep)
_trajectory_writers = {'dcd': DCD.DCDWriter,
                       }
