# $Id$
"""The formats sub module contains code to read coordinates, either
single frames (e.g. the PDB module) or trajectories (such as the DCD
reader). All readers are supposed to deliver a Reader object that
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
_frame_writers = {'pdb': PDB.PrimitivePDBWriter,
                 }
_trajectory_writers = {'dcd': DCD.DCDWriter,
                       }
