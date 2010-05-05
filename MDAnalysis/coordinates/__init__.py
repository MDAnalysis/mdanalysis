# $Id$
"""
:mod:`MDAnalysis.coordinates` --- reading and writing coordinates
=================================================================

The coordinates submodule contains code to read coordinates, either
single frames (e.g. the PDB module) or trajectories (such as the DCD
reader). All readers are supposed to expose a Reader object that
represents a common API to the other code.

The Universe contains the API entry point attribute

  :attr:`Universe.trajectory`

that points to the actual Reader object; all Readers are supposed to
be accessible through this entry point in the same manner (`duck typing`_)

.. _duck typing: http://c2.com/cgi/wiki?DuckTyping


Trajectory API
--------------

(Draft, orbeckst 2010-04-30)


Timestep
~~~~~~~~

A Timestep instance holds data for the current frame. It is updated whenever a
new frame of the trajectory is read.  The :class:`DCD.Timestep` class is the
primary implementation example.

Methods
.......

  __init__(arg)
      depending on the type of *arg*, Timestep instances are created in
      different ways:
        int
            empty Timestep for *arg* atoms (allocate arrays etc)
        Timestep
            makes a deep copy of the *arg*
        numpy.ndarray
            update the Timesteps positions array with the contents of *arg*

      Anything else raises an exception; in particular it is not possible to
      create a "empty" Timestep instance.
  __getitem__(atoms)
      position(s) of atoms; can be a slice or numpy array and then returns
      coordinate array
  __len__()
      number of coordinates (atoms) in the frame
  __iter__()
      iterator over all coordinates
  copy()
      deep copy of the instance

Attributes
..........

  numatoms
      number of atoms in the frame
  frame
      current frame number
  dimensions
      system box dimensions (x, y, z, alpha, beta, gamma)


Trajectory Reader
~~~~~~~~~~~~~~~~~

Methods
.......

The :class:`DCD.DCDReader` class is the primary implementation example.

**Mandatory methods**

The following methods must be implemented in a Reader class.

 __init__(filename)
     open *filename*
 __iter__()
     allow iteration from beginning to end::
        for ts in trajectory:
            print ts.frame
 close_trajectory()
     close the file and cease I/O
 __del__()
     ensure that the trajectory is closed
 next()
     advance to next time step or raise :exc:`IOError` when moving
     past the last frame
 rewind()
     reposition to first frame


**Optional methods**

Not all trajectory formats support the following methods, either because the
data are not available or the methods have not been implemented. Code should
deal with missing methods gracefully.

 __len__()
     number of frames in trajectory
 __getitem__(arg)
     advance to time step arg=*frame* and return Timestep; or if arg is a
     slice, then return an interator over that part of the trajectory
 timeseries(atomGroup, [start[,stop[,skip[,format]]]])
     returns a subset of coordinate data
 correl(timeseriesCollection[,start[,stop[,skip]]])
     populate a :class:`MDAnalysis.core.Timeseries.TimeseriesCollection` object
     with observable timeseries computed from the trajectory

Attributes
..........

 filename
     filename of the trajectory
 numatoms
     number of atoms (coordinate sets) in a frame (constant)
 numframes
     total number of frames (if known) -- None if not known
 fixed
     bool, saying if there are fixed atoms (e.g. dcds)
 skip
     step size for iterating through the trajectory [1]
 skip_timestep
     number of integrator steps between frames + 1 (i.e.
     the stride at which the MD simulation was sampled)
 delta
     integrator time step (in native units); hence the "length" 
     of a trajctory frame is  skip_timestep*delta time units
 periodic
     contains box information for periodic boundary conditions (?)
 ts
     the :class:`~DCD.Timestep` object (typically customized for each
     trajectory format although at the moment all are derived from
     :class:`DCD.Timestep`.
 units
     dictionary with keys *time* and *length* and the appropriate 
     unit (e.g. 'AKMA' and 'Angstroem' for Charmm dcds, 'ps' and 'nm' 
     for Gromacs trajectories, ``None`` and 'Angstroem' for PDB).

Trajectory Writer
~~~~~~~~~~~~~~~~~

Methods
.......

 __init__(filename,numatoms[,start[,step[,delta[,remarks]]]])
     opens *filename* and writes header if required by format
 write_next_timestep([timestep])
     write data in *timestep* to trajectory file
 close_trajectory()
     close file and finish I/O 
 __del__()
     ensure that close_trajectory() is called

Attributes
..........
 filename
     name of the trajectory file
 start, stop, step
     first and last frame and step 
 units
     dictionary with keys *time* and *length* and the appropriate 
     unit (e.g. 'AKMA' and 'Angstroem' for Charmm dcds, 'ps' and 'nm' 
     for Gromacs trajectories, ``None`` and 'Angstroem' for PDB)


**Optional**
 
 ts
     Timestep instance      
"""

__all__ = ['DCD', 'PDB', 'CRD', 'XTC', 'TRR']

import PDB, DCD, CRD, XTC, TRR

# trajectory readers: present unified interface (based on DCD.Timestep)
_trajectory_readers = {'dcd': DCD.DCDReader,
                       'trj': DCD.DCDReader,
                       'xtc': XTC.XTCReader,
                       'trr': TRR.TRRReader,
                       'pdb': PDB.PDBReader,
                       }

# frame writers: export to single frame formats such as PDB, gro, crd
# Signature:
#   W = FrameWriter(filename)
#   W.write(AtomGroup)
_frame_writers = {'pdb': PDB.PrimitivePDBWriter,
                  'crd': CRD.CRDWriter,
                 }

# trajectory writers: export frames, typically only saving coordinates
# (although PDB movies are the exception); maybe look at OpenBabel as
# not to reinvent the wheel.
# Signature:
#   W = TrajectoryWriter(filename,numatoms,**kwargs)
#   W.write_next_timestep(TimeStep)
_trajectory_writers = {'dcd': DCD.DCDWriter,
                       'xtc': XTC.XTCWriter,
                       'trr': TRR.TRRWriter,
                       }

