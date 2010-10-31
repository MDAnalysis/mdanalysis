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

In order to write coordinates, a factory function is provided
(:func:`MDAnalysis.core.writer`) which is made available as
:func:`MDAnalysis.Writer`) that returns a Writer appropriate for the desired
file format.



Trajectory API
--------------

The **Trajectory API** defines how classes have to be structured that allow
reading and writing of coordinate files. By following the API it is possible to
seamlessly enhance the I/O capabilities of MDAnalysis. The actual underlying
I/O code can be written in C or python or a mixture thereof.

Typically, each format resides in its own module, named by the format specifier
(and using upper case by convention).

Reader and Writer classes are derived from base classes in
:mod:`MDAnalysis.coordinates.base`.


History
~~~~~~~

- 2010-04-30 Draft [orbeckst]
- 2010-08-20 added single frame writers to API [orbeckst]
- 2010-10-09 added write() method to Writers [orbeckst]
- 2010-10-19 use close() instead of close_trajectory() [orbeckst]
- 2010-10-30 clarified Writer write() methods (see also `Issue 49`_)

.. _Issue 49: http://code.google.com/p/mdanalysis/issues/detail?id=49


Registry
~~~~~~~~

In various places, MDAnalysis tries to automatically select appropriate formats
(e.g. by looking at file extensions). In order to allow it to choose the
correct format, all I/O classes must be registered in one of three dictionaries
with the format (typically the file extension in upper case):

- Trajectory reader classes must be added to
  :data:`MDAnalysis.coordinates._trajectory_readers`.

- Trajectory writer classes must be added to
  :data:`MDAnalysis.coordinates._trajectory_writers`.

- Single-frame writer classes must be added to to
  :data:`MDAnalysis.coordinates._frame_writers`.


Timestep class
~~~~~~~~~~~~~~

A Timestep instance holds data for the current frame. It is updated whenever a
new frame of the trajectory is read. 

Timestep classes are derived from
:class:`MDAnalysis.coordinates.base.Timestep`, which is the primary
implementation example (and used directly for the DCDReader).

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
      (typically implemented as a property because it needs to translate whatever is in the
      underlying :attr:`Timestep._unitcell` attribute)

Private attributes
..................

These attributes are set directly by the underlying trajectory
readers. Normally the user should not have to directly access tose (although in
some cases it is convenient to directly use :attr:`Timestep._pos`).

  _pos
      raw coordinates, a numpy.float32 array; X = _pos[:,0], Y =
      _pos[:,1], Z = _pos[:,2]

  _unitcell
      unit cell dimensions and angles; the format depends on the underlying
      trajectory format. A user should use :attr:`Timestep.dimensions` to
      access the data in a standard format.


Trajectory Reader class
~~~~~~~~~~~~~~~~~~~~~~~

Trajectory readers are derived from :class:`MDAnalysis.coordinates.base.Reader`.
Typically, many methods and attributes are overriden.

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
 close()
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
     contains box information for periodic boundary conditions (???)
 ts
     the :class:`~base.Timestep` object; typically customized for each
     trajectory format and derived from :class:`base.Timestep`.
 units
     dictionary with keys *time* and *length* and the appropriate 
     unit (e.g. 'AKMA' and 'Angstrom' for Charmm dcds, 'ps' and 'nm' 
     for Gromacs trajectories, ``None`` and 'Angstrom' for PDB).
 format
     string that identifies the file format, e.g. "DCD", "PDB", "CRD", "XTC",
     "TRR"; this is typically the file extension in upper case.
 dt
     time between frames in ps; a managed attribute (read only) that computes
     on the fly skip_timestep * delta and converts to the MDAnalysis base
     unit for time (pico seconds by default)
 totaltime
     total length of the trajectory = numframes * dt

**Optional attributes**

compressed
     string that identifies the compression (e.g. "gz" or "bz2") or ``None``.


Trajectory Writer class
~~~~~~~~~~~~~~~~~~~~~~~

Trajectory readers are derived from
:class:`MDAnalysis.coordinates.base.Writer`. They are use to write multiple
frames to a trajectory file. Every time the write() method is called, another
frame is appended to the trajectory.

Typically, many methods and attributes are overriden.

Signature::
   W = TrajectoryWriter(filename,numatoms,**kwargs)
   W.write_next_timestep(Timestep)
or::
   W.write(AtomGroup)   # write a selection
   W.write(Universe)    # write a whole universe
   W.write(Timestep)    # same as write_next_timestep()


Methods
.......

 __init__(filename,numatoms[,start[,step[,delta[,remarks]]]])
     opens *filename* and writes header if required by format
 write(obj)
     write Timestep data in *obj*
 write_next_timestep([timestep])
     write data in *timestep* to trajectory file
 convert_dimensions_to_unitcell(timestep)
     take the dimensions from the timestep and convert to the native
     unitcell representation of the format
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
     unit (e.g. 'AKMA' and 'Angstrom' for Charmm dcds, 'ps' and 'nm' 
     for Gromacs trajectories, ``None`` and 'Angstrom' for PDB)
 format
     string that identifies the file format, e.g. "DCD", "PDB", "CRD", "XTC",
     "TRR"


**Optional**
 
 ts
     Timestep instance


Single Frame Writer class
~~~~~~~~~~~~~~~~~~~~~~~~~

A single frame writer is a special case of a trajectory writer in that it
writes only a single coordinate frame to a file, for instance, a pdb or gro
file. Unlike trajectory formats, which only contains coordinates, "single
frame" formats contains much more information (e.g. atom and residue names and
numbers) and hence it is possible to write selections of atoms in a meaningful
way.

Signature::
   W = FrameWriter(filename, **kwargs)
   W.write(AtomGroup)
   W.write(Universe)

The blanket kwargs is required so that one can pass the same kind of
arguments (filename and numatoms) as for the Trajectory writers. In
this way, the simple :func:`writer` factory function can be used for
all writers.

Methods
.......
 __init__(filename, **kwargs)
   opens *filename* for writing; kwargs are typically ignored
 write(obj)
   writes the object *obj*, containing a
   :class:`~MDAnalysis.core.AtomGroup.AtomGroup` group of atoms (typically
   obtained from a selection) or :class:`~MDAnalysis.core.AtomGroup.Universe`
   to the file and closes the file

.. Note:: Trajectory and Frame writers can be used in almost exactly the same
   manner with the one difference that Frame writers cannot deal with raw
   :class:`~MDAnalysis.coordinates.base.Timestep` objects.

"""

__all__ = ['reader', 'writer']

import PDB, DCD, CRD, XTC, TRR, GRO, XYZ
import base
from core import reader, writer

#: trajectory readers: present unified interface (based on DCD.Timestep)
_trajectory_readers = {'DCD': DCD.DCDReader,
                       'TRJ': DCD.DCDReader,
                       'XTC': XTC.XTCReader,
                       'XYZ': XYZ.XYZReader,
                       'TRR': TRR.TRRReader,
                       'PDB': PDB.PDBReader,
                       'CRD': CRD.CRDReader,
                       'GRO': GRO.GROReader,
                       'CHAIN': base.ChainReader,
                       }
#: readers of files that contain both topology/atom data and coordinates
#: (currently only the keys are used)
_topology_coordinates_readers = {
                       'PDB': PDB.PrimitivePDBReader,
                       'GRO': GRO.GROReader,
                       'CRD': CRD.CRDReader,
}    

#: hack: readers that ignore most errors (permissive=True)
_trajectory_readers_permissive = _trajectory_readers.copy()
_trajectory_readers_permissive['PDB'] =  PDB.PrimitivePDBReader

#: frame writers: export to single frame formats such as PDB, gro, crd
# Signature:
#   W = FrameWriter(filename)
#   W.write(AtomGroup)
_frame_writers = {'PDB': PDB.PrimitivePDBWriter,
                  'CRD': CRD.CRDWriter,
                  'GRO': GRO.GROWriter,
                 }

# trajectory writers: export frames, typically only saving coordinates
# (although PDB movies are the exception); maybe look at OpenBabel as
# not to reinvent the wheel.
# Signature:
#   W = TrajectoryWriter(filename,numatoms,**kwargs)
#   W.write_next_timestep(TimeStep)
#   W.write(Timestep)
#   W.write(AtomGroup)
#   W.write(Universe)
_trajectory_writers = {'DCD': DCD.DCDWriter,
                       'XTC': XTC.XTCWriter,
                       'TRR': TRR.TRRWriter,
                       }

