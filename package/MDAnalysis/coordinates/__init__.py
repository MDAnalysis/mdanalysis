# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#


"""\
Trajectory Readers and Writers  --- :mod:`MDAnalysis.coordinates`
=================================================================

The coordinates submodule contains code to read, write and store coordinate
information, either single frames (e.g., the :mod:`~MDAnalysis.coordinates.GRO`
format) or trajectories (such as the :mod:`~MDAnalyis.coordinates.DCD` format);
see the :ref:`Supported coordinate formats` for all formats.

MDAnalysis calls the classes that read a coordinate trajectory and make the
data available *"Readers"*. Similarly, classes that write out coordinates are
called *"Writers"*. Readers and Writers provide a common interface to the
underlying coordinate data. This abstraction of coordinate access through an
object-oriented interface is one of the key capabilities of MDAnalysis.


.. _Readers:

Readers
-------

All Readers are based on a :class:`~MDAnalysis.coordinates.base.ProtoReader`
class that defines a common :ref:`Trajectory API` and allows other code to
interface with all trajectory formats in the same way, independent of the
details of the trajectory format itself.

The :class:`~MDAnalysis.core.universe.Universe` contains the API entry point
attribute :attr:`Universe.trajectory` that points to the actual
:class:`~MDAnalysis.coordinates.base.ProtoReader` object; all Readers are accessible
through this entry point in the same manner ("`duck typing`_").

There are three types of base Reader which act as starting points for each
specific format. These are:

:class:`~MDAnalysis.coordinates.base.ReaderBase`
   A standard multi frame Reader which allows iteration over a single
   file to provide multiple frames of data.  This is used by formats
   such as TRR and DCD.

:class:`~MDAnalysis.coordinates.base.SingleFrameReaderBase`
   A simplified Reader which reads a file containing only a single
   frame of information.  This is used with formats such as GRO
   and CRD

:class:`~MDAnalysis.coordinates.chain.ChainReader`
   An advanced Reader designed to read a sequence of files, to
   provide iteration over all the frames in each file seamlessly.
   This Reader can also provide this functionality over a
   sequence of files in different formats.

Normally, one does not explicitly need to select a reader. This is handled
automatically when creating a :class:`~MDAnalysis.core.universe.Universe` and
the appropriate reader for the file type is selected (typically by the file
extension but this choice can be overriden with the ``format`` argument to
:class:`~MDAnalysis.core.universe.Universe`).

If additional simulation data is available, it may be added to and read
alongside a trajectory using
:meth:`~MDAnalysis.coordinates.base.ProtoReader.add_auxiliary` as described in
the :ref:`Auxiliary API`.

.. _duck typing: http://c2.com/cgi/wiki?DuckTyping


.. _writing-trajectories:

Writers
-------

In order to **write coordinates**, a factory function is provided
(:func:`MDAnalysis.coordinates.core.writer`, which is also made available as
:func:`MDAnalysis.Writer`) that returns a *Writer* appropriate for the desired
file format (as indicated by the filename suffix). Furthermore, a trajectory
:class:`~MDAnalysis.coordinates.base.ProtoReader` can also have a method
:meth:`~MDAnalysis.coordinates.base.ProtoReader.Writer` that returns an appropriate
:class:`~MDAnalysis.coordinates.base.WriterBase` for the file format of the
trajectory.

In analogy to :func:`MDAnalysis.coordinates.core.writer`, there is also a
:func:`MDAnalysis.coordinates.core.reader` function available to return a
trajectory :class:`~MDAnalysis.coordinates.base.ProtoReader` instance although this
is often not needed because the :class:`~MDAnalysis.core.universe.Universe`
class can choose an appropriate reader automatically.

A typical approach is to generate a new trajectory from an old one, e.g., to
only keep the protein::

  u = MDAnalysis.Universe(PDB, XTC)
  protein = u.select_atoms("protein")
  with MDAnalysis.Writer("protein.xtc", protein.n_atoms) as W:
      for ts in u.trajectory:
          W.write(protein)

Using the :func:`with` statement will automatically close the trajectory when
the last frame has been written.


Timesteps
---------

Both Readers and Writers use Timesteps as their working object.  A
:class:`~MDAnalysis.coordinates.Timestep` represents all data for a given
frame in a trajectory.  The data inside a
:class:`~MDAnalysis.coordinates.Timestep` is often accessed indirectly
through a :class:`~MDAnalysis.core.groups.AtomGroup` but it is also possible to
manipulate Timesteps directly.

The current :class:`~MDAnalysis.coordinates.Timestep` can be accessed
through the :attr:`~MDAnalysis.coordinates.base.ProtoReader.ts` attribute of
the trajectory attached to the active
:class:`~MDAnalysis.core.universe.Universe`::

   ts = u.trajectory.ts
   ts.positions  # returns a numpy array of positions

Most individual formats have slightly different data available in each Timestep
due to differences in individual simulation packages, but all share in common a
broad set of basic data, detailed in `Timestep API`_


Supported coordinate formats
----------------------------

The table below lists the coordinate file formats understood by MDAnalysis. The
emphasis is on formats that are used in popular molecular dynamics codes. By
default, MDAnalysis figures out formats by looking at the extension. Thus, a
DCD file always has to end with ".dcd" to be recognized as such unless the
format is explicitly specified with the *format* keyword to
:class:`~MDAnalysis.core.universe.Universe` or
:meth:`~MDAnalysis.core.universe.Universe.load_new`.  A number of files are
also recognized when they are compressed with :program:`gzip` or
:program:`bzip2` such as ".xyz.bz2".

.. _Supported coordinate formats:

.. table:: Table of supported coordinate formats

   +---------------+-----------+-------+------------------------------------------------------+
   |Name           | extension |  IO   | remarks                                              |
   +===============+===========+=======+======================================================+
   | CHARMM,       | dcd       |  r/w  | standard CHARMM binary trajectory; endianness is     |
   | NAMD          |           |       | autodetected. Fixed atoms may not be handled         |
   |               |           |       | correctly (requires testing). Module                 |
   |               |           |       | :mod:`MDAnalysis.coordinates.DCD`                    |
   +---------------+-----------+-------+------------------------------------------------------+
   | LAMMPS        | dcd       |  r/w  | CHARMM-style binary trajectory; endianness is        |
   |               |           |       | autodetected. Units are appropriate for LAMMPS.      |
   |               |           |       | Module :mod:`MDAnalysis.coordinates.LAMMPS`          |
   +---------------+-----------+-------+------------------------------------------------------+
   | LAMMPS [#a]_  | data      |  r    | Single frame of coordinates read from .data files    |
   +---------------+-----------+-------+------------------------------------------------------+
   | LAMMPS [#a]_  | lammpsdump|  r    | Ascii trajectory in atom style                       |
   +---------------+-----------+-------+------------------------------------------------------+
   | Gromacs       | xtc       |  r/w  | Compressed (lossy) xtc trajectory format. Module     |
   |               |           |       | :mod:`MDAnalysis.coordinates.XTC`                    |
   +---------------+-----------+-------+------------------------------------------------------+
   | Gromacs       | trr       |  r/w  | Full precision trr trajectory. Coordinates and       |
   |               |           |       | velocities are processed. Module                     |
   |               |           |       | :mod:`MDAnalysis.coordinates.TRR`                    |
   +---------------+-----------+-------+------------------------------------------------------+
   | Gromacs       | tng       |  r    | Variable precision tng trajectory. Coordinates,      |
   |               |           |       | velocities and forces are processed along with any   |
   |               |           |       | `additional tng block data`_ requested for reading.  |
   |               |           |       | Uses the `PyTNG package`_ for tng file reading.      |  
   |               |           |       | Module :mod:`MDAnalysis.coordinates.TNG`             |
   +---------------+-----------+-------+------------------------------------------------------+
   | XYZ [#a]_     |  xyz      |  r/w  | Generic white-space separate XYZ format; can be      |
   |               |           |       | compressed (gzip or bzip2). Module                   |
   |               |           |       | :mod:`MDAnalysis.coordinates.XYZ`                    |
   +---------------+-----------+-------+------------------------------------------------------+
   | TXYZ [#a]_    |  txyz,    |  r    | Tinker XYZ format.                                   |
   |               |  arc      |       | Module :mod:`MDAnalysis.coordinates.TXYZ`            |
   +---------------+-----------+-------+------------------------------------------------------+
   | HOOMD [#a]_   |  gsd      |  r    | HOOMD GSD format (using :mod:`gsd.hoomd`).           |
   |               |           |       | Module :mod:`MDAnalysis.coordinates.GSD`             |
   +---------------+-----------+-------+------------------------------------------------------+
   | GAMESS        |  gms,     |  r    | Generic semi-formatted GAMESS output log; can be     |
   |               |  log,     |       | compressed (gzip or bzip2). Module                   |
   |               |  out      |       | :mod:`MDAnalysis.coordinates.GMS`                    |
   +---------------+-----------+-------+------------------------------------------------------+
   | AMBER         | trj,      |  r    | formatted (ASCII) trajectories; the presence of a    |
   |               | mdcrd     |       | periodic box is autodetected (*experimental*).       |
   |               |           |       | Module :mod:`MDAnalysis.coordinates.TRJ`             |
   +---------------+-----------+-------+------------------------------------------------------+
   | AMBER         | inpcrd,   | r     | formatted (ASCII) coordinate/restart file            |
   |               | restrt    |       | Module :mod:`MDAnalysis.coordinates.INPCRD`          |
   +---------------+-----------+-------+------------------------------------------------------+
   | AMBER         | ncdf, nc  |  r/w  | binary (NetCDF) trajectories are fully supported with|
   |               |           |       | optional `netcdf4-python`_ module (coordinates and   |
   |               |           |       | velocities). Module :mod:`MDAnalysis.coordinates.TRJ`|
   +---------------+-----------+-------+------------------------------------------------------+
   | Brookhaven    | pdb/ent   |  r/w  | a relaxed PDB format (as used in MD simulations)     |
   | [#a]_         |           |       | is read by default; Multiple frames (MODEL)          |
   |               |           |       | are supported but require the *multiframe* keyword.  |
   |               |           |       | Module :mod:`MDAnalysis.coordinates.PDB`             |
   +---------------+-----------+-------+------------------------------------------------------+
   | XPDB          | pdb       |   r   | Extended PDB format (can use 5-digit residue         |
   |               |           |       | numbers). To use, specify the format "XPBD"          |
   |               |           |       | explicitly: ``Universe(..., format="XPDB")``.        |
   |               |           |       | Module :mod:`MDAnalysis.coordinates.PDB`             |
   +---------------+-----------+-------+------------------------------------------------------+
   | PDBQT [#a]_   | pdbqt     | r/w   | file format used by AutoDock with atom types *t*     |
   |               |           |       | and partial charges *q*. Module:                     |
   |               |           |       | :mod:`MDAnalysis.coordinates.PDBQT`                  |
   +---------------+-----------+-------+------------------------------------------------------+
   | PQR [#a]_     | pqr       | r/w   | PDB-like but whitespace-separated files with charge  |
   |               |           |       | and radius information. Module                       |
   |               |           |       | :mod:`MDAnalysis.coordinates.PQR`                    |
   +---------------+-----------+-------+------------------------------------------------------+
   | GROMOS96      | gro       |  r/w  | basic GROMOS96 format (velocities as well).  Only    |
   |               |           |       | the first frame present will be read.                |
   | [#a]_         |           |       | Module :mod:`MDAnalysis.coordinates.GRO`             |
   +---------------+-----------+-------+------------------------------------------------------+
   | GROMOS11      | trc       |  r    | Basic GROMOS11 trajectory format.                    |
   |               |           |       | Can read positions, box-sizes and timesteps.         |
   |               |           |       | Module :mod:`MDAnalysis.coordinates.TRC`             |
   +---------------+-----------+-------+------------------------------------------------------+
   | CHARMM        | crd       |  r/w  | "CARD" coordinate output from CHARMM; deals with     |
   | CARD [#a]_    |           |       | either standard or EXTended format. Module           |
   |               |           |       | :mod:`MDAnalysis.coordinates.CRD`                    |
   +---------------+-----------+-------+------------------------------------------------------+
   | DESRES [#a]_  | dms       |  r    | DESRES Molecular Structure file format reader.       |
   |               |           |       | Module :mod:`MDAnalysis.coordinates.DMS`             |
   +---------------+-----------+-------+------------------------------------------------------+
   | IBIsCO/YASP   | trz       |  r/w  | Binary IBIsCO or YASP trajectories Module            |
   |               |           |       | :mod:`MDAnalysis.coordinates.TRZ`                    |
   +---------------+-----------+-------+------------------------------------------------------+
   | MOL2          | mol2      |  r/w  | Text-based Tripos molecular structure format         |
   |               |           |       | :mod:`MDAnalysis.coordinates.MOL2`                   |
   +---------------+-----------+-------+------------------------------------------------------+
   | DL_Poly [#a]_ | config    |  r    | DL_Poly ascii config file                            |
   |               |           |       | :mod:`MDAnalysis.coordinates.DLPOLY`                 |
   +---------------+-----------+-------+------------------------------------------------------+
   | DL_Poly [#a]_ | history   |  r    | DL_Poly ascii history file                           |
   |               |           |       | :mod:`MDAnalysis.coordinates.DLPOLY`                 |
   +---------------+-----------+-------+------------------------------------------------------+
   | MMTF [#a]_    | mmtf      |  r    | Macromolecular Transmission Format                   |
   |               |           |       | :mod:`MDAnalysis.coordinates.MMTF`                   |
   +---------------+-----------+-------+------------------------------------------------------+
   | NAMD          | coor,     |  r/w  | NAMD binary file format for coordinates              |
   |               | namdbin   |       | :mod:`MDAnalysis.coordinates.NAMDBIN`                |
   +---------------+-----------+-------+------------------------------------------------------+
   | FHIAIMS       | in        |  r/w  | FHI-AIMS file format for coordinates                 |
   |               |           |       | :mod:`MDAnalysis.coordinates.FHIAIMS`                |
   +---------------+-----------+-------+------------------------------------------------------+
   | H5MD          | h5md      |  r    | H5MD_ file format for coordinates                    |
   |               |           |       | :mod:`MDAnalysis.coordinates.H5MD`                   |
   +---------------+-----------+-------+------------------------------------------------------+
   | `chemfiles`_  | CHEMFILES |  r/w  | interface to `chemfiles`_, see the `list of chemfiles|
   | library       |           |       | file formats`_ and                                   |
   |               |           |       | :mod:`MDAnalysis.coordinates.chemfiles`              |
   +---------------+-----------+-------+------------------------------------------------------+

.. [#a] This format can also be used to provide basic *topology*
   information (i.e. the list of atoms); it is possible to create a
   full :mod:`~MDAnalysis.core.universe.Universe` by simply
   providing a file of this format: ``u = Universe(filename)``

.. _`netcdf4-python`: https://github.com/Unidata/netcdf4-python
.. _`H5MD`: https://nongnu.org/h5md/index.html
.. _`chemfiles`: https://chemfiles.org/
.. _`list of chemfiles file formats`: https://chemfiles.org/chemfiles/latest/formats.html
.. _`additional tng block data`: https://www.mdanalysis.org/pytng/documentation_pages/Blocks.html
.. _`PyTNG package`: https://github.com/MDAnalysis/pytng

.. _`Trajectory API`:

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


Registry
~~~~~~~~

In various places, MDAnalysis tries to automatically select appropriate formats
(e.g. by looking at file extensions). In order to allow it to choose the
correct format, all I/O classes must subclass either
:class:`MDAnalysis.coordinates.base.ReaderBase`,
:class:`MDAnalysis.coordinates.base.SingleFrameReaderBase`,
or :class:`MDAnalysis.coordinates.base.WriterBase` and set the
:attr:`~MDAnalysis.coordinates.base.ProtoReader.format` attribute with a string
defining the expected suffix.  To assign multiple suffixes to an I/O class, a
list of suffixes can be given.

In addition to this, a Reader may define a ``_format_hint`` staticmethod, which
returns a boolean of if it can process a given object. E.g. the
:class:`MDAnalysis.coordinates.memory.MemoryReader` identifies itself as
capable of reading numpy arrays.  This functionality is used in
:func:`MDAnalysis.core._get_readers.get_reader_for` when figuring out how to
read an object (which was usually supplied to mda.Universe).

To define that a Writer can write multiple trajectory frames, set the
`multiframe` attribute to ``True``.  The default is ``False``.
To define that a Writer *does not* support single frame writing the
`singleframe` attribute can be set to ``False``.  This is ``True``
by default, ie we assume all Writers can also do a single frame.


.. _Timestep API:

Timestep class
~~~~~~~~~~~~~~

A Timestep instance holds data for the current frame. It is updated whenever a
new frame of the trajectory is read.

Timestep classes are derived from
:class:`MDAnalysis.coordinates.timestep.Timestep`, which is the primary
implementation example (and used directly for the DCDReader).

The discussion on this format is detailed in `Issue 250`_

.. _Issue 250: https://github.com/MDAnalysis/mdanalysis/issues/250


Methods
.......

``__init__(n_atoms, positions=True, velocities=False, forces=False)``
    Define the number of atoms this Timestep will hold and whether or not
    it will have velocity and force information
``__eq__``
    Compares a Timestep with another
``__getitem__(atoms)``
    position(s) of atoms; can be a slice or numpy array and then returns
    coordinate array
``__len__()``
    number of coordinates (atoms) in the frame
``__iter__()``
    iterator over all coordinates
``copy()``
    deep copy of the instance
``_init_unitcell``
    hook that returns empty data structure for the unitcell representation
    of this particular file format; called from within ``__init__()`` to
    initialize :attr:`Timestep._unitcell`.

Attributes
..........

``n_atoms``
    number of atoms in the frame
``frame``
    current frame number (0-based)
``_frame``
    The native frame number of the trajectory.  This can differ from ``frame``
    as that will always count sequentially from 0 on iteration, whilst
    ``_frame`` is taken directly from the trajectory.
``time``
    The current system time in ps.  This value is calculated either from a time
    set as the Timestep attribute, or from `frame` * `dt`.  Either method allows
    allows an offset to be applied to the time.
``dt``
    The change in system time between different frames.  This can be set as an
    attribute, but defaults to 1.0 ps.
``data``
    A dictionary containing all miscellaneous information for the
    current Timestep.
``positions``
    A numpy array of all positions in this Timestep, otherwise raises a
    :class:`~MDAnalysis.exceptions.NoDataError`
``velocities``
    If present, returns a numpy array of velocities, otherwise raises a
    :class:`~MDAnalysis.exceptions.NoDataError`
``forces``
    If present, returns a numpy array of forces, otherwise raises a
    :class:`~MDAnalysis.exceptions.NoDataError`
``has_positions``
    Boolean of whether position data is available
``has_velocities``
    Boolean of whether velocity data is available
``has_forces``
    Boolean of whether force data is available
``dimensions``
    system box dimensions (`x, y, z, alpha, beta, gamma`)
    Also comes with a setter that takes a MDAnalysis box so that one can do ::

        Timestep.dimensions = [A, B, C, alpha, beta, gamma]

    which then converts automatically to the underlying representation and stores it
    in :attr:`Timestep._unitcell`.
``volume``
    system box volume (derived as the determinant of the box vectors of ``dimensions``)
``aux``
    namespace for the representative values of any added auxiliary data.


Private attributes
..................

These attributes are set directly by the underlying trajectory
readers. Normally the user should not have to directly access those,
but instead should use the attribute above.

``_pos``
    raw coordinates, a :class:`numpy.float32` array; ``X = _pos[:,0], Y =
    _pos[:,1], Z = _pos[:,2]``

``_velocities``
    raw velocities, a :class:`numpy.float32` array containing velocities
    (similar to ``_pos``)

``_forces``
    forces, similar to velocities above.

``_unitcell``
    native unit cell description; the format depends on the
    underlying trajectory format. A user should use the
    :class:`~MDAnalysis.coordinates.Timestep.dimensions`
    attribute to access the data in a canonical format instead of
    accessing :class:`Timestep._unitcell` directly.

    The method :meth:`Timestep._init_unitcell` is a hook to initialize
    this attribute.



Trajectory Reader class
~~~~~~~~~~~~~~~~~~~~~~~

Trajectory readers are derived from
:class:`MDAnalysis.coordinates.base.ReaderBase` (or from
:class:`MDAnalysis.coordinates.base.ProtoReader` if they do not required
:meth:`Reader.__del__` method). A special case are *SingleFrame readers* for
formats that contain only a single coordinate frame. These readers are derived
from a subclass of :class:`~MDAnalysis.coordinates.base.ProtoReader` named
:class:`MDAnalysis.coordinates.base.SingleFrameReaderBase`.

Typically, many methods and attributes are overriden but the ones listed below
*must* be implemented.

.. SeeAlso::

   See the section on :ref:`ReadersBase` in :mod:`MDAnalysis.coordinates.base`
   for implementation details.


Methods
.......

The :class:`MDAnalysis.coordinates.DCD.DCDReader` class is the primary
implementation example.

**Mandatory methods**

The following methods must be implemented in a Reader class.

``__init__(filename, **kwargs)``
    open *filename*; other *kwargs* are processed as needed and the
    Reader is free to ignore them. Typically, when MDAnalysis creates
    a Reader from :class:`MDAnalysis.Universe` it supplies as much
    information as possible in `kwargs`; at the moment the following
    data are supplied:

    - *n_atoms*: the number of atoms from the supplied topology.  This is
                not required for all readers and can be ignored if not
                required.

``__iter__()``
    allow iteration from beginning to end::

        for ts in trajectory:
            print(ts.frame)

    Readers will automatically rewind the trajectory to before the initial
    frame (often by re-opening the file) before starting the iteration. *Multi
    frame readers* (see :ref:`ReadersBase`) will also rewind the trajectory
    *after* the iteration so that the current trajectory frame is set to the
    first trajectory frame. *Single frame readers* do not explicitly rewind
    after iteration but simply remain on the one frame in the trajectory.

``close()``
    close the file and cease I/O
``next()``
    advance to next time step or raise :exc:`IOError` when moving
    past the last frame
``rewind()``
    reposition to first frame
``__entry__()``
    entry method of a `Context Manager`_ (returns self)
``__exit__()``
    exit method of a `Context Manager`_, should call ``close()``.

.. _Context Manager: http://docs.python.org/2/reference/datamodel.html#context-managers

.. Note::
   a ``__del__()`` method should also be present to ensure that the
   trajectory is properly closed. However, certain types of Reader can ignore
   this requirement. These include the :class:`SingleFrameReaderBase` (file reading
   is done within a context manager and needs no closing by hand) and the :class:`ChainReader`
   (it is a collection of Readers, each already with its own ``__del__`` method).

**Optional methods**

Not all trajectory formats support the following methods, either because the
data are not available or the methods have not been implemented. Code should
deal with missing methods gracefully.

``__len__()``
    number of frames in trajectory

``__getitem__(arg)``
    advance to time step `arg` = `frame` and return :class:`Timestep`; or if `arg` is a
    slice, then return an iterable over that part of the trajectory.

    The first functionality allows one to randomly access frames in the
    trajectory::

        universe.trajectory[314]

    would load frame 314 into the current :class:`Timestep`.

    Using slices allows iteration over parts of a trajectory ::

        for ts in universe.trajectory[1000:2000]:
            process_frame(ts)   # do some analysis on ts

    or skipping frames ::

        for ts in universe.trajectory[1000::100]:
            process_frame(ts)   # do some analysis on ts

    The last example starts reading the trajectory at frame 1000 and
    reads every 100th frame until the end.

    A sequence of indices or a mask of booleans can also be provided to index
    a trajectory.

    The performance of the ``__getitem__()`` method depends on the underlying
    trajectory reader and if it can implement random access to frames. All
    readers in MDAnalysis should support random access.

    For external custom readers this may not be easily (or reliably)
    implementable and thus one is restricted to sequential iteration.
    If the Reader is not able to provide random access to frames then it
    should raise :exc:`TypeError` on indexing. It is possible to partially
    implement ``__getitem__`` (as done on
    :class:`MDAnalysis.coordinates.base.ProtoReader.__getitem__` where slicing the
    full trajectory is equivalent to
    :class:`MDAnalysis.coordinates.base.ProtoReader.__iter__` (which is always
    implemented) and other slices raise :exc:`TypeError`.

    When indexed with a slice, a sequence of indices, or a mask of booleans,
    the return value is an instance of :class:`FrameIteratorSliced` or
    :class:`FrameIteratorIndices`. See :ref:`FrameIterators` for more details.

``parse_n_atoms(filename, **kwargs)``
    Provide the number of atoms in the trajectory file, allowing the Reader
    to be used to provide an extremely minimal Topology.
    Must be implemented as either a staticmethod or a classmethod.

``Writer(filename, **kwargs)``
    returns a :class:`~MDAnalysis.coordinates.base.WriterBase` which is set up with
    the same parameters as the trajectory that is being read (e.g. time step,
    length etc), which facilitates copying and simple on-the-fly manipulation.

    If no Writer is defined then a :exc:`NotImplementedError` is raised.

    The *kwargs* can be used to customize the Writer as they are typically
    passed through to the init method of the Writer, with sensible defaults
    filled in; the actual keyword arguments depend on the Writer.

``timeseries(atomGroup, [start[,stop[,skip[,format]]]])``
    returns a subset of coordinate data

Attributes
..........

``filename``
    filename of the trajectory
``n_atoms``
    number of atoms (coordinate sets) in a frame (constant)
``n_frames``
    total number of frames (if known) -- ``None`` if not known
``ts``
    the :class:`~timestep.Timestep` object; typically customized for each
    trajectory format and derived from :class:`timestep.Timestep`.
``units``
    dictionary with keys *time*, *length*, *speed*, *force* and the
    appropriate unit (e.g. 'AKMA' and 'Angstrom' for Charmm dcds, 'ps' and
    'nm' for Gromacs trajectories, ``None`` and 'Angstrom' for PDB).
    Any field not used should be set to ``None``.
``format``
    string that identifies the file format, e.g. "DCD", "PDB", "CRD", "XTC",
    "TRR"; this is typically the file extension in upper case.
``dt``
    time between frames in ps; a managed attribute (read only) that computes
    on the fly ``skip_timestep * delta`` and converts to the MDAnalysis base
    unit for time (pico seconds by default)
``totaltime``
    total length of the trajectory = ``n_frames * dt``
``time``
    time of the current time step, in MDAnalysis time units (ps)
``frame``
    frame number of the current time step (0-based)
``aux_list``
    list of the names of any added auxiliary data.
``_auxs``
    dictionary of the :class:`~MDAnalysis.auxiliary.base.AuxReader`
    instances for any added auxiliary data.

**Optional attributes**

``delta``
    integrator time step (in native units); hence the "length"
    of a trajctory frame is  ``skip_timestep*delta`` time units
``compressed``
    string that identifies the compression (e.g. "gz" or "bz2") or ``None``.
``fixed``
    bool, saying if there are fixed atoms (e.g. dcds)
``periodic``
    boolean saying if contains box information for periodic boundary conditions
    unit cell information is stored in attribute `dimensions`
``skip_timestep``
    number of integrator steps between frames + 1 (i.e.
    the stride at which the MD simulation was sampled)


Trajectory Writer class
~~~~~~~~~~~~~~~~~~~~~~~

Trajectory writers are derived from
:class:`MDAnalysis.coordinates.base.WriterBase`. They are used to write
multiple frames to a trajectory file. Every time the
:meth:`~MDAnalysis.coordinates.base.WriterBase.write` method is called,
another frame is appended to the trajectory.

Typically, many methods and attributes are overriden.

Signature::

   with TrajectoryWriter(filename, n_atoms, **kwargs) as w:
       w.write(Universe)    # write a whole universe

or::

   w.write(AtomGroup)  # write a selection of Atoms from Universe

.. SeeAlso::

   See the section on :ref:`WritersBase` in :mod:`MDAnalysis.coordinates.base`
   for implementation details.

Methods
.......

``__init__(filename, n_atoms, **kwargs)``
    Set-up the reader. This *may* open file *filename* and *may*
    write content to it such as headers immediately but the writer is
    allowed to delay I/O up to the first call of ``write()``.

    Any ``**kwargs`` that are not processed by the writer must be
    silently ignored.

``write(obj)``
    write Timestep data in *obj*
``convert_dimensions_to_unitcell(timestep)``
    take the dimensions from the timestep and convert to the native
    unitcell representation of the format
``close()``
    close file and finish I/O
``__del__()``
    ensures that close() is called

Attributes
..........

``filename``
    name of the trajectory file
``units``
    dictionary with keys *time*, *length*, *speed*, *force* and the
    appropriate unit (e.g. 'AKMA' and 'Angstrom' for Charmm dcds, 'ps' and
    'nm' for Gromacs trajectories, ``None`` and 'Angstrom' for PDB).
    Any field not used should be set to ``None``.
``format``
    string that identifies the file format, e.g. "DCD", "PDB", "CRD", "XTC",
    "TRR"


**Optional**

``ts``
    :class:`Timestep` instance


Single Frame Writer class
~~~~~~~~~~~~~~~~~~~~~~~~~

A single frame writer is a special case of a trajectory writer in that it
writes only a single coordinate frame to a file, for instance, a pdb or gro
file. Unlike trajectory formats, which only contains coordinates, *single
frame* formats contain much more information (e.g. atom and residue names and
numbers) and hence it is possible to write selections of atoms in a meaningful
way.

Signature::

   W = FrameWriter(filename, **kwargs)
   W.write(AtomGroup)
   W.write(Universe)

The blanket *kwargs* is required so that one can pass the same kind of
arguments (filename and n_atoms) as for the Trajectory writers. In
this way, the simple :func:`~MDAnalysis.coordinates.core.writer`
factory function can be used for all writers.

Methods
.......

``__init__(filename, **kwargs)``
    opens *filename* for writing; `kwargs` are typically ignored
``write(obj)``
    writes the object *obj*, containing a
    :class:`~MDAnalysis.core.groups.AtomGroup` group of atoms (typically
    obtained from a selection) or :class:`~MDAnalysis.core.universe.Universe`
    to the file and closes the file

.. Note::

   Trajectory and Frame writers can be used in almost exactly the same
   manner with the one difference that Frame writers cannot deal with
   raw :class:`~MDAnalysis.coordinates.Timestep` objects.

"""
__all__ = ['reader', 'writer', 'timestep']

from . import base
from . import timestep
from .core import reader, writer
from . import chain
from . import chemfiles
from . import CRD
from . import DCD
from . import DLPoly
from . import DMS
from . import GMS
from . import GRO
from . import INPCRD
from . import LAMMPS
from . import MOL2
from . import PDB
from . import PDBQT
from . import PQR
from . import TRC
from . import TRJ
from . import TRR
from . import H5MD
from . import TRZ
from . import XTC
from . import XYZ
from . import TXYZ
from . import memory
from . import MMTF
from . import GSD
from . import null
from . import NAMDBIN
from . import FHIAIMS
from . import TNG
from . import MMCIF