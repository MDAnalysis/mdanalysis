.. -*- coding: utf-8 -*-
.. _reading-and-writing:

=========================
Reading and writing files
=========================

Input
=====

Read information from topology and coordinate files to create a :ref:`Universe <universe-loading>`::

    import MDAnalysis as mda
    u = mda.Universe('topology.gro', 'trajectory.xtc')


A topology file is always required for a Universe, whereas coordinate files are optional. :ref:`Some file formats provide both topology and coordinate information <format-overview>`. MDAnalysis supports a number of :ref:`formats <formats>`, which are automatically detected based on the file extension. For example, the above loads a :ref:`GROMACS XTC trajectory <XTC-format>`. Multiple coordinate files can be loaded, as :ref:`described below <chainreader>`; the following code loads two :ref:`CHARMM/NAMD DCD files <DCD-format>` and concatenates them::

    u = mda.Universe('topology.psf', 'trajectory1.dcd', 'trajectory2.dcd')
 
 
Some formats can be loaded with format-specific keyword arguments, such as the :ref:`LAMMPS DATA <DATA-format>` ``atom_style`` specification.

.. seealso::

    See :ref:`universe-loading` for more information on loading data into a Universe from files.


.. _chainreader:

--------------------------------
Reading multiple trajectories
--------------------------------

A Universe can load multiple trajectories, which are concatenated in the order given. One exception to this is with :ref:`XTC <XTC-format>` and :ref:`TRR <TRR-format>` files. If the ``continuous=True`` flag is :ref:`passed to Universe <universe-kwargs>`, MDAnalysis will try to stitch them together so that the trajectory is as time-continuous as possible. This means that there will be no duplicate time-frames, or jumps back in time.  

As an example, the following depicted trajectory is split into four parts. The column represents the time. As you can see, some segments overlap. With ``continuous=True``, only the frames marked with a + will be read.

::

    part01:  ++++--
    part02:      ++++++-
    part03:            ++++++++
    part04:                         ++++

However, there can be gaps in time (i.e. frames appear to be missing). Ultimately it is the user's responsiblity to ensure that the trajectories can be stitched together meaningfully.

.. note::

    While you can set ``continuous=True`` for either XTC or TRR files, you cannot mix different formats.

More information can be found at the API reference for :class:`~MDAnalysis.coordinates.chain.ChainReader`.

Trajectory formats
------------------

If no ``format`` keyword is provided, :class:`~MDAnalysis.coordinates.chain.ChainReader` will try to guess the format for each file from its extension. You can force :class:`~MDAnalysis.coordinates.chain.ChainReader` to use the same format for every file by using the ``format`` keyword. You can also specify which format to use by file, by passing in a sequence of ``(filename, format)`` tuples.

.. ipython:: python
    :okwarning:

    import MDAnalysis as mda
    from MDAnalysis.tests.datafiles import PDB, GRO

    u = mda.Universe(PDB, [(GRO, 'gro'), (PDB, 'pdb'), (GRO, 'gro')])
    u.trajectory




.. _memory-reader:

--------------------------------
In-memory trajectories
--------------------------------

Reading trajectories into memory
--------------------------------

If your device has sufficient memory to load an entire trajectory into memory, then analysis can be sped up substantially by transferring the
trajectory to memory. This makes it possible to
operate on raw coordinates using existing MDAnalysis tools. In
addition, it allows the user to make changes to the coordinates in a
trajectory (e.g. through
:attr:`AtomGroup.positions <MDAnalysis.core.groups.AtomGroup.positions>`) without having
to write the entire state to file.

The most straightforward way to do this is to pass ``in_memory=True`` to :class:`~MDAnalysis.core.universe.Universe`, which
automatically transfers a trajectory to memory:

.. ipython:: python

    from MDAnalysis.tests.datafiles import TPR, XTC

    universe = mda.Universe(TPR, XTC, in_memory=True)


MDAnalysis uses the :class:`~MDAnalysis.coordinates.memory.MemoryReader` class to load this data in.

Transferring trajectories into memory
-------------------------------------

The decision to transfer the trajectory to memory can be made at any
time with the
:meth:`~MDAnalysis.core.universe.Universe.transfer_to_memory` method
of a :class:`~MDAnalysis.core.universe.Universe`:

.. ipython:: python

    universe = mda.Universe(TPR, XTC)
    universe.transfer_to_memory()

This operation may take a while (passing ``verbose=True`` to :meth:`~MDAnalysis.core.universe.Universe.transfer_to_memory` will display a progress bar). However, subsequent operations on the trajectory will be very fast.

Building trajectories in memory
--------------------------------

:class:`~MDAnalysis.coordinates.memory.MemoryReader` can also be used to directly generate a trajectory as a numpy array.

.. ipython:: python

    from MDAnalysisTests.datafiles import PDB
    from MDAnalysis.coordinates.memory import MemoryReader
    import numpy as np

    universe = mda.Universe(PDB)
    universe.atoms.positions

The :meth:`~MDAnalysis.core.universe.Universe.load_new` method can be used to load coordinates into a Universe, replacing the old coordinates:

.. ipython:: python

    coordinates = np.random.rand(len(universe.atoms), 3)
    universe.load_new(coordinates, format=MemoryReader);
    universe.atoms.positions

or they can be directly passed in when creating a Universe.

.. ipython:: python

    universe2 = mda.Universe(PDB, coordinates, format=MemoryReader)
    universe2.atoms.positions


In-memory trajectories of an atom selection
-------------------------------------------

Creating a trajectory of an atom selection requires transferring the appropriate units. This is often needed when using :meth:`~MDAnalysis.core.universe.Merge` to create a new Universe, as coordinates are not automatically loaded in.



Output
======================

------------------------
Frames and trajectories
------------------------

MDAnalysis :class:`~MDAnalysis.core.universe.Universe`\ s can be written out to a :ref:`number of formats <coordinate-readers>` with :meth:`~MDAnalysis.core.groups.AtomGroup.write`. For example, to write the current frame as a PDB:

.. code-block :: python

    from MDAnalysis.tests.datafiles import PDB, TRR
    u = mda.Universe(PDB, TRR)
    ag = u.select_atoms("name CA")
    ag.write("c-alpha.pdb")

Pass in the ``frames`` keyword to write out trajectories.

.. code-block :: python

    ag.write('c-alpha_all.xtc', frames='all')

Slice or index the trajectory to choose which frames to write:

.. code-block :: python

    ag.write('c-alpha_skip2.trr', frames=u.trajectory[::2])
    ag.write('c-alpha_some.dcd', frames=u.trajectory[[0,2,3]])

Alternatively, iterate over the trajectory frame-by-frame with :func:`~MDAnalysis.Writer`. This requires you to pass in the number of atoms to write.

.. code-block :: python

    with mda.Writer('c-alpha.xyz', ag.n_atoms) as w:
        for ts in u.trajectory:
            w.write(ag)

You can pass keyword arguments to some format writers. For example, the :ref:`LAMMPS DATA <DATA-format>` format allows the ``lengthunit`` and ``timeunit`` keywords to specify the output units.  

Pickling
========

MDAnalysis currently supports pickling of AtomGroups and trajectories that have *not* been read in as :ref:`multiple files <chainreader>` or from a PDB file (`Issue 1981`_). Universe cannot be pickled.

.. ipython:: python

    import pickle
    from MDAnalysis.tests.datafiles import PSF, DCD
    psf = mda.Universe(PSF, DCD)
    pickle.loads(pickle.dumps(psf.trajectory))

While *trajectories* from PDB files cannot be pickled, trajectories where only the topology information comes from a PDB file *can*. For example, the universe below loads the trajectory information from a :ref:`TRR <TRR-format>` file.

.. ipython:: python

    u = mda.Universe(PDB, TRR)
    pickle.loads(pickle.dumps(u.trajectory))


.. _`Issue 1981`: https://github.com/MDAnalysis/mdanalysis/issues/1981>