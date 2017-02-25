.. -*- coding: utf-8 -*-

.. Working with coordinates and things specific to MDA
   Maybe a basic overview of numpy arrays, but nothing too much
   Mostly dealing with lib.distances and ways to abuse it!

Coordinate manipulation in MDAnalysis
=====================================

.. reminder of how cool numpy arrays are

Coordinate, velocity and force data is the meat of molecular simulation and
these are made directly available as properties of an AtomGroup.  These
properties will always be available from the AtomGroup, however if the loaded
trajectory does not have this data then a ``NoDataError`` will be raised.

::

     >>> import MDAnalysis as mda
     >>> from MDAnalysisTests.datafiles import PSF, DCD
     >>> ag = u.select_atoms('name CA')
     >>> ag.positions
     array([[ 11.66462231,   8.39347267,  -8.98323059],
            [ 11.41483879,   5.43442154,  -6.51348448],
          ...
     >>> vel = ag.velocities
     NoDataError: Timestep does not contain velocities

These properties return a **copy** of the data for the AtomGroup at the
timestep currently loaded by the trajectory Reader.  This means that the
arrays returned can be freely manipulated without changing the underlying data.
The arrays returned by the ``positions``, ``velocities`` and ``forces`` methods
are numpy arrays.  They can therefore immediately be used with any tool from
the vast existing scientific Python ecosystem.


Periodic boundary conditions
----------------------------

.. coordinate space is a flat circle

One of the main areas where existing tools are not adequate is in calculating
distances between atoms.  This is because of the periodic boundary conditions
(PBC) necessary to perform the molecular simulation.  PBC mean that the
primary unit cell is surrounded by virtual images of itself in all directions,
so that there is no "edge" to the simulation volume.  To properly evaluate
distances between points in our simulation space we must always consider if
there is a shorter distance existing by using the virtual images.

To support this most AtomGroup methods have a ``pbc=True/False`` to control
the consideration of PBC.  When dealing with coordinate arrays the box
information (available as the ``dimensions`` property) must be passed to the
relevent function using the ``box=`` keyword. Functions for handling distance
calculations in MDAnalysis are in ``MDAnalysis.lib.distances``.
::

     >>> import MDAnalysis as mda
     >>> from MDAnalysisTests.datafiles import GRO, TRR
     >>> from MDAnalysis.lib import distances
     >>> u = mda.Universe(GRO, TRR)
     >>> N = u.select_atoms('type N')
     >>> O = u.select_atoms('type O')
     >>> # Checking the effect of PBC on the maximum distance
     >>> distances.distance_array(N.positions, O.positions).max()
     116.73172083366822
     >>> distances.distance_array(N.positions, O.positions, box=u.dimensions).max()
     63.136479983880214


Overview of available functions
-------------------------------

For a full description of the functions see :mod:`MDAnalysis.lib.distances`.
All of these functions accept an optional box argument to take periodic
boundary conditions into account when calculating distances between positions.
Note that this is important when calculating angles too, as the correct
vector between positions must be considered.  Either orthogonal (rectangular)
or a triclinic box can be used.


* :func:`~MDAnalysis.lib.distances.distance_array` --
  For two arrays of coordinates `N` and `M`, calculates all combinations
  of pairwise distances, returning a N x M array.
* :func:`~MDAnalysis.lib.distances.self_distance_array` --
  For a single array of coordinates returns a flattened array of all
  combinations pairwise distances.
* :func:`~MDAnalysis.lib.distances.calc_bonds` --
  For two arrays of coordinates `N` and `M` of equal length, returns the
  distance between the first two items, second two items, etc.
* :func:`~MDAnalysis.lib.distances.calc_angles` --
  For three arrays of coordinates of equal length, calculates the
  calculates the angle between each triplet of positions, with the
  middle array acting as the apex of each angle.
* :func:`~MDAnalysis.lib.distances.calc_torsions` --
  For four arrays of coordinates of equal length, calculates the
  angle between the planes formed by the first three positions
  and last three positions.
* :func:`~MDAnalysis.lib.distances.apply_PBC` --
  Moves a set of coordinates to within the primary unit cell.  For a
  coordinate ``x`` in a dimension of size ``Lx``, this is defined as
  being ``0 <= x < Lx``.  
* :func:`~MDAnalysis.lib.distances.transform_RtoS` and
  :func:`~MDAnalysis.lib.distances.transform_StoR` --
  Transforms coordinates from real to fractional (S) space, and back
  again.


