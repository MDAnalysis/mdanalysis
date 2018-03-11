.. _load_ibisco:

#############################################
Loading IBIsCO and YASP files with MDAnalysis
#############################################

MDAnalysis is able to read `IBIsCO`_ / `YASP`_ TRZ binary trajectories, including
the coordinates and velocities as well as:

 - ``pressure``
 - ``pressure_tensor``
 - ``total_energy``
 - ``potential_energy``
 - ``kinetic_energy``
 - ``temperature``

Which are stored in the ``Timestep.data`` dictionary.  For example to view the kinetic
energy for each timestep::

  for ts in u.trajectory:
      print(ts.data['kinetic_energy'])

Trajectory data are read and written in binary representation but because this depends on
the machine hardware architecture, MDAnalysis *always* reads and writes TRZ
trajectories in *little-endian* byte order.

For the implementation details, see :mod:`MDAnalysis.coordinates.TRZ`.

There is currently no support for reading a topology for this format, however the Reader
is able to determine the number of atoms and provide a minimal topology (ie no names etc,
only indices of atoms).  To get a topology, it is required to convert an input into an
XYZ or similar format.

.. _IBIsCO: http://www.theo.chemie.tu-darmstadt.de/ibisco/IBISCO.html
.. _YASP: http://www.theo.chemie.tu-darmstadt.de/group/services/yaspdoc/yaspdoc.html
