.. _load_tinker:

##################################
Loading Tinker files in MDAnalysis
##################################

MDAnalysis can read the ascii trajectories produced by Tinker_.  Although these
are referred to as XYZ abd ARC files, these deviate from the XYZ standard
(featuring bonds
on the end of each line) and so are referred to as TXYZ files in MDAnalysis.
To read these files, the filename must end in either ``txyz`` or ``arc``, or
alternatively these format hints can be given to :class:`Universe` on creation
using the `format` keyword.

.. _Tinker: https://dasher.wustl.edu/tinker/

.. _load_txyz:


Reading TXYZ files
------------------

When used as a topology file, TXYZ files give the following attributes:
 - atomids (index as given in file)
 - names
 - types (numerical atom type)
 - bonds
 - masses (guessed from atom names)

When used as a coordinate file, MDAnalysis can read both single frame and multi frame
(archive) files.  This will yield information on the positions of all atoms for each
frame.

.. note::
   Tinker XYZ files have no information on the box size for each frame,
   therefore this will be 0 throughout the trajectory within MDAnalysis.

For details of the implementation of the Parser see :mod:`MDAnalysis.topology.TXYZParser`
and for details of the coordinate Reader see :mod:`MDAnalysis.coordinates.TXYZ`.
