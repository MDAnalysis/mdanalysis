.. _loading_gromacs: 

###################################
Loading Gromacs files in MDAnalysis
###################################

MDAnalysis is able to read most files used by Gromacs with the exception
of the TNG format which is currently in development.
The TPR file provides a wide variety of topology attributes
and should be used where possible,
however it is also possible to use either a GRO or similar ASCII file.
For the trajectory either an XTC or TRR file can be used,
with the latter also providing velocities and forces if these were saved.

.. seealso::
   `Loading PDB files <load_databank>`
   `Loading XYZ files <load_xyz>`
   

.. _load_tpr:

TPR files
---------

MDAnalysis supports reading TPR portable topology files created by
version 3.3 of Gromacs onwards; in particular this includes files from
versions 4.x, 5.x, 2016.x, and 2018.x. TPR files are binary files in
XDR format that, despite being binary, are portable between different
computing architectures.

Full topology information is read, including names, bonds, masses, and
even molecule information.

For details see :mod:`MDAnalysis.topology.TPRParser`.


.. _load_gro:


GRO files
---------

GRO files are ASCII files that contain coordinate (and potentially
velocity) information as well as the unit cell.

.. Note:: GRO files can also contain trajectories in the form of
	  multiple coordinate frames but GRO trajectories are
	  currently not supported.


Reading GRO files
~~~~~~~~~~~~~~~~~

MDAnalysis can use a GRO file both as a topology file and a
coordinate/velocity file.

When used as a topology file in a Universe, MDAnalysis will read
``names``, ``ids``, ``resids`` and ``resnames``,
and guess ``masses`` and ``types`` based on the names of atoms.
The ``segid`` for of all atoms is set to "``SYSTEM``".

For implementation details see
:mod:`MDAnalysis.topology.GROParser` and :mod:`MDAnalysis.coordinates.GRO`.



Writing GRO files
^^^^^^^^^^^^^^^^^

MDAnalysis can also write single trajectory frames to GRO files, for
example using the :meth:`~MDAnalysis.core.groups.AtomGroup.write`
method.  It will attempt to use the ``name``, ``resname`` and
``resid`` attribute from atoms if available, otherwise default values
of "``X``", "``UNK``" and "``1``" respectively will be used.  If the
provided atoms have velocities these will also be written, otherwise
only the positions will be written.  The precision is limited to three
decimal places.

By default any written GRO files will renumber the atom ids to move sequentially
from 1.  This can be disabled, and instead the original atom ids kept, by
using the ``reindex=False`` keyword argument.  This is useful when writing a
subsection of a larger Universe while wanting to preserve the original
identities of atoms.

For example::

   >>> u = mda.Universe()`

   >>> u.atoms.write('out.gro', reindex=False)

   # OR
   >>> with mda.Writer('out.gro', reindex=False) as w:
   ...     w.write(u.atoms)




.. Links
.. _Gromacs: http://www.gromacs.org
.. _`Gromacs manual`: http://manual.gromacs.org/documentation/5.1/manual-5.1.pdf
.. _TPR file: http://manual.gromacs.org/current/online/tpr.html
.. _`Issue Tracker`: https://github.com/MDAnalysis/mdanalysis/issues
.. _`Issue 2`: https://github.com/MDAnalysis/mdanalysis/issues/2
.. _`Issue 463`: https://github.com/MDAnalysis/mdanalysis/pull/463
.. _TPRReaderDevelopment: https://github.com/MDAnalysis/mdanalysis/wiki/TPRReaderDevelopment

.. _load_trr:

TRR and XTC files
-----------------

Gromacs TRR and XTC files are binary trajectory files. TRR files can
contain positions, velocities and/or forces at full precision as well
as time information, box dimensions, and values of the free energy
perturbation parameter :math:`\lambda` (which is 0 if it is not
used). XTC files only contain positions (together with time and box
information) and are compressed with a lossy algorithm, which reduces
the precision of coordinates to typically 0.001 nm (although this can
be set when running the simulation). MDAnalysis fully supports reading
and writing XTC and TRR trajectories.

.. SeeAlso::
   MDAnalysis.coordinates.XTC
   MDAnalysis.coordinates.TRR
   MDAnalysis.lib.formats.libmdaxdr


Reading TRR and XTC files
~~~~~~~~~~~~~~~~~~~~~~~~~

All information present in the trajectory is read and made available
in the :class:`Timestep`, either as `Timestep.positions`,
`Timestep.velocities`, or `Timestep.forces` (or the corresponding
attributes of AtomGroups). The FEP lambda is stored as an entry in the
:attr:`Timestep.data` dictionary (``Timestep.data['lambda']``).

Normally, XTC and TRR formats are not random access
formats. MDAnalysis implements an algorithm to jump to any frame in
the trajectory by initially building a list of frames and
offsets. This list is built the first time a trajectory is read and
this frame scanning process may take a while. The list of offsets is
cached and saved to disk as a hidden file that is read the next time
the trajectory file is accessed so that subsequently, loading and
random access in XTC and TRR files is very fast.
       

TODO:

TRR files may not record positions, velocities, and forces at the same
time point so MDAnalysis DOES WHAT?



Writing TRR and XTC files
~~~~~~~~~~~~~~~~~~~~~~~~~

TODO:

Anything interesting about writing these files
