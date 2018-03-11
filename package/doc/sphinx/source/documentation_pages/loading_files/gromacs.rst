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
   `Loading XYZ flies <load_xyz>`
   

.. _load_tpr:

TPR files
---------

MDAnalysis supports reading TPR files created by version 3.3 of Gromacs
onwards.

.. _load_gro:


Reading GRO files
-----------------

Stuff about the GRO format specifically

When used as a topology file in a Universe, MDAnalysis will read
``names``, ``ids``, ``resids`` and ``resnames``,
and guess ``masses`` and ``types`` based on the names of atoms.
The ``segid`` for of all atoms is set to "``SYSTEM``".

For implementation details see
:mod:`MDAnalysis.topology.GROParser`.


Writing GRO files
^^^^^^^^^^^^^^^^^

MDAnalysis can also write GRO files, for example using the
:meth:`~MDAnalysis.core.groups.AtomGroup.write` method.
It will attempt to use the ``name``, ``resname`` and ``resid`` attribute
from atoms if available, otherwise default values of "``X``", "``UNK``"
and "``1``" respectively will be used.
If the provided atoms have velocities these will also be written,
otherwise only the positions will be written.
The precision is limited to three decimal places.

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


Writing TRR and XTC files
^^^^^^^^^^^^^^^^^^^^^^^^^

Anything interesting about writing these files
