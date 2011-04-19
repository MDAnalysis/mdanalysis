.. Contains the formatted docstrings for the core modules located in 'mdanalysis/MDAnalysis/core'

.. _`Selection exporters`:

*******************
Selection exporters
*******************

The classes in this module allow writing a selection to a file that can be read by
*another programme* to select the atoms in a MDAnalysis
:class:`MDAnalysis.core.AtomGroup.AtomGroup`. Such cross-package interoperability
allows a user to combine their favourite tools with MDAnalysis for further
visualization or simulation.

There exist different :class:`~MDAnalysis.selections.base.SelectionWriter` classes
for different package. The :func:`~MDAnalysis.selections.base.get_writer` function
can automatically pick the appropriate one, based on the file name extension in the
table below.

.. _Supported selection exporters:

Supported exporters and recognized file name extensions:

+---------------+-----------+-------+--------------------------------------------+
|Name           | extension |  IO   | remarks                                    |
+===============+===========+=======+============================================+
| VMD_          | tcl       | w     | VMD macros, available in Representations;  |
|               |           |       | module :mod:`MDAnalysis.selections.vmd`    |
+---------------+-----------+-------+--------------------------------------------+
| PyMOL_        | pml       | w     | simple selection string;                   |
|               | 	    | 	    | module :mod:`MDAnalysis.selections.pymol`  |
+---------------+-----------+-------+--------------------------------------------+
| Gromacs_      | ndx       | w     | index file;                                |
|               |           |       | module :mod:`MDAnalysis.selections.gromacs`|
+---------------+-----------+-------+--------------------------------------------+
| CHARMM_       | str       | w     | selection of individual atoms;             |
|               |           |       | module :mod:`MDAnalysis.selections.charmm` |
+---------------+-----------+-------+--------------------------------------------+


.. toctree::
   :maxdepth: 1

   selections/vmd
   selections/pymol
   selections/gromacs
   selections/charmm
   selections/base

.. _CHARMM: http://www.charmm.org
.. _Gromacs: http://www.gromacs.org
.. _VMD: http://www.ks.uiuc.edu/Research/vmd/
.. _PyMOL: http://www.pymol.org
