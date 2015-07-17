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

.. table:: Table of supported exporters and recognized file name extensions.

   +---------------+-----------+-------+--------------------------------------------+
   |Name           | extension |  IO   | remarks                                    |
   +===============+===========+=======+============================================+
   | vmd           | tcl       | w     | VMD_ macros, available in Representations; |
   |               |           |       | module :mod:`MDAnalysis.selections.vmd`    |
   +---------------+-----------+-------+--------------------------------------------+
   | pymol         | pml       | w     | simple PyMOL_ selection string;            |
   |               | 	       |       | module :mod:`MDAnalysis.selections.pymol`  |
   +---------------+-----------+-------+--------------------------------------------+
   | gromacs       | ndx       | w     | Gromacs_ index file;                       |
   |               |           |       | module :mod:`MDAnalysis.selections.gromacs`|
   +---------------+-----------+-------+--------------------------------------------+
   | charmm        | str       | w     | CHARMM_ selection of individual atoms;     |
   |               |           |       | module :mod:`MDAnalysis.selections.charmm` |
   +---------------+-----------+-------+--------------------------------------------+
   | jmol          | spt       | w     | Jmol_ selection commands;                  |
   |               |           |       | module :mod:`MDAnalysis.selections.jmol`   |
   +---------------+-----------+-------+--------------------------------------------+

.. rubric:: Contents

.. toctree::
   :maxdepth: 1

   selections/vmd
   selections/pymol
   selections/gromacs
   selections/charmm
   selections/jmol
   selections/base

.. _CHARMM: http://www.charmm.org
.. _Gromacs: http://www.gromacs.org
.. _VMD: http://www.ks.uiuc.edu/Research/vmd/
.. _PyMOL: http://www.pymol.org
.. _Jmol: http://wiki.jmol.org/
