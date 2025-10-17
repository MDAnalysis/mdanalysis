.. Contains the formatted docstrings for the core modules located in 'mdanalysis/MDAnalysis/core'

.. _`Selection exporters`:

*******************
Selection exporters
*******************

The classes in this module allow writing a selection to a file that can be read by
*another program* to select the atoms in a MDAnalysis
:class:`MDAnalysis.core.groups.AtomGroup`. Such cross-package interoperability
allows a user to combine their favourite tools with MDAnalysis for further
visualization or simulation.


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

How to write selections
=======================

Single AtomGroup
----------------

The typical situation is that one has an
:class:`~MDAnalysis.core.groups.AtomGroup` and wants to work with the
same selection of atoms in a different package, for example, to
visualize the atoms in VMD_. First create an :class:`AtomGroup` (named
``g`` in the example below) and then use its
:class:`~MDAnalysis.core.groups.AtomGroup.write` method with the
appropriate *file extension* (see :ref:`Supported selection
exporters` for the recognized *extension*)::

  g = u.select_atoms('protein")
  g.write("selections.vmd", name="mda_protein")

In VMD_, sourcing the file ``selections.vmd`` (written in Tcl) defines
the "macro" ``mda_protein`` that contains the atom indices to select

.. code-block:: tcl

   source selection.vmd
   set sel [atomselect top mda_mdanalysis]

and in the GUI the macro appears in the :menuselection:`Graphics -->
Representations` window in the list :guilabel:`Selections: Singlewords` as
"mda_protein".


Multiple selections
-------------------

The :class:`~MDAnalysis.core.groups.AtomGroup.write` method can take
additional keyword arguments, including ``mode``. The default is
``mode="w"``, which will overwrite the provided file. If ``mode="a"``
then the selection is *appended* to the file. 

Alternatively, one may use the
:class:`~MDAnalysis.selections.base.SelectionWriter` itself as a
context manager and write each :class:`AtomGroup` inside the
context. For example, to write multiple groups that were selected to
mark different parts of a lipid bilayer to Gromacs_ index file named
"leaflets.ndx"::

   with mda.selections.gromacs.SelectionWriter('leaflets.ndx', mode='w') as ndx:
       ndx.write(upper_saturated, name='upper_sat')
       ndx.write(lower_saturated, name='lower_sat')
       ndx.write(upper_unsaturated, name='upper_unsat')
       ndx.write(lower_unsaturated, name='lower_unsat')

There is a separate :class:`SelectionWriter` for each format, as
described next.


Selection writers
=================

There exist different :class:`~MDAnalysis.selections.base.SelectionWriterBase`
classes for different packages. The
:func:`~MDAnalysis.selections.get_writer` function can automatically pick
the appropriate one, based on the file name extension in the :ref:`Supported
selection exporters`.

.. autofunction:: MDAnalysis.selections.get_writer


.. rubric:: Formats

Each module implements a :class:`SelectionWriter` for a specific format.

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
.. _Jmol: https://jmol.sourceforge.net/
