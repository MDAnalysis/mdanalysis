
.. -*- coding: utf-8 -*-
.. _selection-exporters:

===================
Selection exporters
===================

Selection exporters allow you to write a selection of atoms to a file that can be read by another program. 

.. table:: Supported selection exporters

    .. include:: selection_exporter_formats.txt

Writing selections
==================

Single AtomGroup
----------------

The typical situation is that one has an
:class:`~MDAnalysis.core.groups.AtomGroup` and wants to work with the
same selection of atoms in a different package, for example, to
visualize the atoms in VMD_. 

.. ipython:: python

    import MDAnalysis as mda
    from MDAnalysis.tests.datafiles import PDB
    u = mda.Universe(PDB)
    ag = u.select_atoms('resname ALA')

As with a normal structure file, use :meth:`AtomGroup.write <MDAnalysis.core.groups.AtomGroup.write>` method with the appropriate file extension.

.. code-block :: python

    ag.write("ala_selection.vmd", name="alanine")

In VMD_, sourcing the file ``ala_selection.vmd`` (written in Tcl) defines the "macro" ``alanine`` that contains the atom indices to select.

.. code-block:: tcl

    source ala_selection.vmd
    set sel [atomselect top alanine]

and in the GUI the macro appears in the :menuselection:`Graphics -->
Representations` window in the list :guilabel:`Selections: Singlewords` as
"alanine".

Names are not always required; if ``name`` is not passed to  :meth:`AtomGroup.write <MDAnalysis.core.groups.AtomGroup.write>`, MDAnalysis defaults to "mdanalysis001", "mdanalysis002", and so on.

Multiple selections
-------------------

:meth:`AtomGroup.write <MDAnalysis.core.groups.AtomGroup.write>` can take additional keyword arguments, including ``mode``. The default is ``mode='w'``, which will overwrite the provided filename. If ``mode='a'``, the selection is appended to the file.

.. code-block :: python

    u.select_atoms('resname T*').write('residues.ndx',
                                       name='TYR_THR',
                                       mode='a')
    u.select_atoms('resname GLY').write('residues.ndx', 
                                        name='GLY', 
                                        mode='a')
    u.select_atoms('resname PRO').write('residues.ndx', 
                                        name='PRO', 
                                        mode='a')

Looking at this GROMACS index file, we see:

.. code-block:: console

    $ gmx make_ndx -n residues.ndx

    Command line:
    gmx make_ndx -n residues.ndx

    Going to read 1 old index file(s)
    Counted atom numbers up to 3341 in index file

    0 TYR_THR             :   301 atoms
    1 GLY                 :   141 atoms
    2 PRO                 :   140 atoms

    nr : group      '!': not  'name' nr name   'splitch' nr    Enter: list groups
    'a': atom       '&': and  'del' nr         'splitres' nr   'l': list residues
    't': atom type  '|': or   'keep' nr        'splitat' nr    'h': help
    'r': residue              'res' nr         'chain' char
    "name": group             'case': case sensitive           'q': save and quit
    'ri': residue index

Alternatively, you can direcly use the selection writer itself as a `context manager`_ and write each :class:`~MDAnalysis.core.groups.AtomGroup` inside the context. For example:

.. code-block :: python

    with mda.selections.gromacs.SelectionWriter('residues.ndx', mode='w') as ndx:
        ndx.write(u.select_atoms('resname T*'), 
                  name='TYR_THR')
        ndx.write(u.select_atoms('resname GLY'),
                  name='GLY')

And again, you can append to the file with ``mode='a'``:

.. code-block :: python

    with mda.selections.gromacs.SelectionWriter('residues.ndx', mode='a') as ndx:
        ndx.write(u.select_atoms('resname PRO'), 
                  name='PRO')


Reading in selections
=====================

Currently, MDAnalysis doesn't support reading in atom selections. However, there are other tools that can read files from other programs, such as `GromacsWrapper`_. 

.. _CHARMM: http://www.charmm.org
.. _Gromacs: http://www.gromacs.org
.. _VMD: http://www.ks.uiuc.edu/Research/vmd/
.. _PyMol: http://www.pymol.org
.. _Jmol: http://wiki.jmol.org/
.. _GromacsWrapper: https://gromacswrapper.readthedocs.io/en/latest/
.. _`context manager`: https://docs.python.org/3/reference/datamodel.html#context-managers
