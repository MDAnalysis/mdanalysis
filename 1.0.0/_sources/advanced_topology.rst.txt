.. -*- coding: utf-8 -*-
.. _advanced-topology:

==========================
Advanced topology concepts
==========================

.. _adding-residue:

Adding a Residue or Segment to a Universe
=========================================

To add a :class:`~MDAnalysis.core.groups.Residue` or :class:`~MDAnalysis.core.groups.Segment` to a topology, use the :meth:`Universe.add_Residue <MDAnalysis.core.universe.Universe.add_Residue>` or :meth:`Universe.add_Segment <MDAnalysis.core.universe.Universe.add_Segment>` methods.

.. code-block::

    >>> u = mda.Universe(PSF, DCD)
    >>> u.segments
    <SegmentGroup with 1 segment>
    >>> u.segments.segids
    array(['4AKE'], dtype=object)
    >>> newseg = u.add_Segment(segid='X')
    >>> u.segments.segids
    array(['4AKE', 'X'], dtype=object)
    >>> newseg.atoms
    <AtomGroup with 0 atoms>

To assign the last 100 residues from the :class:`~MDAnalysis.core.universe.Universe` to this new Segment:

.. code-block::

    >>> u.residues[-100:].segments = newseg
    >>> newseg.atoms
    <AtomGroup with 1600 atoms>

Another example is `creating custom segments for protein domains <examples/constructing_universe.ipynb#Adding-a-new-segment>`_.


.. _molecule:

Molecules
=========

In MDAnalysis, a molecule is a GROMACS-only concept that is relevant in some analysis methods. A group of atoms is considered a "molecule" if it is defined by the :code:`[ moleculetype ]` section in a `GROMACS topology <http://manual.gromacs.org/documentation/2019/reference-manual/file-formats.html#top>`_. Molecules are only defined if a Universe is created from a GROMACS topology file (i.e. with a .tpr extension). Unlike fragments, they are not accessible directly from atoms.

.. code-block::

    >>> tpr = mda.Universe(TPR)
    >>> tpr.atoms.molecules
    Traceback (most recent call last):
    File "<stdin>", line 1, in <module>
    File "MDAnalysis/core/groups.py", line 2278, in __getattr__
        cls=self.__class__.__name__, attr=attr))
    AttributeError: AtomGroup has no attribute molecules

However, the order (:code:`molnum`) and name (:code:`moltype`) of each molecule is accessible as :ref:`topology attributes <topology-attributes>`::

    >>> tpr.atoms.molnums
    array([    0,     0,     0, ..., 11086, 11087, 11088])
    >>> tpr.atoms.moltypes
    array(['AKeco', 'AKeco', 'AKeco', ..., 'NA+', 'NA+', 'NA+'], dtype=object)


.. _custom-topologyattrs:

Adding custom TopologyAttrs
===========================