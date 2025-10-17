.. _overview-label:

==========================
 Overview of MDAnalysis
==========================

MDAnalysis is a Python package for the analysis of molecular dynamics simulations. 
It provides an object-oriented interface to molecular structures and trajectories, 
with direct access to atomic coordinates as :class:`numpy.ndarray` objects for seamless 
integration with `NumPy`_ and `SciPy`_.

This page gives a high-level overview of the most important public classes and modules. 
For usage examples and tutorials, refer to the `User Guide`_.

.. _User Guide: https://userguide.mdanalysis.org/stable/index.html
.. _NumPy: https://numpy.org/
.. _SciPy: https://scipy.org/

Key Classes
===========

The core of MDAnalysis revolves around the :class:`~MDAnalysis.core.universe.Universe` class, 
which serves as the central data structure that loads and connects topology and coordinate data.
From a :class:`~MDAnalysis.core.universe.Universe`, users typically interact with :class:`~MDAnalysis.core.groups.AtomGroup` 
objects — flexible collections of atoms that support structural selections and analysis operations. These selections 
are created using `CHARMM-style`_ selection syntax via the :meth:`~MDAnalysis.core.groups.AtomGroup.select_atoms` method, 
allowing users to query atoms based on names, residue numbers, segments, and more.

Individual atoms are represented by the :class:`~MDAnalysis.core.groups.Atom` class, while residues and segments (or chains) are modeled using the 
:class:`~MDAnalysis.core.groups.Residue` and :class:`~MDAnalysis.core.groups.Segment` classes, respectively. Together, these 
classes form an intuitive, object-oriented hierarchy that makes it easy to navigate and analyze molecular systems.

.. _CHARMM-style: http://www.charmm.org/documentation/c37b1/select.html

Core modules
============

MDAnalysis is organized into several core modules that provide specialized functionality for 
handling and analyzing molecular dynamics data. The :mod:`MDAnalysis.core` module defines the 
essential data structures such as :class:`~MDAnalysis.core.universe.Universe`, :class:`~MDAnalysis.core.groups.AtomGroup`, 
and related objects. The :mod:`MDAnalysis.analysis` module contains a collection of analysis tools for tasks like RMSD calculation, 
diffusion analysis, contact maps, and more. The :mod:`MDAnalysis.selections` module implements the flexible selection language used 
to query atoms based on structural properties. Finally, :mod:`MDAnalysis.topology` manages topology parsing and representation, 
supporting a wide range of file formats for loading molecular structures.




