.. Contains the formatted docstrings for the core modules located in 'mdanalysis/MDAnalysis/core'

**************************
Core modules
**************************

The :mod:`MDAnalysis.core` modules contain functionality essential for
MDAnalysis, such as the central data structures in
:mod:`MDAnalysis.core.universe` and :mod:`MDAnalysis.core.groups` or
the selection definitions and parsing in
:mod:`MDAnalysis.core.selection`.

.. toctree::
   :maxdepth: 1

   core/init

Important objects for users
===========================

All users of MDAnalysis need to understand the two most important
classes in this section, namely the
:class:`~MDAnalysis.core.universe.Universe` and the
:class:`~MDAnalysis.core.groups.AtomGroup`.

.. toctree::
   :maxdepth: 1

   core/universe
   core/groups


.. _topology-system-label:

Topology system
===============

The topology system is primarily of interest to developers.

.. toctree::
   :maxdepth: 1

   core/topology
   core/topologyobjects
   core/topologyattrs

.. SeeAlso:: :ref:`Developer notes for Topology
             Parsers <topology-parsers-developer-notes>`

Selection system
================

The selection system is primarily of interest to developers.

.. toctree::
   :maxdepth: 1

   core/selection

