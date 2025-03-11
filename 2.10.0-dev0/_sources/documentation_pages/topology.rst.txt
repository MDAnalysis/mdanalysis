.. -*- coding: utf-8 -*-
.. _topology-label:

=====================
 The topology system
=====================

As shown briefly in :ref:`overview-label`, the :class:`~MDAnalysis.core.universe.Universe` class is the primary object and core interface to molecular dynamics data in MDAnalysis.
When loading topology information from a file, as with ::

  >>> from MDAnalysis import Universe
  >>> from MDAnalysis.tests.datafiles import PSF
  >>> u = Universe(PSF)

the file is read, the contents parsed, and a :class:`~MDAnalysis.core.topology.Topology` object is constructed from these contents.
