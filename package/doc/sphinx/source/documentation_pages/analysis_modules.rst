.. Contains the formatted docstrings from the analysis modules located in 'mdanalysis/MDAnalysis/analysis', although in some cases the documentation imports functions and docstrings from other files which are also curated to reStructuredText markup.

****************
Analysis modules
****************

The :mod:`MDAnalysis.analysis` module contains code to carry out
specific analysis functionality. It is based on the core functionality
(i.e. trajectory I/O, selections etc). The analysis modules can be
used as examples for how to use MDAnalysis but also as working code
for research projects; typically all contributed code has been used by
the authors in their own work.

Please see the individual module documentation for additional
references and citation information.

These modules are not imported by default; in order to use them one
has to import them from :mod:`MDAnalysis.analysis`, for instance ::

    import MDAnalysis.analysis.align

.. Note::

  Some of the modules require additional Python packages such as :mod:`scipy`
  from the SciPy_ package. These package are *not automatically installed*
  (although one can add the ``[analysis]`` requirement to the :program:`pip`
  command line to force their installation.

.. _scipy: http://www.scipy.org/


Building blocks for Analysis
============================

.. toctree::
   :maxdepth: 1

   analysis/base


Distances and contacts
======================

.. toctree::
   :maxdepth: 1

   analysis/align
   analysis/contacts
   analysis/distances
   analysis/rms
   analysis/psa

Hydrogen bonding
================

.. toctree::
   :maxdepth: 1

   analysis/hbond_analysis
   analysis/hbond_autocorrel

Membranes and membrane proteins
===============================

.. toctree::
   :maxdepth: 1

   analysis/hole
   analysis/leaflet

Nucleic acids
=============

.. toctree::
   :maxdepth: 1

   analysis/nuclinfo
   analysis/legacy/x3dna

Polymers
========

.. toctree::
   :maxdepth: 1

   analysis/polymer


Structure
=========

.. toctree::
   :maxdepth: 1

   analysis/gnm
   analysis/helanal
   analysis/rdf


Volumetric analysis
===================

.. toctree::
   :maxdepth: 1

   analysis/density
   analysis/lineardensity
   analysis/waterdynamics

Dimensionality Reduction
========================
.. toctree::
   :maxdepth: 1

   analysis/diffusionmap
   analysis/pca


Legacy analysis modules
=======================

The :mod:`MDAnalysis.analysis.legacy` module contains code that for a
range of reasons is not as well maintained and tested as the other
analysis modules. *Use with care.*

.. toctree::
   :maxdepth: 1

   analysis/legacy_module
