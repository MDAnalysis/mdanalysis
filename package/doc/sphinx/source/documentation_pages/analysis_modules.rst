.. Contains the formatted docstrings from the analysis modules located in 'mdanalysis/MDAnalysis/analysis', although in some cases the documentation imports functions and docstrings from other files which are also curated to reStructuredText markup.

****************
Analysis modules
****************

The :mod:`MDAnalysis.analysis` module contains code to carry out specific
analysis functionality. It is based on the core functionality (i.e. trajectory
I/O, selections etc). The analysis modules can be used as examples for how to
use MDAnalysis but also as working code for research projects; typically all
contributed code has been used by the authors in their own work.

Please see the individual module documentation for additional references and
citation information.

These modules are not imported by default; in order to use them one has to
import them from :mod:`MDAnalysis.analysis`, for instance ::

    import MDAnalysis.analysis.align

.. rubric:: Additional dependencies

Some of the modules in :mod:`MDAnalysis.analysis` require additional Python
packages to enable full functionality. For example,
:mod:`MDAnalysis.analysis.encore` provides more options if `scikit-learn`_ is
installed. These package are *not automatically installed* with
:program:`pip`(although one can add the ``[analysis]`` requirement to the
:program:`pip` command line to force their installation). If you install
MDAnalysis with :program:`conda` (see :ref:`installation-instructions`) then a
*full set of dependencies* is automatically installed.

Other modules require external programs. For instance, the
:mod:`MDAnalysis.analysis.hole` module requires an installation of the HOLE_
suite of programs. You will need to install these external dependencies by
following their installation instructions before you can use the corresponding
MDAnalysis module.

.. _scikit-learn: http://scikit-learn.org/
.. _HOLE: http://www.holeprogram.org/


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
   analysis/encore

Hydrogen bonding
================

.. toctree::
   :maxdepth: 1

   analysis/hbond_analysis
   analysis/hbond_autocorrel
   analysis/wbridge_analysis

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
   analysis/dihedrals


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

   analysis/legacy_modules
