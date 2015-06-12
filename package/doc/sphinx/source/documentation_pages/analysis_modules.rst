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

  Some of the modules require additional Python packages such as
  :mod:`scipy` from the SciPy_ package or :mod:`networkx` from
  NetworkX_. These package are *not automatically installed* (although
  one can add the ``[analysis]`` requirement to the
  :program:`easy_install` command line to force their installation.

.. _scipy: http://www.scipy.org/
.. _networkx: http://networkx.lanl.gov/



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
   analysis/x3dna


Structure
=========

.. toctree::
   :maxdepth: 1

   analysis/gnm
   analysis/helanal


Volumetric analysis
===================   

.. toctree::
   :maxdepth: 1

   analysis/density
   analysis/waterdynamics


