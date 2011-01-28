.. Contains the formatted docstrings from the analysis modules located in 'mdanalysis/MDAnalysis/analysis', although in some cases the documentation imports functions and docstrings from other files which are also curated to reStructuredText markup.

**************************
Specific analysis modules
**************************

These modules are located in the :mod:`MDAnalysis.analysis` module (and the
code lives in the ``MDAnalysis/analysis/`` directory). They contain abstract
functionality for specific analysis tasks (for example, bilayer
:mod:`~MDAnalysis.analysis.leaflet` selection). In many cases it will be
preferable to use these provided tools over manually-coded alternatives because
some components may be written in C for time-critical routines.

These modules are not imported by default; in order to use them one has to ::

  import MDAnalysis.analysis

.. _scipy: http://www.scipy.org/
.. _networkx: http://networkx.lanl.gov/


.. toctree::
   :maxdepth: 1

   analysis/align
   analysis/contacts
   analysis/distances
   analysis/density
   analysis/leaflet

.. Note:: Some of the modules require additional Python packages such as
  :mod:`scipy` from the SciPy_ package or :mod:`networkx` from
  NetworkX_. These package are *not automatically installed* (although one can
  add the ``[analysis]`` requirement to the :program:`easy_install` command
  line to force their installation.

.. TODO: write a INSTALLATION page and link to it
   



