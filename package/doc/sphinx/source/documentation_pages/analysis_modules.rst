.. Contains the formatted docstrings from the analysis modules located in 'mdanalysis/MDAnalysis/analysis', although in some cases the documentation imports functions and docstrings from other files which are also curated to reStructuredText markup.

****************
Analysis modules
****************

The :mod:`MDAnalysis.analysis` module contains code to carry out specific
analysis functionality for MD trajectories.  It is based on the core
functionality (i.e. trajectory I/O, selections etc). The analysis modules can
be used as examples for how to use MDAnalysis but also as working code for
research projects; typically all contributed code has been used by the authors
in their own work.

Getting started with analysis
=============================

.. SeeAlso::

   The `User Guide: Analysis`_ contains extensive documentation of the analysis
   capabilities with user-friendly examples.

.. _`User Guide: Analysis`:
   https://userguide.mdanalysis.org/stable/examples/analysis/README.html   

Using the analysis classes
--------------------------

Most analysis tools in MDAnalysis are written as a single class. An analysis
usually follows the same pattern:

#. Import the desired module, since analysis modules are not imported 
   by default.
#. Initialize the analysis class instance from the previously imported module.   
#. Run the analysis, optionally for specific trajectory slices.
#. Access the analysis from the :attr:`results` attribute

.. code-block:: python

   from MDAnalysis.analysis import ExampleAnalysisModule  # (e.g. RMSD)

   analysis_obj = ExampleAnalysisModule.AnalysisClass(universe, ...)
   analysis_obj.run(start=start_frame, stop=stop_frame, step=step)
   print(analysis_obj.results)


Please see the individual module documentation for any specific caveats 
and also read and cite the reference papers associated with these algorithms.


Using parallelization for built-in analysis runs
------------------------------------------------

.. versionadded:: 2.8.0

:class:`~MDAnalysis.analysis.base.AnalysisBase` subclasses can run on a backend
that supports parallelization (see :mod:`MDAnalysis.analysis.backends`). All
analysis runs use ``backend='serial'`` by default, i.e., they do not use
parallelization by default, which has been standard before release 2.8.0 of
MDAnalysis.

Without any dependencies, only one backend is supported -- built-in
:mod:`multiprocessing`, that processes parts of a trajectory running separate
*processes*, i.e. utilizing multi-core processors properly.

.. Note::

   For now, parallelization has only been added to
   :class:`MDAnalysis.analysis.RMS.RMSD`, but by release 3.0 version it will be
   introduced to all subclasses that can support it.

In order to use that feature, simply add ``backend='multiprocessing'`` to your
run, and supply it with proper ``n_workers`` (use ``multiprocessing.cpu_count()``
for maximum available on your machine):

.. code-block:: python

   import multiprocessing
   import MDAnalysis as mda
   from MDAnalysisTests.datafiles import PSF, DCD
   from MDAnalysis.analysis.rms import RMSD
   from MDAnalysis.analysis.align import AverageStructure

   # initialize the universe
   u = mda.Universe(PSF, DCD)

   # calculate average structure for reference
   avg = AverageStructure(mobile=u).run()
   ref = avg.results.universe

   # initialize RMSD run
   rmsd = RMSD(u, ref, select='backbone')
   rmsd.run(backend='multiprocessing', n_workers=multiprocessing.cpu_count())

For now, you have to be explicit and specify both ``backend`` and ``n_workers``,
since the feature is new and there are no good defaults for it. For example,
if you specify a too big `n_workers`, and your trajectory frames are big,
you might get and out-of-memory error when executing your run.

You can also implement your own backends -- see :mod:`MDAnalysis.analysis.backends`.


Additional dependencies
-----------------------

Some of the modules in :mod:`MDAnalysis.analysis` require additional Python
packages to enable full functionality. For example,
:mod:`MDAnalysis.analysis.encore` provides more options if `scikit-learn`_ is
installed. If you installed MDAnalysis with :program:`pip` (see
:ref:`installation-instructions`) these packages are *not automatically
installed*.  Although, one can add the ``[analysis]`` tag to the :program:`pip`
command to force their installation. If you installed MDAnalysis with
:program:`conda` then a *full set of dependencies* is automatically installed.

Other modules require external programs. For instance, the
:mod:`MDAnalysis.analysis.hole2` module requires an installation of the HOLE_
suite of programs. You will need to install these external dependencies by
following their installation instructions before you can use the corresponding
MDAnalysis module.

.. _scikit-learn: http://scikit-learn.org/
.. _HOLE: http://www.holeprogram.org/


Building blocks for Analysis
============================

The building block for the analysis modules is
:class:`MDAnalysis.analysis.base.AnalysisBase`.
To build your own analysis class start by reading the documentation.

.. toctree::
   :maxdepth: 1

   analysis/base
   analysis/backends
   analysis/results
   analysis/parallelization

Distances and contacts
======================

.. toctree::
   :maxdepth: 1

   analysis/align
   analysis/contacts
   analysis/distances
   analysis/atomicdistances
   analysis/rms
   analysis/psa
   analysis/encore
   analysis/bat

Hydrogen bonding
================

.. toctree::
   :maxdepth: 1

   analysis/hydrogenbonds
   analysis/hbond_autocorrel
   analysis/wbridge_analysis

Deprecated modules:

.. toctree::
   :maxdepth: 1

   analysis/hbond_autocorrel_deprecated	      

Membranes and membrane proteins
===============================

.. toctree::
   :maxdepth: 1

   analysis/hole2
   analysis/leaflet

Nucleic acids
=============

.. toctree::
   :maxdepth: 1

   analysis/nuclinfo
   analysis/nucleicacids

Polymers
========

.. toctree::
   :maxdepth: 1

   analysis/polymer


Structure
=========

Macromolecules
--------------

.. toctree::
   :maxdepth: 1

   analysis/gnm
   analysis/helix_analysis
   analysis/dihedrals

Liquids
-------

.. toctree::
   :maxdepth: 1

   analysis/rdf
   analysis/msd

Volumetric analysis
===================

.. toctree::
   :maxdepth: 1

   analysis/density
   analysis/lineardensity
   analysis/waterdynamics
   analysis/dielectric

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

Data
====

.. toctree::
   :maxdepth: 1

   analysis/data
