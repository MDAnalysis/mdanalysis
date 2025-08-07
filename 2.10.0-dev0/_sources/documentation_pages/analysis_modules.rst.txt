.. Contains the formatted docstrings from the analysis modules located in 'mdanalysis/MDAnalysis/analysis', although in some cases the documentation imports functions and docstrings from other files which are also curated to reStructuredText markup.

.. module:: MDAnalysis.analysis

.. _analysis-label:

****************
Analysis modules
****************

The :mod:`MDAnalysis.analysis` module provides a wide collection of analysis tools for 
molecular dynamics trajectories. These modules build upon MDAnalysis core functionality 
(trajectory I/O, selections, etc.) and are designed both for reuse in research workflows 
and as examples of using the MDAnalysis API. Each module typically defines an analysis 
class that follows a standard interface.

See the `User Guide Analysis section`_ for interactive examples and additional context.

.. _User Guide Analysis section: https://userguide.mdanalysis.org/stable/examples/analysis/README.html   


Getting started with analysis
=============================

General usage pattern
---------------------

Most analysis tools are implemented as single classes and follow this usage pattern:

#. Import the module (e.g., :mod:`MDAnalysis.analysis.rms`).
#. Initialize the analysis class with the required arguments.
#. Run the analysis with :meth:`~MDAnalysis.analysis.base.AnalysisBase.run`.
#. Access results via the :attr:`~MDAnalysis.analysis.base.AnalysisBase.results` attribute.

.. code-block:: python

   from MDAnalysis.analysis import ExampleAnalysisModule  # (e.g. RMSD)
   analysis_obj = ExampleAnalysisModule.AnalysisClass(universe, ...)
   analysis_obj.run(start=start_frame, stop=stop_frame, step=step)
   print(analysis_obj.results)

Please see the individual module documentation for any specific caveats 
and also read and cite the reference papers associated with these algorithms.


Using parallelization for analysis tools
----------------------------------------

.. versionadded:: 2.8.0

Many analysis tools (based on :class:`~MDAnalysis.analysis.base.AnalysisBase`)
can be :ref:`run in parallel <parallel-analysis>` using a simple
split-apply-combine scheme whereby slices of the trajectory ("split") are analyzed in
parallel ("apply" the analysis function) and the data from the parallel executions
are "combined" at the end.

MDAnalysis supports different :ref:`backends <backends>` for the parallel execution such as
:mod:`multiprocessing` or `dask`_ (see :mod:`MDAnalysis.analysis.backends`).
As a special case, serial execution is handled by the default  ``backend='serial'``, i.e.,
by default, none of the analysis tools run in parallel and one has to explicitly request 
parallel execution. Without any additionally installed dependencies, only one parallel backend
is supported -- Python :mod:`multiprocessing` (which is available in the Python standard 
library), which processes each slice of a trajectory by running a separate *process* on a 
different core of a multi-core CPU.

.. _dask: https://dask.org/

.. Note::

   Not all analysis tools in MDAnalysis can be parallelized and others have 
   not yet been updated to make use of the :ref:`parallelization framework <parallel-analysis>`,
   which was introduced in release 2.8.0. MDAnalysis aims to have parallelization enabled for
   all analysis tools that support it by release 3.0.

In order to use parallelization, add ``backend='multiprocessing'`` to the arguments of the
:meth:`~MDAnalysis.analysis.base.AnalysisBase.run` method together with  ``n_workers=N`` where
``N`` is the number of CPUs that you want to use for parallelization. 
(You can use ``multiprocessing.cpu_count()`` to get the maximum available number of CPUs on your 
machine but this may not always lead to the best performance because of computational overheads and
the fact that parallel access to a single trajectory file is often a performance bottleneck.) As an
example we show how to run an RMSD calculation in parallel:

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

Be explicit and specify both ``backend`` and ``n_workers``. Choosing too many 
workers or using large trajectory frames may lead to an out-of-memory error.

You can also implement your own backends -- see :mod:`MDAnalysis.analysis.backends`.

.. SeeAlso::
   :ref:`parallel-analysis` for technical details
   


Additional dependencies
-----------------------

Some of the modules in :mod:`MDAnalysis.analysis` require additional Python
packages to enable full functionality. For example,
:mod:`MDAnalysis.analysis.encore` provides more options if `scikit-learn`_ is
installed. If you installed MDAnalysis with :program:`pip` (see
:ref:`installation-instructions`) these packages are *not automatically
installed* although one can add the ``[analysis]`` tag to the :program:`pip`
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
   analysis/dssp

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
