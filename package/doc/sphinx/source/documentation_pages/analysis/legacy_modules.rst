============================================================
 :mod:`MDAnalysis.analysis.legacy` --- Legacy analysis code
============================================================

.. versionadded:: 0.16.0

The :mod:`MDAnalysis.analysis.legacy` package contains analysis
modules that are not or only incompletely tested and not regularly
maintained. They nevertheless still provide useful and sometimes
unique analysis capabilities and are therefore provided **as
is**. (For further discussion, see `Issue 743`_.)

.. warning::

   Code in the :mod:`~MDAnalysis.analysis.legacy` package is not
   regularly maintained. Please use it very carefully.

If you want to use modules from this package then you will have to import
them explicitly. For example, ::

   from MDAnalysis.analysis.legacy import x3dna


.. _Issue 743: https://github.com/MDAnalysis/mdanalysis/issues/743


Legacy modules
==============

.. toctree::
   :maxdepth: 1

   legacy/x3dna


