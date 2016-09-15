============================================================
 :mod:`MDAnalysis.analysis.legacy` --- Legacy analysis code
============================================================

.. versionadded:: 0.16.0

The :mod:`MDAnalysis.analysis.legacy` package contains analysis
modules that are not or only incompletely tested and not regularly
maintained. They nevertheless still provide useful and sometimes
unique analysis capabilities and are therefore provided **as is**.

.. warning::

   Code in this module is not regularly maintained. Please use it very
   carefully.

If you want to use code from this module then you will have to import
it explicitly. For example, ::

   from MDAnalysis.analysis.legacy import x3dna

(For further discussion, see `Issue 743`_.)


.. _Issue 743: https://github.com/MDAnalysis/mdanalysis/issues/743


Legacy modules
==============

.. toctree::
   :maxdepth: 1

   legacy/x3dna


