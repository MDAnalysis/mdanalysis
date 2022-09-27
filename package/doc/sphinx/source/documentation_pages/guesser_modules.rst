.. Contains the formatted docstrings from the guesser modules located in 'mdanalysis/package/MDAnalysis/guesser'

**************************
guesser modules
**************************
This module contains the context-specific guessers. Context-specific guessers main purpose
is to be a more tailored guesser classes that serve specific file formats or force fields. Having
such guessers makes attribute guessing more accurate and reliable than having one generic guessing
class that is used with all topologies.

Example uses of guessers
------------------------

Guessing using :ref:`guess_TopologyAttributes <guess_TopologyAttributes>` API
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Guessing can be done through the universe's :func:`guess_TopologyAttributes` as following::

  import MDAnalysis as mda
  from MDAnalysisTests.datafiles import PDB

  u = mda.Universe(PDB)
  print(hasattr(u.atoms, 'elements'))  # returns False
  u.guess_TopologyAttributes(to_guess=['elements'])
  print(u.atoms.elements) # print ['N' 'H' 'H' ... 'NA' 'NA' 'NA']

In the above example, we passed ``elements`` as the attribute we want to guess, and
:func:`guess_TopologyAttributes` guess then add it as a topology
attribute to the ``AtomGroup`` of the universe.

.. rubric:: available guessers
.. toctree::
   :maxdepth: 1

   guesser_modules/init
   guesser_modules/default_guesser


.. rubric:: guesser core modules

The remaining pages are primarily of interest to developers as they
contain functions and classes that are used in the implementation of
the context-specific guessers.

.. toctree::
   :maxdepth: 1

   guesser_modules/base
   guesser_modules/tables
