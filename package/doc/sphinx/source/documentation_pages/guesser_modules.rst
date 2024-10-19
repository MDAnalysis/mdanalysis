.. Contains the formatted docstrings from the guesser modules located in 'mdanalysis/package/MDAnalysis/guesser'

.. _Guessers:

**************************
Guesser modules
**************************
This module contains the context-aware guessers, which are used by the :meth:`~MDAnalysis.core.universe.Universe.guess_TopologyAttrs` API. Context-aware guessers' main purpose
is to be tailored guesser classes that target specific file format or force field (eg. PDB file format, or Martini forcefield).
Having such guessers makes attribute guessing more accurate and reliable than having generic guessing methods that don't fit all scenarios.

Example uses of guessers
------------------------

Default behavior
~~~~~~~~~~~~~~~~

By default, MDAnalysis will guess the "mass" and "type" (atom type) attributes for all particles in the Universe
using the :class:`~MDAnalysis.guesser.default_guesser.DefaultGuesser` at the time of Universe creation,
if they are not read from the input file.
Please see the :ref:`DefaultGuesser`_ for more information.



Guessing using :meth:`~MDAnalysis.core.universe.Universe.guess_TopologyAttrs` API
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Guessing can be done through the Universe's :meth:`~MDAnalysis.core.universe.Universe.guess_TopologyAttrs` as following::

  import MDAnalysis as mda
  from MDAnalysisTests.datafiles import PDB

  u = mda.Universe(PDB)
  print(hasattr(u.atoms, 'elements'))  # returns False
  u.guess_TopologyAttrs(to_guess=['elements'])
  print(u.atoms.elements) # print ['N' 'H' 'H' ... 'NA' 'NA' 'NA']

In the above example, we passed ``elements`` as the attribute we want to guess
:meth:`~MDAnalysis.core.universe.Universe.guess_TopologyAttrs` guess then add it as a topology
attribute to the ``AtomGroup`` of the universe.

If the attribute already exists in the universe, passing the attribute of interest to the ``to_guess`` parameter will only fill the empty values of the attribute if any exists.
To override all the attribute values, you can pass the attribute to the ``force_guess`` parameter instead of ``to_guess`` as following::

   import MDAnalysis as mda
   from MDAnalysisTests.datafiles import PRM12
 
   u = mda.Universe(PRM12, context='default', to_guess=()) # types ['HW', 'OW', ..]

   u.guess_TopologyAttrs(force_guess=['types']) # types ['H', 'O', ..]


.. note::
   The default ``context`` will use the :class:`~MDAnalysis.guesser.default_guesser.DefaultGuesser`




.. rubric:: Available guessers
.. toctree::
   :maxdepth: 1

   guesser_modules/init
   guesser_modules/default_guesser


.. rubric:: Guesser core modules

The remaining pages are primarily of interest to developers as they
contain functions and classes that are used in the implementation of
the context-specific guessers.

.. toctree::
   :maxdepth: 1

   guesser_modules/base
   guesser_modules/tables
