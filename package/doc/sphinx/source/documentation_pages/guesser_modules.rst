.. Contains the formatted docstrings from the guesser modules located in 'mdanalysis/package/MDAnalysis/guesser'

**************************
guesser modules
**************************
This module contains the context-aware guessers, which are used by the :ref:`guess_TopologyAttributes <guess_TopologyAttributes>` API. Context-aware guessers' main purpose
is to be tailored guesser classes that target specific file format or force field (eg. PDB file format, or Martini forcefield).
Having such guessers makes attribute guessing more accurate and reliable than having generic guessing methods that doesn't  fit all topologies.

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

If the attribute already exist in the universe, passing the attribute of interest to the ``to_guess`` parameter will only fill the empty values of the attribute if any exists.
To override all the attribute values, you can pass the attribute to the ``force_guess`` parameter instead of the to_guess one as following::

   import MDAnalysis as mda
   from MDAnalysisTests.datafiles import PRM12
 
   u = mda.Universe(PRM12, context='default', to_guess=()) # types ['HW', 'OW', ..]

   u.guess_TopologyAttributes(force_guess=['types']) # types ['H', 'O', ..]

N.B.: If you didn't pass any ``context`` to the API, it will use the :ref:`DefaultGuesser <DefaultGuesser>`

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
