================================
  MDAnalysis Repository README
================================

|numfocus| |build| |cov| [*]_

|docs| |devdocs| |usergroup| |developergroup| |anaconda| |mybinder|

MDAnalysis_ is a Python library for the analysis of computer simulations of many-body systems at the molecular scale, spanning use cases from interactions of drugs with proteins to novel materials. It is widely used in the scientific community and is written by scientists for scientists. 

It works with a a wide range of popular simulation packages including Gromacs, Amber, NAMD, CHARMM, DL_Poly, HooMD, LAMMPS and many other — see the lists of  supported `trajectory formats`_ and `topology formats`_.
MDAnalysis also includes widely used analysis algorithms in the `MDAnalysis.analysis`_ module.

.. _numfocus-fiscal-sponsor-attribution:

The MDAnalysis project uses an `open governance model`_ and is fiscally sponsored by `NumFOCUS`_. Consider making 
a `tax-deductible donation`_ to help the project pay for developer time, professional services, travel, workshops, and a variety of other needs.

.. image:: https://raw.githubusercontent.com/numfocus/templates/master/images/numfocus-logo.png
  :height: 60px
  :target: https://numfocus.org/project/mdanalysis
  :align: center
  :alt: NumFOCUS
  
This project is bound by a `Code of Conduct`_.


Example analysis script
=======================

.. code:: python

   import MDAnalysis as mda

   # Load simulation results with a single line
   u = mda.Universe('topol.tpr','traj.trr')

   # Select atoms
   ag = u.select_atoms('name OH')

   # Atom data made available as Numpy arrays
   ag.positions
   ag.velocities
   ag.forces

   # Iterate through trajectories
   for ts in u.trajectory:
       print(ag.center_of_mass())
 
There are a number of tutorials_ on the MDAnalysis homepage that explain
how to conduct RMSD calculations, Alignment and many more features of MDAnalysis.

Source code
===========

Source code is hosted in a git repository at

https://github.com/MDAnalysis/mdanalysis

and is available under the GNU General Public License, version 2 (see
the file LICENSE_).

This is the top level of the master repository. It contains

1. the MDAnalysis toolkit source files in the directory ::

      package/

2. the unit tests together with any input files required for
   running those tests in the directory ::

      testsuite/

The directory ``maintainer`` contains scripts only needed for
maintaining releases and are not generally useful for the user or the
typical developer.

(For more details on the directory layout see `Issue 87`_ on the
MDAnalysis issue tracker.)

Guide for Developers
====================

To setup a development environment and run the testsuite you can use this
guide_. If you are a new developer who would like to start contributing to
MDAnalysis as a start you can increase our code coverage, the guides explain how
to find uncovered code.



.. Footnotes

.. [*] **build**: Unit testing is for the whole package; **coverage** is
       shown for the core library modules and the analysis modules.

.. _NumFOCUS: https://numfocus.org/
.. _open governance model: https://www.mdanalysis.org/about/#governance
.. _tax-deductible donation: https://numfocus.org/donate-to-mdanalysis
.. _`Code of Conduct`: https://www.mdanalysis.org/pages/conduct/
.. _trajectory formats: https://docs.mdanalysis.org/documentation_pages/coordinates/init.html#id1
.. _topology formats: https://docs.mdanalysis.org/documentation_pages/topology/init.html#supported-topology-formats
.. _Issue 87: https://github.com/MDAnalysis/mdanalysis/issues/87
.. _MDAnalysis: https://www.mdanalysis.org
.. _LICENSE: https://github.com/MDAnalysis/mdanalysis/blob/master/LICENSE
.. _`#286`: https://github.com/MDAnalysis/mdanalysis/issues/286
.. _`MDAnalysis.analysis`: https://docs.mdanalysis.org/documentation_pages/analysis_modules.html
.. _`tutorials`: https://www.mdanalysis.org/pages/learning_MDAnalysis/
.. _`guide`: https://github.com/MDAnalysis/mdanalysis/wiki/Guide-for-Developers

.. |usergroup| image:: https://img.shields.io/badge/Google%20Group-Users-lightgrey.svg
   :alt: User Google Group
   :target: http://users.mdanalysis.org

.. |developergroup| image:: https://img.shields.io/badge/Google%20Group-Developers-lightgrey.svg
   :alt: Developer Google Group
   :target: http://developers.mdanalysis.org

.. |docs| image:: https://img.shields.io/badge/docs-latest-brightgreen.svg
   :alt: Documentation (latest release)
   :target: https://docs.mdanalysis.org

.. |devdocs| image:: https://img.shields.io/badge/docs-development-yellow.svg
   :alt: Documentation (development version)
   :target: https://www.mdanalysis.org/mdanalysis/

.. |numfocus| image:: https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A
   :alt: Powered by NumFOCUS
   :target: https://www.numfocus.org/

.. |build| image:: https://travis-ci.com/MDAnalysis/mdanalysis.svg?branch=develop
   :alt: Build Status
   :target: https://travis-ci.com/MDAnalysis/mdanalysis

.. |cov|   image:: https://codecov.io/gh/MDAnalysis/mdanalysis/branch/develop/graph/badge.svg
   :alt: Coverage Status
   :target: https://codecov.io/gh/MDAnalysis/mdanalysis

.. |anaconda| image:: https://anaconda.org/conda-forge/mdanalysis/badges/version.svg
   :alt: Anaconda
   :target: https://anaconda.org/conda-forge/mdanalysis

.. |mybinder| image:: https://mybinder.org/badge.svg
   :alt: My Binder
   :target: https://mybinder.org/v2/gh/MDAnalysis/binder-notebook/master
