================================
  MDAnalysis Repository README
================================

|numfocus| |build| |travis| |cov| [*]_

|docs| |devdocs| |usergroup| |developergroup| |anaconda| |mybinder|

MDAnalysis_ is a Python library for the analysis of computer simulations of many-body systems at the molecular scale, spanning use cases from interactions of drugs with proteins to novel materials. It is widely used in the scientific community and is written by scientists for scientists. 

It works with a wide range of popular simulation packages including Gromacs, Amber, NAMD, CHARMM, DL_Poly, HooMD, LAMMPS and many others — see the lists of supported `trajectory formats`_ and `topology formats`_.
MDAnalysis also includes widely used analysis algorithms in the `MDAnalysis.analysis`_ module.

.. _numfocus-fiscal-sponsor-attribution:

The MDAnalysis project uses an `open governance model`_ and is fiscally sponsored by `NumFOCUS`_. Consider making 
a `tax-deductible donation`_ to help the project pay for developer time, professional services, travel, workshops, and a variety of other needs.

.. image:: https://www.mdanalysis.org/public/images/numfocus-sponsored-small.png
   :alt: NumFOCUS (Fiscally Sponsored Project)
   :target: https://numfocus.org/project/mdanalysis
   :align: center

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


Documentation
=============

**New users** should read the `Quickstart Guide`_ and might want to
look at our videos_, in which core developers explain various aspects
of MDAnalysis.

**All users** should read the `User Guide`_.

**Developers** may also want to refer to the `MDAnalysis API docs`_.

A growing number of `tutorials`_ are available that explain how to
conduct RMSD calculations, structural alignment, distance and contact
analysis, and many more.


Installation and availability
=============================

The latest release can be **installed via ``pip`` or ``conda``** as
described in the `Installation Quick Start`_.

**Source code** is hosted in a git repository at
https://github.com/MDAnalysis/mdanalysis and is available under the
GNU General Public License, version 2 (see the file LICENSE_).


Contributing
============

Please report **bugs** or **enhancement requests** through the `Issue
Tracker`_. Questions can also be asked on the `user mailing list`_.

If you are a **new developer** who would like to start contributing to
MDAnalysis get in touch on the `developer mailing list`_. To set up a
development environment and run the test suite read the `developer
guide`_.


Citation
========

When using MDAnalysis in published work, please cite the following
two papers:

*   R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy,
    M. N. Melo, S. L. Seyler, D. L. Dotson, J. Domanski,
    S. Buchoux, I. M. Kenney, and O. Beckstein. MDAnalysis:
    A Python package for the rapid analysis of molecular
    dynamics simulations. In S. Benthall and S. Rostrup,
    editors, Proceedings of the 15th Python in Science
    Conference, pages 102-109, Austin, TX, 2016. SciPy.
    doi:`10.25080/Majora-629e541a-00e`_    

*   N. Michaud-Agrawal, E. J. Denning, T. B. Woolf,
    and O. Beckstein. MDAnalysis: A Toolkit for the Analysis of Molecular
    Dynamics Simulations. *J. Comput. Chem.* **32** (2011), 2319--2327.
    doi:`10.1002/jcc.21787`_

For citations of included algorithms and sub-modules please see the references_.



.. Footnotes

.. [*] **build**: Unit testing is for the whole package; **coverage** is
       shown for the core library modules and the analysis modules.

.. _NumFOCUS: https://numfocus.org/
.. _open governance model: https://www.mdanalysis.org/about/#governance
.. _tax-deductible donation: https://numfocus.org/donate-to-mdanalysis
.. _`Code of Conduct`: https://www.mdanalysis.org/pages/conduct/
.. _trajectory formats: https://docs.mdanalysis.org/documentation_pages/coordinates/init.html#id1
.. _topology formats: https://docs.mdanalysis.org/documentation_pages/topology/init.html#supported-topology-formats
.. _MDAnalysis: https://www.mdanalysis.org
.. _LICENSE:
   https://github.com/MDAnalysis/mdanalysis/blob/master/LICENSE
.. _`Installation Quick Start`:
   https://www.mdanalysis.org/pages/installation_quick_start/
.. _`MDAnalysis.analysis`: https://docs.mdanalysis.org/documentation_pages/analysis_modules.html
.. _`tutorials`: https://userguide.mdanalysis.org/examples/README.html
.. _`videos`: https://www.mdanalysis.org/pages/learning_MDAnalysis/#videos
.. _`Quickstart Guide`:
   https://userguide.mdanalysis.org/examples/quickstart.html
.. _`User Guide`: https://userguide.mdanalysis.org
.. _`MDAnalysis API docs`:
   https://docs.mdanalysis.org
.. _`Issue Tracker`: https://github.com/mdanalysis/mdanalysis/issues
.. _`user mailing list`:
   https://groups.google.com/group/mdnalysis-discussion
.. _`developer guide`:
   https://userguide.mdanalysis.org/contributing.html
.. _`developer mailing list`:
   https://groups.google.com/group/mdnalysis-devel
.. _`10.1002/jcc.21787`: https://dx.doi.org/10.1002/jcc.21787
.. _`10.25080/Majora-629e541a-00e`: https://doi.org/10.25080/Majora-629e541a-00e
.. _references: https://docs.mdanalysis.org/documentation_pages/references.html


.. |usergroup| image:: https://img.shields.io/badge/Google%20Group-Users-lightgrey.svg
   :alt: User Google Group
   :target: https://groups.google.com/group/mdnalysis-discussion

.. |developergroup| image:: https://img.shields.io/badge/Google%20Group-Developers-lightgrey.svg
   :alt: Developer Google Group
   :target: https://groups.google.com/group/mdnalysis-devel

.. |docs| image:: https://img.shields.io/badge/docs-latest-brightgreen.svg
   :alt: Documentation (latest release)
   :target: https://docs.mdanalysis.org

.. |devdocs| image:: https://img.shields.io/badge/docs-development-yellow.svg
   :alt: Documentation (development version)
   :target: https://docs.mdanalysis.org/dev

.. |numfocus| image:: https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A
   :alt: Powered by NumFOCUS
   :target: https://www.numfocus.org/

.. |build| image:: https://github.com/MDAnalysis/mdanalysis/actions/workflows/gh-ci.yaml/badge.svg
   :alt: Github Actions Build Status
   :target: https://github.com/MDAnalysis/mdanalysis/actions/workflows/gh-ci.yaml

.. |travis| image:: https://img.shields.io/travis/MDAnalysis/mdanalysis/develop?label=Travis%20CI
   :alt: Travis CI Build Status
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
