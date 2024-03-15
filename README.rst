================================
  MDAnalysis Repository README
================================

|numfocus| 

|build| |cron| |cirruscron| |linters| |cov|

|docs| |devdocs| |discussions|

|anaconda| |asv|


This is my forked version of MDAnalysis. I will modify some code for special use.  

--------------------------------


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

|powered_by_MDA|

If you use MDAnalysis_ in your project consider letting your users and the world know about it by displaying the MDAnalysis_ badge! `Embedding code`_ is available for different markups.

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

The latest release can be **installed via pip or conda** as
described in the `Installation Quick Start`_.

**Source code** is hosted in a git repository at
https://github.com/MDAnalysis/mdanalysis and is packaged under the
GNU General Public License, version 3 or any later version. Invidiual
source code components are provided under a mixture of GPLv3+ compatible
licenses, including LGPLv2.1+ and GPLv2+. Please see the file LICENSE_
for more information.


Contributing
============

Please report **bugs** or **enhancement requests** through the `Issue
Tracker`_. Questions can also be asked on `GitHub Discussions`_.

If you are a **new developer** who would like to start contributing to
MDAnalysis get in touch on `GitHub Discussions`_. To set up a
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


.. _NumFOCUS: https://numfocus.org/
.. _open governance model: https://www.mdanalysis.org/about/#governance
.. _tax-deductible donation: https://numfocus.org/donate-to-mdanalysis
.. _`Code of Conduct`: https://www.mdanalysis.org/pages/conduct/
.. _trajectory formats: https://docs.mdanalysis.org/documentation_pages/coordinates/init.html#id1
.. _topology formats: https://docs.mdanalysis.org/documentation_pages/topology/init.html#supported-topology-formats
.. _MDAnalysis: https://www.mdanalysis.org
.. _LICENSE:
   https://github.com/MDAnalysis/mdanalysis/blob/develop/LICENSE
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
.. _`GitHub Discussions`:
   https://github.com/MDAnalysis/mdanalysis/discussions
.. _`developer guide`:
   https://userguide.mdanalysis.org/contributing.html
.. _`10.1002/jcc.21787`: https://dx.doi.org/10.1002/jcc.21787
.. _`10.25080/Majora-629e541a-00e`: https://doi.org/10.25080/Majora-629e541a-00e
.. _references: https://docs.mdanalysis.org/documentation_pages/references.html
.. _Embedding code: https://www.mdanalysis.org/pages/citations/#powered-by-mdanalysis

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

.. |cron| image:: https://github.com/MDAnalysis/mdanalysis/actions/workflows/gh-ci-cron.yaml/badge.svg
   :alt: Github Actions Cron Job Status
   :target: https://github.com/MDAnalysis/mdanalysis/actions/workflows/gh-ci-cron.yaml

.. |cirruscron| image:: https://img.shields.io/cirrus/github/MDAnalysis/mdanalysis/develop
   :alt: Cirrus CI - Cron job status
   :target: https://cirrus-ci.com/github/MDAnalysis/mdanalysis/develop

.. |linters| image:: https://github.com/MDAnalysis/mdanalysis/actions/workflows/linters.yaml/badge.svg
   :alt: Github Actions Linters Status
   :target: https://github.com/MDAnalysis/mdanalysis/actions/workflows/linters.yaml

.. |cov|   image:: https://codecov.io/gh/MDAnalysis/mdanalysis/branch/develop/graph/badge.svg
   :alt: Coverage Status
   :target: https://codecov.io/gh/MDAnalysis/mdanalysis

.. |anaconda| image:: https://anaconda.org/conda-forge/mdanalysis/badges/version.svg
   :alt: Anaconda
   :target: https://anaconda.org/conda-forge/mdanalysis
   
.. |asv| image:: https://img.shields.io/badge/benchmarked%20by-asv-blue.svg
   :alt: ASV Benchmarks
   :target:  https://www.mdanalysis.org/benchmarks/

.. |powered_by_MDA| image:: https://img.shields.io/badge/Powered%20by-MDAnalysis-orange.svg?logoWidth=15&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA
   :alt: Powered by MDAnalysis
   :target: https://www.mdanalysis.org

.. |discussions| image:: https://img.shields.io/github/discussions/MDAnalysis/MDAnalysis
   :alt: GitHub Discussions
   :target: https://github.com/MDAnalysis/mdanalysis/discussions
