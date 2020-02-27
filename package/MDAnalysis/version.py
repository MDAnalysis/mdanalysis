# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""
Version information for MDAnalysis - :mod:`MDAnalysis.version`
==============================================================

The version information in :mod:`MDAnalysis.version` indicates the
release of MDAnalysis. MDAnalysis uses `semantic versioning`_ as
described in the wiki page on `versioning of MDAnalysis`_.

In brief:

Given a version number MAJOR.MINOR.PATCH, we increment the

1. **MAJOR** version when we make **incompatible API changes**,
2. **MINOR** version when we **add functionality** in a
   **backwards-compatible** manner, and
3. **PATCH** version when we make backwards-compatible **bug fixes**.

However, as long as the **MAJOR** number is **0** (i.e. the API has
not stabilized), even **MINOR** increases *may* introduce incompatible
API changes. As soon as we have a 1.0.0 release, the public API can
only be changed in a backward-incompatible manner with an increase in
MAJOR version.

Additional labels for pre-release and build metadata are available as
extensions to the MAJOR.MINOR.PATCH format.


.. Note:: Development versions and pre-releases have a suffix after
          the release number, such as ``0.11.0-dev``. If you have
          problems, try out a full release (e.g. ``0.11.0``) first.

.. _`semantic versioning`: http://semver.org/
.. _`versioning of MDAnalysis`:
   https://github.com/MDAnalysis/mdanalysis/wiki/SemanticVersioning

Data
----

.. autodata:: __version__

"""

# keep __version__ in separate file to avoid circular imports
# e.g. with lib.log

#: Release of MDAnalysis as a string, using `semantic versioning`_.
__version__ = "0.20.2-dev0"  # NOTE: keep in sync with RELEASE in setup.py
