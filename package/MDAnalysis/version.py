# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""
Version information for MDAnalysis — :mod:`MDAnalysis.version`
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
   http://wiki.mdanalysis.org/SemanticVersioning

Data
----

.. autodata:: __version__

"""

# keep __version__ in separate file to avoid circular imports
# e.g. with lib.log

#: Release of MDAnalysis as a string, using `semantic versioning`_.
__version__ = "0.11.1-dev"  # NOTE: keep in sync with RELEASE in setup.py
