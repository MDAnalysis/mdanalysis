# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
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

"""=========================
Test cases for MDAnalysis
=========================

The test cases and the test data are kept in this package,
MDAnalysisTests. They will only run when MDAnalysis is also
installed. MDAnalysis and MDAnalysisTests *must* have the same release
number, which can be found in :data:`MDAnalysis.__version__` and
:data:`MDAnalysisTests.__version__`. If the versions don't match then
an :exc:`ImportError` is raised.

We are using the NumPy_ testing frame work; thus, :mod:`numpy` *must* be
installed for the tests to run at all.

Run all the tests with ::

    pytest --pyargs MDAnalysisTests

If you have the `pytest-xdist`_ plugin installed then you can run the
tests in parallel

.. code-block:: bash

   pytest -n 4 --pyargs MDAnalysisTests


.. _`pytest-xdist`:
   https://github.com/pytest-dev/pytest-xdist


Data
====

The simulation data used in some tests are from [Beckstein2009]_ (``adk.psf``,
``adk_dims.dcd``) or unpublished simulations (O. Beckstein).

   adk_dims
      Trajectory of a macromolecular transition of the enzyme adenylate kinase
      between a closed and an open conformation. The simulation was run in
      Charmm_ c35a1.

   adk_oplsaa
      Ten frames from the first 1 ns of a equilibrium trajectory of AdK in
      water with Na+ counter ions. The OPLS/AA forcefield is used with the
      TIP4P water model. The simulation was run with Gromacs_ 4.0.2.


[Beckstein2009] O. Beckstein, E.J. Denning, J.R. Perilla and T.B. Woolf,
                Zipping and Unzipping of Adenylate Kinase: Atomistic Insights
                into the Ensemble of Open ↔ Closed Transitions. J Mol Biol 394
                (2009), 160–176, doi:10.1016/j.jmb.2009.09.009


Writing test cases
==================

The unittests use the :mod:`pytest <https://docs.pytest.org/en/latest/>`_ module. See the
examples in the ``MDAnalysisTests`` directory.

The `SciPy testing guidelines`_ are a good howto for writing test cases,
especially as we are directly using this framework (imported from numpy).

.. _NumPy: http://www.numpy.org/
.. _SciPy testing guidelines:
   http://projects.scipy.org/numpy/wiki/TestingGuidelines#id11
.. _Charmm: http://www.charmm.org
.. _Gromacs: http://www.gromacs.org

"""
import logging

import pytest

logger = logging.getLogger("MDAnalysisTests.__init__")

# keep in sync with RELEASE in setup.py
__version__ = "2.8.0-dev0"


# Do NOT import MDAnalysis at this level. Tests should do it themselves.
# If MDAnalysis is imported here coverage accounting might fail because all the import
#  code won't be run again under coverage's watch. See Issue 344.

import os
import sys

# To test duecredit in utils/test_duecredits.py.
#
# Note that the test environment on travis should have duecredit installed
# (put it in PIP_DEPENDENCIES). Setting DUECREDIT_ENABLE to yes triggers
# collection of citations on the first `import MDAnalysis` so the environment
# variable *must* come before MDAnalysis is imported the first time. See
# issue #412 https://github.com/MDAnalysis/mdanalysis/issues/412 and PR #1822.
os.environ['DUECREDIT_ENABLE'] = 'yes'

# Any tests that plot with matplotlib need to run with the simple agg backend
# because on Travis there is no DISPLAY set.
#
# Instead of using matplotlib.use() we set MPLBACKEND=agg in the CI environment.
# See https://matplotlib.org/3.2.1/tutorials/introductory/usage.html#backends

from MDAnalysisTests.util import (
    block_import,
    executable_not_found,
    no_deprecated_call,
    in_dir,
    assert_nowarns,
)
from MDAnalysisTests.dummy import make_Universe


def run(*args, **kwargs):
    pytest.main()
