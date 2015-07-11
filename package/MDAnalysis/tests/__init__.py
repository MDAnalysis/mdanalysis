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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#


"""
=========================
Test cases for MDAnalysis
=========================

Test cases are stored in a separate packacge (:mod:`MDAnalysisTests`)
that has to be installed separately from
http://www.mdanalysis.org. For more details on the tests see
http://wiki.mdanalysis.org/UnitTests .

We are using the NumPy_ testing frame work; thus, numpy *must* be
installed for the tests to run at all.

Run all the tests with

   >>> import MDAnalysis.tests
   >>> MDAnalysis.tests.test()

Some tests can take a few seconds; in order to skip the slow tests run

   >>> MDAnalysis.tests.test(label='fast')

Additional information is displayed at a higher verbosity level (the default is
1):

   >>> MDAnalysis.tests.test(label='fast', argv=['--verbosity=3'])

Note that if no tests are being run then one might have to run the
tests with the ``--exe`` flag

   >>> MDAnalysis.tests.test(label='fast', argv=['--exe'])

(This happens when python files are installed with the executable bit set. By
default the nose_ testing framework refuses to use those files and must be
encouraged to do so with the ``--exe`` switch.)

See `nose commandline options`_ for additional options that can be used.

For the particular case of code coverage MDAnalysis mustn't be imported prior
to testing, and testing must be invoked directly from `:meth:MDAnalysisTests.run`:

  >>> import MDAnalysisTests
  >>> MDAnalysisTests.run(argv=['--exe', '--with-coverage', '--cover-package=MDAnalysis'])


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


.. [Beckstein2009] O. Beckstein, E.J. Denning, J.R. Perilla and
   T.B. Woolf, Zipping and Unzipping of Adenylate Kinase: Atomistic
   Insights into the Ensemble of Open ↔ Closed Transitions. J Mol Biol
   394 (2009), 160–176, doi:10.1016/j.jmb.2009.09.009


Writing test cases
==================

The unittests use the :mod:`unittest` module together with nose_. See the
examples provided alongside the ``MDAnalysisTests`` module.

The `SciPy testing guidelines`_ are also a good howto for writing test cases.


.. _nose:
   http://somethingaboutorange.com/mrl/projects/nose/0.11.3/index.html
.. _nose commandline options:
   http://somethingaboutorange.com/mrl/projects/nose/0.11.3/usage.html#extended-usage
.. _SciPy testing guidelines:
   http://projects.scipy.org/numpy/wiki/TestingGuidelines#id11
.. _Charmm: http://www.charmm.org
.. _Gromacs: http://www.gromacs.org

"""

try:
    from MDAnalysisTests import run as test
except ImportError:
    print("Install MDAnalysisTests first. The source package is available from")
    print("http://pypi.python.org/pypi/MDAnalysisTests")
    raise ImportError("Package MDAnalysisTests required!")
