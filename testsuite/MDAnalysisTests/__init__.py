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

The test cases and the test data are kept in this package,
MDAnalysisTests. They will only run when MDAnalysis is also
installed. MDAnalysis and MDAnalysisTests *must* have the same release
number, which can be found in :data:`MDAnalysis.__version__` and
:data:`MDAnalysisTests.__version__`. If the versions don't match then
an :exc:`ImportError` is raised.

We are using the NumPy_ testing frame work; thus, :mod:`numpy` *must* be
installed for the tests to run at all.

Run all the tests with

   >>> import MDAnalysisTests.tests
   >>> MDAnalysisTests.tests.test(label='full')

Some tests can take a few seconds; in order to skip the slow tests run

   >>> MDAnalysisTests.tests.test(label='fast')

Additional information is displayed at a higher verbosity level (the default is
1):

   >>> MDAnalysisTests.tests.test(label='fast', verbose=3)

Note that if no tests are being run then one might have to run the
tests with the ``--exe`` flag

   >>> MDAnalysisTests.tests.test(label='fast', extra_argv=['--exe'])

(This happens when python files are installed with the executable bit set. By
default the nose_ testing framework refuses to use those files and must be
encouraged to do so with the ``--exe`` switch.)

See `nose commandline options`_ for additional options that can be used; for
instance, code coverage can also be checked:

  >>> MDAnalysisTests.tests.test(label='full', extra_argv=['--exe', '--with-coverage'])


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

The unittests use the :mod:`unittest` module together with nose_. See the
examples in the ``MDAnalysisTests`` directory.

The `SciPy testing guidelines`_ are a good howto for writing test cases,
especially as we are directly using this framework (imported from numpy).

A number of plugins external to nose are automatically loaded. the `knownfailure`
plugin provides the `@knownfailure()` decorator, which can be used to mark tests
that are expected to fail. Beware that even if used with default arguments the
parentheses must be included.

.. _NumPy: http://www.numpy.org/
.. _nose:
   http://somethingaboutorange.com/mrl/projects/nose/0.11.3/index.html
.. _nose commandline options:
   http://somethingaboutorange.com/mrl/projects/nose/0.11.3/usage.html#extended-usage
.. _SciPy testing guidelines:
   http://projects.scipy.org/numpy/wiki/TestingGuidelines#id11
.. _Charmm: http://www.charmm.org
.. _Gromacs: http://www.gromacs.org

"""

__version__ = "0.11.0-dev"  # keep in sync with RELEASE in setup.py

try:
    from numpy.testing import assert_
except ImportError:
    raise ImportError("""numpy>=1.5  is required to run the test suite. Please install it first. """
                      """(For example, try "easy_install 'numpy>=1.5'").""")

try:
    import nose
except ImportError:
    raise ImportError('nose is required to run the test suite. Please install it first. '
                      '(For example, try "pip install nose").')

import nose.plugins.multiprocess
_multiprocess_ok = hasattr(nose.plugins.multiprocess, "_instantiate_plugins")
if not _multiprocess_ok:
    raise ImportWarning("nose >= 1.1.0 is needed for multiprocess testing with external plugins, "
                        "and your setup doesn't meet this requirement. If you're running "
                        "tests in parallel external plugins will be disabled.")

try:
    import MDAnalysis
except ImportError:
    raise ImportError('Install MDAnalysis first to run the tests, '
                      'e.g. "pip install mdanalysis"')

if MDAnalysis.__version__ != __version__:
    raise ImportError("MDAnalysis release {0} must be installed to run the tests, not {1}".format(
            __version__, MDAnalysis.__version__))

from os.path import dirname
from MDAnalysisTests.plugins import loaded_plugins
import sys

def run(*args, **kwargs):
    """Test-running function that loads plugins, sets up arguments, and calls `nose.main()`"""
    try:
        kwargs['argv'] = sys.argv + kwargs['argv'] #sys.argv takes precedence
    except KeyError:
        kwargs['argv'] = sys.argv
    # We emulate numpy's treament of the 'fast' label.
    label = kwargs.get("label", None)
    if label == "fast":
        kwargs['argv'].extend(['-A','not slow'])
    # We keep accepting numpy's 'extra_argv'
    if 'extra_argv' in kwargs:
        kwargs['argv'].extend(kwargs.pop('extra_argv'))
    try:
        kwargs['addplugins'].extend(loaded_plugins.values())
    except KeyError:
        kwargs['addplugins'] = loaded_plugins.values()
    # By default, test our testsuite
    kwargs['defaultTest'] = dirname(__file__)
    return nose.main(*args, **kwargs)
# to keep compatibility with MDAnalysis.tests.test()
test = run # if we define the function directly as 'test' nose picks it up as a test to be run and recurses.

import MDAnalysis.core.util
def executable_not_found(*args):
    """Return ``True`` if not at least one of the executables in args can be found.

    ``False`` otherwise (i.e. at least one was found).
    """
    for name in args:
        found = MDAnalysis.lib.util.which(name) is not None
        if found:
            break
    return not found


def executable_not_found_runtime(*args):
    """Factory function that returns a :func:`executable_not_found`.

    The returned function has its input set to *args* but is only
    evaluated at run time.

    To be used as the argument of::

      @dec.skipif(executable_not_found_runtime("binary_name"), msg="skip test because binary_name not available")
      ...
    """
    return lambda: executable_not_found(*args)

