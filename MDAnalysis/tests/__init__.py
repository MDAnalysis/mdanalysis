# -*- coding: utf-8 -*-

"""
=========================
Test cases for MDAnalysis
=========================

We are using the NumPy_ testing frame work; thus, numpy *must* be
installed for the tests to run at all.

Run the tests with
   >>> import MDAnalysis.tests
   >>> MDAnalysis.tests.test()

Note that if no tests are being run then one has to remove the executable bit
from all *.py files in the test directory::
  chmod a-X ..../MDAnalysis/tests/*.py
as described in `nose, subversion and executable bits`_.

.. _nose, subversion and executable bits:
   http://lacostej.blogspot.com/2009/09/nose-subversion-and-exectubale-bits.html

The simulation data used in some tests are from [Beckstein2009]_ (adk.psf,
adk_dims.dcd) or unpublished simulations (O. Beckstein).

[Beckstein2009] O. Beckstein, E.J. Denning, J.R. Perilla and T.B. Woolf,
                Zipping and Unzipping of Adenylate Kinase: Atomistic Insights
                into the Ensemble of Open ↔ Closed Transitions. J Mol Biol 394
                (2009), 160–176, doi:10.1016/j.jmb.2009.09.009
"""

try:
    from numpy.testing import Tester
    test = Tester().test
except ImportError:
    raise ImportError("""numpy>=1.3  is required to run the test suite. Please install it first. """
                      """(For example, try "easy_install 'numpy>=1.3'".""")

