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


The simulation data used in some tests are from [Beckstein2009]_.

[Beckstein2009] O. Beckstein, E.J. Denning, J.R. Perilla and T.B. Woolf,
                Zipping and Unzipping of Adenylate Kinase: Atomistic Insights
                into the Ensemble of Open ↔ Closed Transitions. J Mol Biol 394
                (2009), 160–176, doi:10.1016/j.jmb.2009.09.009
"""

try:
    from numpy.testing import Tester
    test = Tester().test
except ImportError:
    raise ImportError("NumPy is required to run the test suite. Please install it first.")

