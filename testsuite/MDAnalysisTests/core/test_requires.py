"""Tests for core.groups.requires decorator

"""
import numpy as np
from numpy.testing import (
    assert_,
    assert_raises,
)

from MDAnalysis.core.groups import requires
from MDAnalysis import NoDataError

from MDAnalysisTests.core.groupbase import make_Universe


class TestRequires(object):
    def test_requires_failure_singular(self):
        @requires('masses')
        def mass_multiplier(ag1, ag2, scalar):
            return (ag1.masses + ag2.masses) * scalar

        u = make_Universe('charges')

        assert_raises(NoDataError, mass_multiplier, u.atoms[:10], u.atoms[20:30], 4.0)

    def test_requires_failure_multiple(self):
        @requires('masses', 'charges')
        def mass_multiplier(ag1, ag2, scalar):
            return (ag1.masses + ag2.charges) * scalar
        

        u = make_Universe('masses', 'types')

        assert_raises(NoDataError, mass_multiplier, u.atoms[:10], u.atoms[20:30], 4.0)

    def test_requires_success(self):
        @requires('masses')
        def mass_multiplier(ag1, ag2, scalar):
            return (ag1.masses + ag2.masses) * scalar

        u = make_Universe('masses')

        result = mass_multiplier(u.atoms[:10], u.atoms[20:30], 4.0)

        assert_(isinstance(result, np.ndarray))
