"""Tests for core.groups.requires decorator

"""
from __future__ import absolute_import
import numpy as np
from numpy.testing import (
    assert_,
    assert_raises,
)

from MDAnalysis.core.groups import requires
from MDAnalysis import NoDataError

from MDAnalysisTests import make_Universe


class TestRequires(object):
    def test_requires_failure_singular(self):
        @requires('masses')
        def mass_multiplier(ag1, ag2, scalar):
            return (ag1.masses + ag2.masses) * scalar

        u = make_Universe(('charges',))

        assert_raises(NoDataError, mass_multiplier, u.atoms[:10], u.atoms[20:30], 4.0)

    def test_requires_failure_multiple(self):
        @requires('masses', 'charges')
        def mass_multiplier(ag1, ag2, scalar):
            return (ag1.masses + ag2.charges) * scalar
        

        u = make_Universe(('masses', 'types'))

        assert_raises(NoDataError, mass_multiplier, u.atoms[:10], u.atoms[20:30], 4.0)

    def test_requires_success(self):
        @requires('masses')
        def mass_multiplier(ag1, ag2, scalar):
            return (ag1.masses + ag2.masses) * scalar

        u = make_Universe(('masses',))

        result = mass_multiplier(u.atoms[:10], u.atoms[20:30], 4.0)

        assert_(isinstance(result, np.ndarray))

    def test_failure_errormessage(self):
        # failures should list all required attributes not
        # just the first one
        @requires('cats', 'dogs', 'frogs')
        def animal_print(ag):
            return len(ag.cats), len(ag.dogs), len(ag.frogs)

        u = make_Universe()
        try:
            animal_print(u.atoms)
        except NoDataError as e:
            message = e.args[0]
            # Test function name gets returned (for debug)
            assert_('animal_print' in message)
            assert_('cats' in message)
            assert_('dogs' in message)
            assert_('frogs' in message)
        else:
            raise AssertionError("Should raise NoDataError")
