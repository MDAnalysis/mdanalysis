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
from __future__ import absolute_import

from six.moves import zip
from collections import OrderedDict

import pytest
from numpy.testing import (
    assert_equal,
)

import MDAnalysis as mda

from MDAnalysisTests.datafiles import CRD
from MDAnalysisTests import make_Universe


class TestCRDWriter(object):
    @pytest.fixture()
    def u(self):
        return mda.Universe(CRD)

    @pytest.fixture()
    def outfile(self, tmpdir):
        return str(tmpdir) + '/out.crd'

    def test_write_atoms(self, u, outfile):
        # Test that written file when read gives same coordinates
        u.atoms.write(outfile)

        u2 = mda.Universe(outfile)

        assert_equal(u.atoms.positions,
                     u2.atoms.positions)

    def test_roundtrip(self, u, outfile):
        # Write out a copy of the Universe, and compare this against the original
        # This is more rigorous than simply checking the coordinates as it checks
        # all formatting
        u.atoms.write(outfile)

        def CRD_iter(fn):
            with open(fn, 'r') as inf:
                for line in inf:
                    if not line.startswith('*'):
                        yield line

        for ref, other in zip(CRD_iter(CRD), CRD_iter(outfile)):
            assert ref == other

    def test_write_EXT(self):
        # TODO: Write tests that use EXT output format
        # Must have *lots* of atoms, maybe fake the system
        # to make tests faster
        pass


class TestCRDWriterMissingAttrs(object):
    # All required attributes with the default value
    req_attrs = OrderedDict([
        ('resnames', 'UNK'),
        ('resids', 1),
        ('names', 'X'),
        ('tempfactors', 0.0),
    ])

    @pytest.mark.parametrize('missing_attr', req_attrs)
    def test_warns(self, missing_attr, tmpdir):
        attrs = list(self.req_attrs.keys())
        attrs.remove(missing_attr)
        u = make_Universe(attrs, trajectory=True)

        outfile = str(tmpdir) + '/out.crd'
        with pytest.warns(UserWarning):
            u.atoms.write(outfile)

    @pytest.mark.parametrize('missing_attr', req_attrs)
    def test_write(self, missing_attr, tmpdir):
        attrs = list(self.req_attrs.keys())
        attrs.remove(missing_attr)
        u = make_Universe(attrs, trajectory=True)

        outfile = str(tmpdir) + '/out.crd'
        u.atoms.write(outfile)
        u2 = mda.Universe(outfile)

        # Check all other attrs aren't disturbed
        for attr in attrs:
            assert_equal(getattr(u.atoms, attr),
                         getattr(u2.atoms, attr))
        # Check missing attr is as expected
        assert_equal(getattr(u2.atoms, missing_attr),
                     self.req_attrs[missing_attr])
