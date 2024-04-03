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
from collections import OrderedDict

import pytest
from numpy.testing import (
    assert_equal,assert_allclose
)

import MDAnalysis as mda
import os

from MDAnalysisTests.datafiles import CRD
from MDAnalysisTests import make_Universe


class TestCRDWriter(object):
    @pytest.fixture()
    def u(self):
        return mda.Universe(CRD)

    @pytest.fixture()
    def outfile(self, tmpdir):
        return os.path.join(str(tmpdir), 'test.crd')

    @pytest.mark.parametrize('testfile',
        ['test.crd', 'test.crd.bz2', 'test.crd.gz'])
    def test_write_atoms(self, u, testfile, tmpdir):
        # Test that written file when read gives same coordinates
        with tmpdir.as_cwd():
            u.atoms.write(testfile)

            u2 = mda.Universe(testfile)

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

    def test_write_EXT(self, u, outfile):
        # Use the `extended` keyword to force the EXT format
        u.atoms.write(outfile, extended=True)

        with open(outfile, 'r') as inf:
            format_line = inf.readlines()[2]
        assert 'EXT' in format_line, "EXT format expected"

    def test_write_EXT_read(self, u, outfile):
        # Read EXT format and check atom positions
        u.atoms.write(outfile, extended=True)

        u2 = mda.Universe(outfile)

        sel1 = u.select_atoms('all')
        sel2 = u2.select_atoms('all')

        cog1 = sel1.center_of_geometry()
        cog2 = sel2.center_of_geometry()

        assert_equal(len(u.atoms), len(u2.atoms)), 'Equal number of '\
                'atoms expected in both CRD formats'
        assert_equal(len(u.atoms.residues),
            len(u2.atoms.residues)), 'Equal number of residues expected in'\
                        'both CRD formats'
        assert_equal(len(u.atoms.segments), 
            len(u2.atoms.segments)), 'Equal number of segments expected in'\
                        'both CRD formats'
        assert_allclose(cog1, cog2, rtol=1e-6, atol=0), 'Same centroid expected for both CRD formats'


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
