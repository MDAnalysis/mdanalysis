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

import pytest
from six.moves import zip

import MDAnalysis as mda
import numpy as np
from numpy.testing import (
    assert_almost_equal,
)

from MDAnalysis.coordinates.XYZ import XYZWriter

from MDAnalysisTests.datafiles import COORDINATES_XYZ, COORDINATES_XYZ_BZ2
from MDAnalysisTests.coordinates.base import (MultiframeReaderTest, BaseReference,
                                              BaseWriterTest)
from MDAnalysisTests import make_Universe


class XYZReference(BaseReference):
    def __init__(self):
        super(XYZReference, self).__init__()
        self.trajectory = COORDINATES_XYZ
        # XYZ is it's own Topology
        self.topology = COORDINATES_XYZ
        self.reader = mda.coordinates.XYZ.XYZReader
        self.writer = mda.coordinates.XYZ.XYZWriter
        self.ext = 'xyz'
        self.volume = 0
        self.dimensions = np.zeros(6)
        self.container_format = True


class TestXYZReader(MultiframeReaderTest):
    @staticmethod
    @pytest.fixture()
    def ref():
        return XYZReference()

    def test_double_open(self):
        with pytest.raises(Exception):
            self.reader.open_trajectory()
            self.reader.open_trajectory()


class TestXYZWriter(BaseWriterTest):
    @staticmethod
    @pytest.fixture()
    def ref():
        return XYZReference()

    def test_write_selection(self, ref, reader, tempdir):
        uni = mda.Universe(ref.topology, ref.trajectory)
        sel_str = 'name CA'
        sel = uni.select_atoms(sel_str)
        outfile = self.tmp_file('write-selection-test', ref, tempdir)

        with ref.writer(outfile, sel.n_atoms) as W:
            for ts in uni.trajectory:
                W.write(sel.atoms)

        copy = ref.reader(outfile)
        for orig_ts, copy_ts in zip(uni.trajectory, copy):
            assert_almost_equal(
                copy_ts._pos, sel.atoms.positions, ref.prec,
                err_msg="coordinate mismatch between original and written "
                        "trajectory at frame {} (orig) vs {} (copy)".format(
                    orig_ts.frame, copy_ts.frame))

    def test_write_different_models_in_trajectory(self, ref, reader, tempdir):
        outfile = self.tmp_file('write-models-in-trajectory', ref, tempdir)
        # n_atoms should match for each TimeStep if it was specified
        with ref.writer(outfile, n_atoms=4) as w:
            with pytest.raises(ValueError):
                w.write(reader.ts)

    def test_no_conversion(self, ref, reader, tempdir):
        outfile = self.tmp_file('write-no-conversion', ref, tempdir)
        with ref.writer(outfile, convert_units=False) as w:
            for ts in reader:
                w.write(ts)
        self._check_copy(outfile, ref, reader)


class XYZ_BZ_Reference(XYZReference):
    def __init__(self):
        super(XYZ_BZ_Reference, self).__init__()
        self.trajectory = COORDINATES_XYZ_BZ2
        self.ext = 'xyz.bz2'


class Test_XYZBZReader(TestXYZReader):
    @staticmethod
    @pytest.fixture()
    def ref():
        return XYZ_BZ_Reference()


class Test_XYZBZWriter(TestXYZWriter):
    @staticmethod
    @pytest.fixture()
    def ref():
        return XYZ_BZ_Reference()


class TestXYZWriterNames(object):
    @pytest.fixture()
    def outfile(self, tmpdir):
        return str(tmpdir.join('/outfile.xyz'))

    def test_no_names(self, outfile):
        u = make_Universe(trajectory=True)

        w = XYZWriter(outfile)
        w.write(u.trajectory.ts)
        w.close()

        u2 = mda.Universe(outfile)
        assert all(u2.atoms.names == 'X')

    def test_single_name(self, outfile):
        u = make_Universe(trajectory=True)

        w = XYZWriter(outfile, atoms='ABC')
        w.write(u.trajectory.ts)
        w.close()

        u2 = mda.Universe(outfile)
        assert all(u2.atoms.names == 'ABC')

    def test_list_names(self, outfile):
        u = make_Universe(trajectory=True)

        names = ['A', 'B', 'C', 'D', 'E'] * 25

        w = XYZWriter(outfile, atoms=names)
        w.write(u.trajectory.ts)
        w.close()

        u2 = mda.Universe(outfile)
        assert all(u2.atoms.names == names)
