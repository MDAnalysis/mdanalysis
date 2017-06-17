# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.mdanalysis.org
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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import absolute_import
from six.moves import range, zip

import MDAnalysis as mda
import numpy as np
import os
from numpy.testing import (
    assert_array_almost_equal,
    raises,
    assert_array_equal,
    assert_,
)

from MDAnalysis.coordinates.XYZ import XYZWriter

from MDAnalysisTests.datafiles import COORDINATES_XYZ, COORDINATES_XYZ_BZ2
from MDAnalysisTests.coordinates.base import (MultiframeReaderTest, BaseReference,
                                              BaseWriterTest)
from MDAnalysisTests import tempdir, make_Universe


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
    def __init__(self, reference=None):
        if reference is None:
            reference = XYZReference()
        super(TestXYZReader, self).__init__(reference)

    @raises
    def test_double_open(self):
        self.reader.open_trajectory()
        self.reader.open_trajectory()


class TestXYZWriter(BaseWriterTest):
    def __init__(self, reference=None):
        if reference is None:
            reference = XYZReference()
        super(TestXYZWriter, self).__init__(reference)

    def test_write_selection(self):
        uni = mda.Universe(self.ref.topology, self.ref.trajectory)
        sel_str = 'name CA'
        sel = uni.select_atoms(sel_str)
        outfile = self.tmp_file('write-selection-test')

        with self.ref.writer(outfile, sel.n_atoms) as W:
            for ts in uni.trajectory:
                W.write(sel.atoms)

        copy = self.ref.reader(outfile)
        for orig_ts, copy_ts in zip(uni.trajectory, copy):
            assert_array_almost_equal(
                copy_ts._pos, sel.atoms.positions, self.ref.prec,
                err_msg="coordinate mismatch between original and written "
                "trajectory at frame {} (orig) vs {} (copy)".format(
                    orig_ts.frame, copy_ts.frame))


    @raises(ValueError)
    def test_write_different_models_in_trajectory(self):
        outfile = self.tmp_file('write-models-in-trajectory')
        # n_atoms should match for each TimeStep if it was specified
        with self.ref.writer(outfile, n_atoms=4) as w:
            w.write(self.reader.ts)

    def test_no_conversion(self):
        outfile = self.tmp_file('write-no-conversion')
        with self.ref.writer(outfile, convert_units=False) as w:
            for ts in self.reader:
                w.write(ts)
        self._check_copy(outfile)


class XYZ_BZ_Reference(XYZReference):
    def __init__(self):
        super(XYZ_BZ_Reference, self).__init__()
        self.trajectory = COORDINATES_XYZ_BZ2
        self.ext = 'xyz.bz2'


class Test_XYZBZReader(TestXYZReader):
    def __init__(self):
        super(Test_XYZBZReader, self).__init__(XYZ_BZ_Reference())


class Test_XYZBZWriter(TestXYZWriter):
    def __init__(self):
        super(Test_XYZBZWriter, self).__init__(XYZ_BZ_Reference())


class TestXYZWriterNames(object):
    def setUp(self):
        self.tmpdir = tempdir.TempDir()
        self.outfile = self.tmpdir.name + '/outfile.xyz'

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except OSError:
            pass
        del self.outfile
        del self.tmpdir

    def test_no_names(self):
        u = make_Universe(trajectory=True)

        w = XYZWriter(self.outfile)
        w.write(u.trajectory.ts)
        w.close()

        u2 = mda.Universe(self.outfile)
        assert_(all(u2.atoms.names == 'X'))

    def test_single_name(self):
        u = make_Universe(trajectory=True)

        w = XYZWriter(self.outfile, atoms='ABC')
        w.write(u.trajectory.ts)
        w.close()

        u2 = mda.Universe(self.outfile)
        assert_(all(u2.atoms.names == 'ABC'))

    def test_list_names(self):
        u = make_Universe(trajectory=True)

        names = ['A', 'B', 'C', 'D', 'E'] * 25

        w = XYZWriter(self.outfile, atoms=names)
        w.write(u.trajectory.ts)
        w.close()

        u2 = mda.Universe(self.outfile)
        assert_(all(u2.atoms.names == names))
