# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
# Beckstein and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from six.moves import cPickle

import numpy as np
from numpy.testing import (
    dec,
    assert_,
    assert_allclose,
    assert_almost_equal,
    assert_equal,
    assert_raises,
)
from nose.plugins.attrib import attr

from MDAnalysisTests.core.groupbase import make_Universe
from MDAnalysisTests.datafiles import (
    PSF, DCD,
    PSF_BAD,
    PDB_small,
)
from MDAnalysisTests import parser_not_found

import MDAnalysis as mda
from MDAnalysis.topology.base import TopologyReader


class IOErrorParser(TopologyReader):
    def parse(self):
        raise IOError("Useful information")


class TestUniverseCreation(object):
    # tests concerning Universe creation and errors encountered
    @staticmethod
    def test_load():
        # Universe(top, trj)
        u = mda.Universe(PSF, PDB_small)
        assert_equal(len(u.atoms), 3341, "Loading universe failed somehow")

    @staticmethod
    def test_make_universe_no_args():
        # universe creation without args should work
        u = mda.Universe()

        assert_(isinstance(u, mda.Universe))
        assert_(u.atoms == None)

    @staticmethod
    def test_Universe_no_trajectory_AE():
        # querying trajectory without a trajectory loaded (only topology)
        u = make_Universe()

        assert_raises(AttributeError, getattr, u, 'trajectory')

    @staticmethod
    def test_Universe_topology_IE():
        assert_raises(IOError,
                      mda.Universe, 'thisfile', topology_format=IOErrorParser)

    @staticmethod
    def test_Universe_topology_IE_msg():
        # should get the original error, as well as Universe error
        try:
            mda.Universe('thisfile', topology_format=IOErrorParser)
        except IOError as e:
            assert_('Failed to load from the topology file' in e.args[0])
            assert_('Useful information' in e.args[0])
        else:
            raise AssertionError

    @staticmethod
    def test_load_new_VE():
        u = mda.Universe()

        assert_raises(TypeError,
                      u.load_new, 'thisfile', format='soup')

    @staticmethod
    def test_universe_kwargs():
        u = mda.Universe(PSF, PDB_small, fake_kwarg=True)
        assert_equal(len(u.atoms), 3341, "Loading universe failed somehow")

        assert_(u.kwargs['fake_kwarg'] is True)

        # initialize new universe from pieces of existing one
        u2 = mda.Universe(u.filename, u.trajectory.filename,
                          **u.kwargs)

        assert_(u2.kwargs['fake_kwarg'] is True)
        assert_equal(u.kwargs, u2.kwargs)


class TestUniverse(object):
    # older tests, still useful
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def test_load_bad_topology(self):
        # tests that Universe builds produce the right error message
        def bad_load():
            return mda.Universe(PSF_BAD, DCD)

        assert_raises(ValueError, bad_load)

    @attr('issue')
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def test_load_new(self):
        u = mda.Universe(PSF, DCD)
        u.load_new(PDB_small)
        assert_equal(len(u.trajectory), 1, "Failed to load_new(PDB)")

    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def test_load_new_TypeError(self):
        u = mda.Universe(PSF, DCD)

        def bad_load(uni):
            return uni.load_new('filename.notarealextension')

        assert_raises(TypeError, bad_load, u)

    def test_load_structure(self):
        # Universe(struct)
        ref = mda.Universe(PSF, PDB_small)
        u = mda.Universe(PDB_small)
        assert_equal(len(u.atoms), 3341, "Loading universe failed somehow")
        assert_almost_equal(u.atoms.positions, ref.atoms.positions)

    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def test_load_multiple_list(self):
        # Universe(top, [trj, trj, ...])
        ref = mda.Universe(PSF, DCD)
        u = mda.Universe(PSF, [DCD, DCD])
        assert_equal(len(u.atoms), 3341, "Loading universe failed somehow")
        assert_equal(u.trajectory.n_frames, 2 * ref.trajectory.n_frames)

    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def test_load_multiple_args(self):
        # Universe(top, trj, trj, ...)
        ref = mda.Universe(PSF, DCD)
        u = mda.Universe(PSF, DCD, DCD)
        assert_equal(len(u.atoms), 3341, "Loading universe failed somehow")
        assert_equal(u.trajectory.n_frames, 2 * ref.trajectory.n_frames)

    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def test_pickle_raises_NotImplementedError(self):
        u = mda.Universe(PSF, DCD)
        assert_raises(NotImplementedError, cPickle.dumps, u, protocol=cPickle.HIGHEST_PROTOCOL)

    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def test_set_dimensions(self):
        u = mda.Universe(PSF, DCD)
        box = np.array([10, 11, 12, 90, 90, 90])
        u.dimensions = np.array([10, 11, 12, 90, 90, 90])
        assert_allclose(u.dimensions, box)
