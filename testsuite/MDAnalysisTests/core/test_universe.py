# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.mdanalysis.org
# Copyright (c) 2006-2016 The MDAnalysis Development Team and contributors
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

from MDAnalysisTests import make_Universe
from MDAnalysisTests.datafiles import (
    PSF, DCD,
    PSF_BAD,
    PDB_small,
    GRO, TRR,
    two_water_gro, two_water_gro_nonames,
    TRZ, TRZ_psf,
)
from MDAnalysisTests import parser_not_found

import MDAnalysis as mda
import MDAnalysis.coordinates
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


class TestGuessBonds(object):
    """Test the AtomGroup methed guess_bonds

    This needs to be done both from Universe creation (via kwarg) and AtomGroup

    It needs to:
     - work if all atoms are in vdwradii table
     - fail properly if not
     - work again if vdwradii are passed.
    """
    def setUp(self):
        self.vdw = {'A':1.05, 'B':0.4}

    def tearDown(self):
        del self.vdw

    def _check_universe(self, u):
        """Verify that the Universe is created correctly"""
        assert_equal(len(u.bonds), 4)
        assert_equal(len(u.angles), 2)
        assert_equal(len(u.dihedrals), 0)
        assert_equal(len(u.atoms[0].bonds), 2)
        assert_equal(len(u.atoms[1].bonds), 1)
        assert_equal(len(u.atoms[2].bonds), 1)
        assert_equal(len(u.atoms[3].bonds), 2)
        assert_equal(len(u.atoms[4].bonds), 1)
        assert_equal(len(u.atoms[5].bonds), 1)
        assert_('guess_bonds' in u.kwargs)

    def test_universe_guess_bonds(self):
        """Test that making a Universe with guess_bonds works"""
        u = mda.Universe(two_water_gro, guess_bonds=True)
        self._check_universe(u)
        assert_(u.kwargs['guess_bonds'] is True)

    def test_universe_guess_bonds_no_vdwradii(self):
        """Make a Universe that has atoms with unknown vdwradii."""
        assert_raises(ValueError, mda.Universe, two_water_gro_nonames, guess_bonds=True)

    def test_universe_guess_bonds_with_vdwradii(self):
        """Unknown atom types, but with vdw radii here to save the day"""
        u = mda.Universe(two_water_gro_nonames, guess_bonds=True,
                                vdwradii=self.vdw)
        self._check_universe(u)
        assert_(u.kwargs['guess_bonds'] is True)
        assert_equal(self.vdw, u.kwargs['vdwradii'])

    def test_universe_guess_bonds_off(self):
        u = mda.Universe(two_water_gro_nonames, guess_bonds=False)

        for attr in ('bonds', 'angles', 'dihedrals'):
            assert_(not hasattr(u, attr))
        assert_(u.kwargs['guess_bonds'] is False)

    def _check_atomgroup(self, ag, u):
        """Verify that the AtomGroup made bonds correctly,
        and that the Universe got all this info
        """
        assert_equal(len(ag.bonds), 2)
        assert_equal(len(ag.angles), 1)
        assert_equal(len(ag.dihedrals), 0)
        assert_equal(len(u.bonds), 2)
        assert_equal(len(u.angles), 1)
        assert_equal(len(u.dihedrals), 0)
        assert_equal(len(u.atoms[0].bonds), 2)
        assert_equal(len(u.atoms[1].bonds), 1)
        assert_equal(len(u.atoms[2].bonds), 1)
        assert_equal(len(u.atoms[3].bonds), 0)
        assert_equal(len(u.atoms[4].bonds), 0)
        assert_equal(len(u.atoms[5].bonds), 0)

    def test_atomgroup_guess_bonds(self):
        """Test an atomgroup doing guess bonds"""
        u = mda.Universe(two_water_gro)

        ag = u.atoms[:3]
        ag.guess_bonds()
        self._check_atomgroup(ag, u)

    def test_atomgroup_guess_bonds_no_vdwradii(self):
        u = mda.Universe(two_water_gro_nonames)

        ag = u.atoms[:3]
        assert_raises(ValueError, ag.guess_bonds)

    def test_atomgroup_guess_bonds_with_vdwradii(self):
        u = mda.Universe(two_water_gro_nonames)

        ag = u.atoms[:3]
        ag.guess_bonds(vdwradii=self.vdw)
        self._check_atomgroup(ag, u)


class TestInMemoryUniverse(object):
    @staticmethod
    @dec.skipif(parser_not_found('DCD'),
               'DCD parser not available. Are you using python 3?')
    def test_reader_w_timeseries():
        universe = mda.Universe(PSF, DCD, in_memory=True)
        assert_equal(universe.trajectory.timeseries(universe.atoms).shape,
                     (3341, 98, 3),
                     err_msg="Unexpected shape of trajectory timeseries")

    @staticmethod
    def test_reader_wo_timeseries():
        universe = mda.Universe(GRO, TRR, in_memory=True)
        assert_equal(universe.trajectory.timeseries(universe.atoms).shape,
                     (47681, 10, 3),
                     err_msg="Unexpected shape of trajectory timeseries")

    @staticmethod
    @dec.skipif(parser_not_found('DCD'),
               'DCD parser not available. Are you using python 3?')
    def test_reader_w_timeseries_frame_interval():
        universe = mda.Universe(PSF, DCD, in_memory=True,
                                       in_memory_frame_interval=10)
        assert_equal(universe.trajectory.timeseries(universe.atoms).shape,
                     (3341, 10, 3),
                     err_msg="Unexpected shape of trajectory timeseries")

    @staticmethod
    def test_reader_wo_timeseries_frame_interval():
        universe = mda.Universe(GRO, TRR, in_memory=True,
                                       in_memory_frame_interval=3)
        assert_equal(universe.trajectory.timeseries(universe.atoms).shape,
                     (47681, 4, 3),
                     err_msg="Unexpected shape of trajectory timeseries")

    @staticmethod
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def test_existing_universe():
        universe = mda.Universe(PDB_small, DCD)
        universe.transfer_to_memory()
        assert_equal(universe.trajectory.timeseries(universe.atoms).shape,
                     (3341, 98, 3),
                     err_msg="Unexpected shape of trajectory timeseries")

    @staticmethod
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def test_frame_interval_convention():
        universe1 = mda.Universe(PSF, DCD)
        array1 = universe1.trajectory.timeseries(skip=10)
        universe2 = mda.Universe(PSF, DCD, in_memory=True,
                                        in_memory_frame_interval=10)
        array2 = universe2.trajectory.timeseries()
        assert_equal(array1, array2,
                     err_msg="Unexpected differences between arrays.")


class TestCustomReaders(object):
    """
    Can pass a reader as kwarg on Universe creation
    """
    @dec.skipif(parser_not_found('TRZ'),
                'TRZ parser not available. Are you using python 3?')
    def test_custom_reader(self):
        # check that reader passing works
        u = mda.Universe(TRZ_psf, TRZ, format=MDAnalysis.coordinates.TRZ.TRZReader)
        assert_equal(len(u.atoms), 8184)

    def test_custom_reader_singleframe(self):
        T = MDAnalysis.topology.GROParser.GROParser
        R = MDAnalysis.coordinates.GRO.GROReader
        u = mda.Universe(two_water_gro, two_water_gro,
                                topology_format=T, format=R)
        assert_equal(len(u.atoms), 6)

    def test_custom_reader_singleframe_2(self):
        # Same as before, but only one argument to Universe
        T = MDAnalysis.topology.GROParser.GROParser
        R = MDAnalysis.coordinates.GRO.GROReader
        u = mda.Universe(two_water_gro,
                                topology_format=T, format=R)
        assert_equal(len(u.atoms), 6)

    @dec.skipif(parser_not_found('TRZ'),
                'TRZ parser not available. Are you using python 3?')
    def test_custom_parser(self):
        # topology reader passing works
        u = mda.Universe(TRZ_psf, TRZ, topology_format=MDAnalysis.topology.PSFParser.PSFParser)
        assert_equal(len(u.atoms), 8184)

    @dec.skipif(parser_not_found('TRZ'),
                'TRZ parser not available. Are you using python 3?')
    def test_custom_both(self):
        # use custom for both
        u = mda.Universe(TRZ_psf, TRZ, format=MDAnalysis.coordinates.TRZ.TRZReader,
                         topology_format=MDAnalysis.topology.PSFParser.PSFParser)
        assert_equal(len(u.atoms), 8184)


