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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import absolute_import, print_function


from six.moves import cPickle

import os
import subprocess

try:
    from cStringIO import StringIO
except:
    from io import StringIO
from MDAnalysisTests.tempdir import TempDir

import numpy as np
from numpy.testing import (
    assert_allclose,
    assert_almost_equal,
    assert_equal,
    assert_array_equal,
)
import pytest

from MDAnalysisTests import make_Universe
from MDAnalysisTests.datafiles import (
    PSF, DCD,
    PSF_BAD,
    PDB_small,
    PDB_chainidrepeat,
    GRO, TRR,
    two_water_gro, two_water_gro_nonames,
    TRZ, TRZ_psf,
)

import MDAnalysis as mda
import MDAnalysis.coordinates
from MDAnalysis.topology.base import TopologyReaderBase
from MDAnalysis.transformations import translate
from MDAnalysisTests import assert_nowarns


class IOErrorParser(TopologyReaderBase):
    def parse(self, **kwargs):
        raise IOError("Useful information")

# This string is not in the `TestUniverseCreation` class or its method because of problems
# with whitespace. Extra indentations make the string unreadable.
CHOL_GRO = """\
Single cholesterol molecule
8
  153CHOL   ROH 1793   6.558   2.936   4.005 -0.1044 -0.1252  0.0966
  153CHOL    R1 1794   6.591   2.999   4.279  0.0998  0.1297  0.0158
  153CHOL    R2 1795   6.657   2.810   4.469  0.0780  0.2817  0.1592
  153CHOL    R3 1796   6.859   2.983   4.524  0.2233  0.2112 -0.1215
  153CHOL    R4 1797   6.804   2.849   4.779  0.4156  0.3232  0.0001
  153CHOL    R5 1798   6.810   3.064   4.744  0.4811  0.3182 -0.0905
  153CHOL    C1 1799   7.080   3.034   5.012  0.7486  0.2811 -0.3559
  153CHOL    C2 1800   6.993   3.163   5.284  0.3677 -0.2104 -0.0829
10 10 10
"""

class TestUniverseCreation(object):
    # tests concerning Universe creation and errors encountered
    def test_load(self):
        # Universe(top, trj)
        u = mda.Universe(PSF, PDB_small)
        assert_equal(len(u.atoms), 3341, "Loading universe failed somehow")

    def test_load_topology_stringio(self):
        u = mda.Universe(StringIO(CHOL_GRO), format='GRO')
        assert_equal(len(u.atoms), 8, "Loading universe from StringIO failed somehow")
        assert_equal(u.trajectory.ts.positions[0], np.array([65.580002, 29.360001, 40.050003], dtype=np.float32))

    def test_load_trajectory_stringio(self):
        u = mda.Universe(StringIO(CHOL_GRO), StringIO(CHOL_GRO),  format='GRO', topology_format='GRO')
        assert_equal(len(u.atoms), 8, "Loading universe from StringIO failed somehow")

    def test_make_universe_no_args(self):
        # universe creation without args should work
        u = mda.Universe()

        assert isinstance(u, mda.Universe)
        assert u.atoms is None

    def test_make_universe_stringio_no_format(self):
        # Loading from StringIO without format arg should raise TypeError
        with pytest.raises(TypeError):
            mda.Universe(StringIO(CHOL_GRO))

    def test_Universe_no_trajectory_AE(self):
        # querying trajectory without a trajectory loaded (only topology)
        u = make_Universe()
        with pytest.raises(AttributeError):
            getattr(u, 'trajectory')

    def test_Universe_topology_unrecognizedformat_VE(self):
        with pytest.raises(ValueError):
            mda.Universe('some.file.without.parser_or_coordinate_extension')

    def test_Universe_topology_unrecognizedformat_VE_msg(self):
        try:
            mda.Universe('some.file.without.parser_or_coordinate_extension')
        except ValueError as e:
            assert 'isn\'t a valid topology format' in e.args[0]
        else:
            raise AssertionError

    def test_Universe_topology_IE(self):
        with pytest.raises(IOError):
            mda.Universe('thisfile', topology_format = IOErrorParser)

    def test_Universe_topology_IE_msg(self):
        # should get the original error, as well as Universe error
        try:
            mda.Universe('thisfile', topology_format=IOErrorParser)
        except IOError as e:
            assert 'Failed to load from the topology file' in e.args[0]
            assert 'Useful information' in e.args[0]
        else:
            raise AssertionError

    def test_Universe_filename_IE_msg(self):
        # check for non existent file
        try:
            mda.Universe('thisfile.xml')
        except IOError as e:
            assert_equal('No such file or directory', e.strerror)
        else:
            raise AssertionError

    def test_Universe_invalidfile_IE_msg(self):
        # check for invalid file (something with the wrong content)
        temp_dir = TempDir()
        with open(os.path.join(temp_dir.name, 'invalid.file.tpr'), 'w') as temp_file:
            temp_file.write('plop')
        try:
            mda.Universe(os.path.join(temp_dir.name, 'invalid.file.tpr'))
        except IOError as e:
            assert 'file or cannot be recognized' in e.args[0]
        else:
            raise AssertionError
        finally:
            temp_dir.dissolve()

    def test_Universe_invalidpermissionfile_IE_msg(self):
        # check for file with invalid permissions (eg. no read access)
        temp_dir = TempDir()
        temp_file = os.path.join(temp_dir.name, 'permission.denied.tpr')
        with open(temp_file, 'w'):
            pass

        if os.name == 'nt':
            subprocess.call("icacls {filename} /deny Users:RX".format(filename=temp_file),
                            shell=True)
        else:
            os.chmod(temp_file, 0o200)
        try:
            mda.Universe(os.path.join(temp_dir.name, 'permission.denied.tpr'))
        except IOError as e:
            assert 'Permission denied' in str(e.strerror)
        else:
            raise AssertionError
        finally:
            temp_dir.dissolve()

    def test_load_new_VE(self):
        u = mda.Universe()

        with pytest.raises(TypeError):
            u.load_new('thisfile', format = 'soup')

    def test_load_new_memory_reader_success(self):
        u = mda.Universe(GRO)
        prot = u.select_atoms('protein')
        u2 = mda.Merge(prot)
        assert u2.load_new( [ prot.positions ], format=mda.coordinates.memory.MemoryReader) is u2

    def test_load_new_memory_reader_fails(self):
        def load():
            u = mda.Universe(GRO)
            prot = u.select_atoms('protein')
            u2 = mda.Merge(prot)
            u2.load_new( [[ prot.positions ]], format=mda.coordinates.memory.MemoryReader)

        with pytest.raises(TypeError):
            load()

    def test_universe_kwargs(self):
        u = mda.Universe(PSF, PDB_small, fake_kwarg=True)
        assert_equal(len(u.atoms), 3341, "Loading universe failed somehow")

        assert u.kwargs['fake_kwarg']

        # initialize new universe from pieces of existing one
        u2 = mda.Universe(u.filename, u.trajectory.filename,
                          **u.kwargs)

        assert u2.kwargs['fake_kwarg']
        assert_equal(u.kwargs, u2.kwargs)

    def test_universe_topology_class_with_coords(self):
        u = mda.Universe(PSF, PDB_small)
        u2 = mda.Universe(u._topology, PDB_small)
        assert isinstance(u2.trajectory, type(u.trajectory))
        assert_equal(u.trajectory.n_frames, u2.trajectory.n_frames)
        assert u2._topology is u._topology


class TestUniverse(object):
    # older tests, still useful
    def test_load_bad_topology(self):
        # tests that Universe builds produce the right error message
        def bad_load():
            return mda.Universe(PSF_BAD, DCD)

        with pytest.raises(ValueError):
            bad_load()

    def test_load_new(self):
        u = mda.Universe(PSF, DCD)
        u.load_new(PDB_small)
        assert_equal(len(u.trajectory), 1, "Failed to load_new(PDB)")

    def test_load_new_returns_Universe(self):
        u = mda.Universe(PSF)
        result = u.load_new(PDB_small)
        assert result is u

    def test_load_new_None_returns_Universe(self):
        u = mda.Universe(PSF)
        result = u.load_new(None)
        assert result is u

    def test_load_new_TypeError(self):
        u = mda.Universe(PSF, DCD)

        def bad_load(uni):
            return uni.load_new('filename.notarealextension')

        with pytest.raises(TypeError):
            bad_load(u)

    def test_load_structure(self):
        # Universe(struct)
        ref = mda.Universe(PSF, PDB_small)
        u = mda.Universe(PDB_small)
        assert_equal(len(u.atoms), 3341, "Loading universe failed somehow")
        assert_almost_equal(u.atoms.positions, ref.atoms.positions)

    def test_load_multiple_list(self):
        # Universe(top, [trj, trj, ...])
        ref = mda.Universe(PSF, DCD)
        u = mda.Universe(PSF, [DCD, DCD])
        assert_equal(len(u.atoms), 3341, "Loading universe failed somehow")
        assert_equal(u.trajectory.n_frames, 2 * ref.trajectory.n_frames)

    def test_load_multiple_args(self):
        # Universe(top, trj, trj, ...)
        ref = mda.Universe(PSF, DCD)
        u = mda.Universe(PSF, DCD, DCD)
        assert_equal(len(u.atoms), 3341, "Loading universe failed somehow")
        assert_equal(u.trajectory.n_frames, 2 * ref.trajectory.n_frames)

    def test_pickle_raises_NotImplementedError(self):
        u = mda.Universe(PSF, DCD)
        with pytest.raises(NotImplementedError):
            cPickle.dumps(u, protocol = cPickle.HIGHEST_PROTOCOL)

    def test_set_dimensions(self):
        u = mda.Universe(PSF, DCD)
        box = np.array([10, 11, 12, 90, 90, 90])
        u.dimensions = np.array([10, 11, 12, 90, 90, 90])
        assert_allclose(u.dimensions, box)


# remove for 1.0
def test_chainid_quick_select():
    # check that chainIDs get grouped together when making the quick selectors
    # this pdb file has 2 segments with chainID A
    u = mda.Universe(PDB_chainidrepeat)

    with pytest.warns(DeprecationWarning):
        for sg in (u.A, u.B):
            assert isinstance(sg, mda.core.groups.SegmentGroup)
        for seg in (u.C, u.D):
            assert isinstance(seg, mda.core.groups.Segment)
        assert len(u.A.atoms) == 10
        assert len(u.B.atoms) == 10
        assert len(u.C.atoms) == 5
        assert len(u.D.atoms) == 7

class TestTransformations(object):
    """Tests the transformations keyword
    """
    def test_callable(self):
        u = mda.Universe(PSF,DCD, transformations=translate([10,10,10]))
        uref = mda.Universe(PSF,DCD)
        ref = translate([10,10,10])(uref.trajectory.ts)
        assert_almost_equal(u.trajectory.ts.positions, ref, decimal=6)

    def test_list(self):
        workflow = [translate([10,10,0]), translate([0,0,10])]
        u = mda.Universe(PSF,DCD, transformations=workflow)
        uref = mda.Universe(PSF,DCD)
        ref = translate([10,10,10])(uref.trajectory.ts)
        assert_almost_equal(u.trajectory.ts.positions, ref, decimal=6)

class TestGuessMasses(object):
    """Tests the Mass Guesser in topology.guessers
    """
    def test_universe_loading_no_warning(self):
        assert_nowarns(UserWarning, lambda x: mda.Universe(x), GRO)


class TestGuessBonds(object):
    """Test the AtomGroup methed guess_bonds

    This needs to be done both from Universe creation (via kwarg) and AtomGroup

    It needs to:
     - work if all atoms are in vdwradii table
     - fail properly if not
     - work again if vdwradii are passed.
    """
    @pytest.fixture()
    def vdw(self):
        return {'A': 1.05, 'B': 0.4}

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
        assert 'guess_bonds' in u.kwargs

    def test_universe_guess_bonds(self):
        """Test that making a Universe with guess_bonds works"""
        u = mda.Universe(two_water_gro, guess_bonds=True)
        self._check_universe(u)
        assert u.kwargs['guess_bonds']

    def test_universe_guess_bonds_no_vdwradii(self):
        """Make a Universe that has atoms with unknown vdwradii."""
        with pytest.raises(ValueError):
            mda.Universe(two_water_gro_nonames, guess_bonds = True)

    def test_universe_guess_bonds_with_vdwradii(self, vdw):
        """Unknown atom types, but with vdw radii here to save the day"""
        u = mda.Universe(two_water_gro_nonames, guess_bonds=True,
                                vdwradii=vdw)
        self._check_universe(u)
        assert u.kwargs['guess_bonds']
        assert_equal(vdw, u.kwargs['vdwradii'])

    def test_universe_guess_bonds_off(self):
        u = mda.Universe(two_water_gro_nonames, guess_bonds=False)

        for attr in ('bonds', 'angles', 'dihedrals'):
            assert not hasattr(u, attr)
        assert not u.kwargs['guess_bonds']

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
        with pytest.raises(ValueError):
            ag.guess_bonds()

    def test_atomgroup_guess_bonds_with_vdwradii(self, vdw):
        u = mda.Universe(two_water_gro_nonames)

        ag = u.atoms[:3]
        ag.guess_bonds(vdwradii=vdw)
        self._check_atomgroup(ag, u)


class TestInMemoryUniverse(object):
    def test_reader_w_timeseries(self):
        universe = mda.Universe(PSF, DCD, in_memory=True)
        assert_equal(universe.trajectory.timeseries(universe.atoms).shape,
                     (3341, 98, 3),
                     err_msg="Unexpected shape of trajectory timeseries")

    def test_reader_wo_timeseries(self):
        universe = mda.Universe(GRO, TRR, in_memory=True)
        assert_equal(universe.trajectory.timeseries(universe.atoms).shape,
                     (47681, 10, 3),
                     err_msg="Unexpected shape of trajectory timeseries")

    def test_reader_w_timeseries_frame_interval(self):
        universe = mda.Universe(PSF, DCD, in_memory=True,
                                       in_memory_step=10)
        assert_equal(universe.trajectory.timeseries(universe.atoms).shape,
                     (3341, 10, 3),
                     err_msg="Unexpected shape of trajectory timeseries")

    def test_reader_wo_timeseries_frame_interval(self):
        universe = mda.Universe(GRO, TRR, in_memory=True,
                                       in_memory_step=3)
        assert_equal(universe.trajectory.timeseries(universe.atoms).shape,
                     (47681, 4, 3),
                     err_msg="Unexpected shape of trajectory timeseries")

    def test_existing_universe(self):
        universe = mda.Universe(PDB_small, DCD)
        universe.transfer_to_memory()
        assert_equal(universe.trajectory.timeseries(universe.atoms).shape,
                     (3341, 98, 3),
                     err_msg="Unexpected shape of trajectory timeseries")

    def test_frame_interval_convention(self):
        universe1 = mda.Universe(PSF, DCD)
        array1 = universe1.trajectory.timeseries(step=10)
        universe2 = mda.Universe(PSF, DCD, in_memory=True,
                                 in_memory_step=10)
        array2 = universe2.trajectory.timeseries()
        assert_equal(array1, array2,
                     err_msg="Unexpected differences between arrays.")

    def test_slicing_with_start_stop(self):
        universe = MDAnalysis.Universe(PDB_small, DCD)
        # Skip only the last frame
        universe.transfer_to_memory(start=10, stop=20)
        assert_equal(universe.trajectory.timeseries(universe.atoms).shape,
                     (3341, 10, 3),
                     err_msg="Unexpected shape of trajectory timeseries")

    def test_slicing_without_start(self):
        universe = MDAnalysis.Universe(PDB_small, DCD)
        # Skip only the last frame
        universe.transfer_to_memory(stop=10)
        assert_equal(universe.trajectory.timeseries(universe.atoms).shape,
                     (3341, 10, 3),
                     err_msg="Unexpected shape of trajectory timeseries")

    def test_slicing_without_stop(self):
        universe = MDAnalysis.Universe(PDB_small, DCD)
        # Skip only the last frame
        universe.transfer_to_memory(start=10)
        print(universe.trajectory.timeseries(universe.atoms).shape)
        assert_equal(universe.trajectory.timeseries(universe.atoms).shape,
                     (3341, 88, 3),
                     err_msg="Unexpected shape of trajectory timeseries")

    def test_slicing_step_without_start_stop(self):
        universe = MDAnalysis.Universe(PDB_small, DCD)
        # Skip only the last frame
        universe.transfer_to_memory(step=2)
        print(universe.trajectory.timeseries(universe.atoms).shape)
        assert_equal(universe.trajectory.timeseries(universe.atoms).shape,
                     (3341, 49, 3),
                     err_msg="Unexpected shape of trajectory timeseries")

    def test_slicing_step_with_start_stop(self):
        universe = MDAnalysis.Universe(PDB_small, DCD)
        # Skip only the last frame
        universe.transfer_to_memory(start=10, stop=30, step=2)
        print(universe.trajectory.timeseries(universe.atoms).shape)
        assert_equal(universe.trajectory.timeseries(universe.atoms).shape,
                     (3341, 10, 3),
                     err_msg="Unexpected shape of trajectory timeseries")

    def test_slicing_step_dt(self):
        universe = MDAnalysis.Universe(PDB_small, DCD)
        dt = universe.trajectory.dt
        universe.transfer_to_memory(step=2)
        assert_almost_equal(dt * 2, universe.trajectory.dt,
                            err_msg="Unexpected in-memory timestep: "
                            + "dt not updated with step information")

    def test_slicing_negative_start(self):
        universe = MDAnalysis.Universe(PDB_small, DCD)
        # Skip only the last frame
        universe.transfer_to_memory(start=-10)
        print(universe.trajectory.timeseries(universe.atoms).shape)
        assert_equal(universe.trajectory.timeseries(universe.atoms).shape,
                     (3341, 10, 3),
                     err_msg="Unexpected shape of trajectory timeseries")

    def test_slicing_negative_stop(self):
        universe = MDAnalysis.Universe(PDB_small, DCD)
        # Skip only the last frame
        universe.transfer_to_memory(stop=-20)
        print(universe.trajectory.timeseries(universe.atoms).shape)
        assert_equal(universe.trajectory.timeseries(universe.atoms).shape,
                     (3341, 78, 3),
                     err_msg="Unexpected shape of trajectory timeseries")


class TestCustomReaders(object):
    """
    Can pass a reader as kwarg on Universe creation
    """
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

    def test_custom_parser(self):
        # topology reader passing works
        u = mda.Universe(TRZ_psf, TRZ, topology_format=MDAnalysis.topology.PSFParser.PSFParser)
        assert_equal(len(u.atoms), 8184)

    def test_custom_both(self):
        # use custom for both
        u = mda.Universe(TRZ_psf, TRZ, format=MDAnalysis.coordinates.TRZ.TRZReader,
                         topology_format=MDAnalysis.topology.PSFParser.PSFParser)
        assert_equal(len(u.atoms), 8184)


class TestAddTopologyAttr(object):
    @pytest.fixture()
    def universe(self):
        return make_Universe()

    def test_add_TA_fail(self, universe):
        with pytest.raises(ValueError):
            universe.add_TopologyAttr('silly')

    def test_nodefault_fail(self, universe):
        with pytest.raises(NotImplementedError):
            universe.add_TopologyAttr('bonds')

    @pytest.mark.parametrize(
        'toadd,attrname,default', (
            ['charge', 'charges', 0.0], ['charges', 'charges', 0.0],
            ['name', 'names', ''], ['names', 'names', ''],
            ['type', 'types', ''], ['types', 'types', ''],
            ['element', 'elements', ''], ['elements', 'elements', ''],
            ['radius', 'radii', 0.0], ['radii', 'radii', 0.0],
            ['chainID', 'chainIDs', ''], ['chainIDs', 'chainIDs', ''],
            ['tempfactor', 'tempfactors', 0.0],
            ['tempfactors', 'tempfactors', 0.0],
            ['mass', 'masses', 0.0], ['masses', 'masses', 0.0],
            ['charge', 'charges', 0.0], ['charges', 'charges', 0.0],
            ['bfactor', 'bfactors', 0.0], ['bfactors', 'bfactors', 0.0],
            ['occupancy', 'occupancies', 0.0],
            ['occupancies', 'occupancies', 0.0],
            ['altLoc', 'altLocs', ''], ['altLocs', 'altLocs', ''],
            ['resid', 'resids', 1], ['resids', 'resids', 1],
            ['resname', 'resnames', ''], ['resnames', 'resnames', ''],
            ['resnum', 'resnums', 1], ['resnums', 'resnums', 1],
            ['icode', 'icodes', ''], ['icodes', 'icodes', ''],
            ['segid', 'segids', ''], ['segids', 'segids', ''],
        )
    )
    def test_add_charges(self, universe, toadd, attrname, default):
        universe.add_TopologyAttr(toadd)

        assert hasattr(universe.atoms, attrname)
        assert getattr(universe.atoms, attrname)[0] == default


class TestAllCoordinatesKwarg(object):
    @pytest.fixture(scope='class')
    def u_GRO_TRR(self):
        return mda.Universe(GRO, TRR)

    @pytest.fixture(scope='class')
    def u_GRO_TRR_allcoords(self):
        return mda.Universe(GRO, TRR, all_coordinates=True)

    @pytest.fixture(scope='class')
    def u_GRO(self):
        return mda.Universe(GRO)

    def test_all_coordinates_length(self, u_GRO_TRR, u_GRO_TRR_allcoords):
        # length with all_coords should be +1
        assert (len(u_GRO_TRR.trajectory) + 1 ==
                len(u_GRO_TRR_allcoords.trajectory))

    def test_all_coordinates_frame(self, u_GRO_TRR_allcoords, u_GRO):
        # check that first frame in u(GRO, TRR, allcords)
        # are the coordinates from GRO
        assert_array_equal(u_GRO_TRR_allcoords.atoms.positions,
                           u_GRO.atoms.positions)

    def test_second_frame(self, u_GRO_TRR_allcoords, u_GRO_TRR):
        # check that second frame in u(GRO, TRR, allcoords)
        # are the coordinates from TRR[0]
        u_GRO_TRR_allcoords.trajectory[1]

        assert_array_equal(u_GRO_TRR_allcoords.atoms.positions,
                           u_GRO_TRR.atoms.positions)


class TestEmpty(object):
    def test_empty(self):
        u = mda.Universe.empty(10)

        assert len(u.atoms) == 10
        assert len(u.residues) == 1
        assert len(u.segments) == 1

    def test_empty_extra(self):
        u = mda.Universe.empty(
            n_atoms=12, n_residues=3, n_segments=2,
            atom_resindex=np.array([0, 0, 0, 0, 0,
                                    1, 1, 1, 1, 1,
                                    2, 2]),
            residue_segindex=np.array([0, 0, 1]),
        )

        assert len(u.atoms) == 12

        assert len(u.residues) == 3
        assert len(u.residues[0].atoms) == 5
        assert len(u.residues[1].atoms) == 5
        assert len(u.residues[2].atoms) == 2

        assert len(u.segments) == 2
        assert len(u.segments[0].atoms) == 10
        assert len(u.segments[1].atoms) == 2

    def test_no_resindex_warning(self):
        with pytest.warns(UserWarning):
            u = mda.Universe.empty(n_atoms=10, n_residues=2, n_segments=1)

    def test_no_segindex_warning(self):
        res = np.array([0, 0, 0, 0, 0,
                        1, 1, 1, 1, 1])

        with pytest.warns(UserWarning):
            u = mda.Universe.empty(n_atoms=10, n_residues=2, n_segments=1,
                                   atom_resindex=res)
