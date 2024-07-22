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
import pickle

import os
import subprocess
import errno
from collections import defaultdict
from io import StringIO
import warnings

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
    PDB, MMTF,
)

import MDAnalysis as mda
import MDAnalysis.coordinates
from MDAnalysis.topology.base import TopologyReaderBase
from MDAnalysis.transformations import translate
from MDAnalysisTests import assert_nowarns
from MDAnalysis.exceptions import NoDataError
from MDAnalysis.core.topologyattrs import AtomStringAttr
from MDAnalysisTests.util import get_userid


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

    def test_Universe_invalidfile_IE_msg(self, tmpdir):
        # check for invalid file (something with the wrong content)
        with tmpdir.as_cwd():
            with open('invalid.file.tpr', 'w') as temp_file:
                temp_file.write('plop')
            try:
                mda.Universe('invalid.file.tpr')
            except IOError as e:
                assert 'file or cannot be recognized' in e.args[0]
            else:
                raise AssertionError

    @pytest.mark.skipif(get_userid() == 0,
                        reason="cannot permisssionerror as root")
    def test_Universe_invalidpermissionfile_IE_msg(self, tmpdir):
        # check for file with invalid permissions (eg. no read access)
        with tmpdir.as_cwd():
            temp_file = 'permission.denied.tpr'
            with open(temp_file, 'w'):
                pass

            if os.name == 'nt':
                subprocess.call("icacls {filename} /deny Users:RX".format(filename=temp_file),
                                shell=True)
            else:
                os.chmod(temp_file, 0o200)

            # Issue #3221 match by PermissionError and error number instead
            with pytest.raises(PermissionError, match=f"Errno {errno.EACCES}"):
                mda.Universe('permission.denied.tpr')

    def test_load_new_VE(self):
        u = mda.Universe.empty(0)

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

    def test_universe_empty_ROMol(self):
        Chem = pytest.importorskip("rdkit.Chem")
        mol = Chem.Mol()
        u = mda.Universe(mol, format="RDKIT")
        assert len(u.atoms) == 0


class TestUniverseFromSmiles(object):
    def setup_class(self):
        pytest.importorskip("rdkit.Chem")

    def test_default(self):
        smi = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
        u = mda.Universe.from_smiles(smi, format='RDKIT')
        assert u.atoms.n_atoms == 24
        assert len(u.bonds.indices) == 25

    def test_from_bad_smiles(self):
        with pytest.raises(SyntaxError) as e:
            u = mda.Universe.from_smiles("J", format='RDKIT')
            assert "Error while parsing SMILES" in str(e.value)

    def test_no_Hs(self):
        smi = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
        u = mda.Universe.from_smiles(smi, addHs=False, 
            generate_coordinates=False, format='RDKIT')
        assert u.atoms.n_atoms == 14
        assert len(u.bonds.indices) == 15

    def test_gencoords_without_Hs_error(self):
        with pytest.raises(ValueError) as e:
            u = mda.Universe.from_smiles("CCO", addHs=False,
                generate_coordinates=True, format='RDKIT')
            assert "requires adding hydrogens" in str (e.value)

    def test_generate_coordinates_numConfs(self):
        with pytest.raises(SyntaxError) as e:
            u = mda.Universe.from_smiles("CCO", numConfs=0, format='RDKIT')
            assert "non-zero positive integer" in str (e.value)
        with pytest.raises(SyntaxError) as e:
            u = mda.Universe.from_smiles("CCO", numConfs=2.1, format='RDKIT')
            assert "non-zero positive integer" in str (e.value)

    def test_rdkit_kwargs(self):
        # test for bad kwarg:
        # Unfortunately, exceptions from Boost cannot be passed to python,
        # we cannot `from Boost.Python import ArgumentError` and use it with 
        # pytest.raises(ArgumentError), so "this is the way"
        try:
            u = mda.Universe.from_smiles("CCO", rdkit_kwargs=dict(abc=42))
        except Exception as e:
            assert "did not match C++ signature" in str(e)
        else:
            raise AssertionError("RDKit should have raised an ArgumentError "
                                 "from Boost")
        # good kwarg
        u1 = mda.Universe.from_smiles("C", rdkit_kwargs=dict(randomSeed=42))
        u2 = mda.Universe.from_smiles("C", rdkit_kwargs=dict(randomSeed=51))
        with pytest.raises(AssertionError) as e:
            assert_equal(u1.trajectory.coordinate_array, 
                         u2.trajectory.coordinate_array)
            assert "Mismatched elements: 15 / 15 (100%)" in str(e.value)


    def test_coordinates(self):
        u = mda.Universe.from_smiles("C", numConfs=2, 
                                     rdkit_kwargs=dict(randomSeed=42))
        assert u.trajectory.n_frames == 2
        expected = np.array([
        [[-0.02209686,  0.00321505,  0.01651974],
         [-0.6664637 ,  0.8884155 , -0.10135844],
         [-0.37778795, -0.8577519 , -0.58829606],
         [ 0.09642092, -0.3151253 ,  1.0637809 ],
         [ 0.96992755,  0.2812466 , -0.39064613]],
        [[-0.0077073 ,  0.00435363,  0.01834692],
         [-0.61228824, -0.83705765, -0.38619974],
         [-0.41925883,  0.9689095 , -0.3415968 ],
         [ 0.03148226, -0.03256683,  1.1267245 ],
         [ 1.0077721 , -0.10363862, -0.41727486]]], dtype=np.float32)
        assert_allclose(u.trajectory.coordinate_array, expected, rtol=1e-5)


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

    def test_pickle(self):
        u = mda.Universe(PSF, DCD)
        s = pickle.dumps(u, protocol=pickle.HIGHEST_PROTOCOL)
        new_u = pickle.loads(s)
        assert_equal(u.atoms.names, new_u.atoms.names)


    @pytest.mark.parametrize('dtype', (int, np.float32, np.float64))
    def test_set_dimensions(self, dtype):
        u = mda.Universe(PSF, DCD)
        box = np.array([10, 11, 12, 90, 90, 90], dtype=dtype)
        u.dimensions = box
        assert_allclose(u.dimensions, box)
        assert u.dimensions.dtype == np.float32


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
    @pytest.fixture(scope='module')
    def vdw(self):
        return {'A': 1.4, 'B': 0.5}

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
            mda.Universe(two_water_gro_nonames, guess_bonds=True)

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

    def test_universe_guess_bonds_arguments(self):
        """Test if 'fudge_factor', and 'lower_bound' parameters
        are being passed correctly.
        """
        u = mda.Universe(two_water_gro, guess_bonds=True)
        
        self._check_universe(u)
        assert u.kwargs["guess_bonds"]
        assert u.kwargs["fudge_factor"]
        assert u.kwargs["lower_bound"]

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

    @pytest.mark.parametrize(
        'ff, lb, nbonds',
        [
            (0.55, 0.1, 2), (0.9, 1.6, 1),
            (0.5, 0.2, 2), (0.1, 0.1, 0)
        ]
    )
    def test_atomgroup_guess_bonds(self, ff, lb, nbonds):
        """Test an atomgroup doing guess bonds"""
        u = mda.Universe(two_water_gro)

        ag = u.atoms[:3]
        ag.guess_bonds(fudge_factor=ff, lower_bound=lb)
        # Apply '_check_atomgroup()' only in case of default values for
        # 'fudge_factor', and 'lower_bound' arguments to avoid unnecessary
        # errors.
        if ff == 0.55 and lb == 0.1:
            self._check_atomgroup(ag, u)
        assert len(ag.bonds) == nbonds

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

    def test_guess_bonds_periodicity(self):
        u = mda.Universe(two_water_gro)

        ag = u.atoms[:3]
        ag[0].position += u.dimensions[:3] * 10
        ag.guess_bonds()

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

    def test_transfer_to_memory_kwargs(self):
        u = mda.Universe(PSF, DCD)
        u.transfer_to_memory(example_kwarg=True)
        assert(u.trajectory._kwargs['example_kwarg'])

    def test_in_memory_kwargs(self):
        u = mda.Universe(PSF, DCD, in_memory=True, example_kwarg=True)
        assert(u.trajectory._kwargs['example_kwarg'])

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

    @pytest.mark.parametrize(
        'attr,values', (
            ('bonds', [(1, 0), (1, 2)]),
            ('bonds', [[1, 0], [1, 2]]),
            ('bonds', set([(1, 0), (1, 2)])),
            ('angles', [(1, 0, 2), (1, 2, 3), (2, 1, 4)]),
            ('dihedrals', [[1, 2, 3, 1], (3, 1, 5, 2)]),
            ('impropers', [[1, 2, 3, 1], (3, 1, 5, 2)]),
        )
    )
    def test_add_connection(self, universe, attr, values):
        universe.add_TopologyAttr(attr, values)
        assert hasattr(universe, attr)
        attrgroup = getattr(universe, attr)
        assert len(attrgroup) == len(values)
        for x in attrgroup:
            ix = x.indices
            assert ix[0] <= ix[-1]

    @pytest.mark.parametrize(
        'attr,values', (
            ('bonds', [(1, 0, 0), (1, 2)]),
            ('bonds', [['x', 'y'], [1, 2]]),
            ('bonds', 'rubbish'),
            ('bonds', [[1.01, 2.0]]),
            ('angles', [(1, 0), (1, 2)]),
            ('angles', 'rubbish'),
            ('dihedrals', [[1, 1, 1, 0.1]]),
            ('impropers', [(1, 2, 3)]),
        )
    )
    def add_connection_error(self, universe, attr, values):
        with pytest.raises(ValueError):
            universe.add_TopologyAttr(attr, values)


class TestDelTopologyAttr(object):
    @pytest.fixture()
    def universe(self):
        u = make_Universe(("masses", "charges", "resnames"))
        u.add_TopologyAttr("bonds", [(0, 1)])
        return u

    def test_del_TA_fail(self, universe):
        with pytest.raises(ValueError, match="Unrecognised"):
            universe.del_TopologyAttr('silly')

    def test_absent_fail(self, universe):
        with pytest.raises(ValueError, match="not in Universe"):
            universe.del_TopologyAttr("angles")

    def test_wrongtype_fail(self, universe):
        with pytest.raises(ValueError, match="must be str or TopologyAttr"):
            universe.del_TopologyAttr(list)

    @pytest.mark.parametrize(
        'todel,attrname', [
            ("charge", "charges"),
            ("charges", "charges"),
            ("bonds", "bonds"),
        ]
    )
    def test_del_str(self, universe, todel, attrname):
        assert hasattr(universe.atoms, attrname)
        universe.del_TopologyAttr(todel)
        assert not hasattr(universe.atoms, attrname)

    def test_del_attr(self, universe):
        assert hasattr(universe.atoms, "resnames")
        assert hasattr(universe.residues, "resnames")
        assert hasattr(universe.segments, "resnames")
        assert hasattr(universe.atoms[0], "resname")
        assert hasattr(universe.residues[0], "resname")
        universe.del_TopologyAttr(universe._topology.resnames)
        assert not hasattr(universe.atoms, "resnames")
        assert not hasattr(universe.residues, "resnames")
        assert not hasattr(universe.segments, "resnames")
        assert not hasattr(universe.atoms[0], "resname")
        assert not hasattr(universe.residues[0], "resname")

    def test_del_transplants(self, universe):
        atom = universe.atoms[0]
        assert hasattr(atom, "fragindex")
        universe.del_TopologyAttr(universe.bonds)
        assert not hasattr(atom, "fragindex")
        assert not hasattr(universe.atoms, "fragindices")

    def test_del_attr_error(self, universe):
        assert not hasattr(universe._topology, "elements")
        with pytest.raises(AttributeError):
            universe._topology.del_TopologyAttr("elements")

    def test_del_attr_from_ag(self, universe):
        ag = universe.atoms[[0]]
        ag.residues.resnames = "xyz"
        universe.del_TopologyAttr("resnames")
        with pytest.raises(NoDataError):
            ag.resnames

    def test_del_func_from_universe(self, universe):
        class RootVegetable(AtomStringAttr):
            attrname = "tubers"
            singular = "tuber"
            transplants = defaultdict(list)

            def potatoes(self):
                """🥔
                """
                return "potoooooooo"

            transplants["Universe"].append(("potatoes", potatoes))
        
        universe.add_TopologyAttr("tubers")
        assert universe.potatoes() == "potoooooooo"
        universe.del_TopologyAttr("tubers")
        with pytest.raises(AttributeError):
            universe.potatoes()


def _a_or_reversed_in_b(a, b):
    """
    Check if smaller array ``a`` or ``a[::-1]`` is in b
    """
    return (a==b).all(1).any() or (a[::-1]==b).all(1).any()

class TestAddTopologyObjects(object):

    small_atom_indices = (
        ('bonds', [[0, 1], [2, 3]]),
        ('angles', [[0, 1, 2], [3, 4, 5]]),
        ('dihedrals', [[8, 22, 1, 3], [4, 5, 6, 7], [11, 2, 3, 13]]),
        ('impropers', [[1, 6, 7, 2], [5, 3, 4, 2]]),
    )

    large_atom_indices = (
        ('bonds', [[0, 111], [22, 3]]),
        ('angles', [[0, 111, 2], [3, 44, 5]]),
        ('dihedrals', [[8, 222, 1, 3], [44, 5, 6, 7], [111, 2, 3, 13]]),
        ('impropers', [[1, 6, 771, 2], [5, 3, 433, 2]]),
    )

    @pytest.fixture()
    def empty(self):
        return make_Universe()

    @pytest.fixture()
    def universe(self):
        return mda.Universe(PSF)

    def _check_valid_added_to_empty(self, u, attr, values, to_add):
        assert not hasattr(u, attr)
        _add_func = getattr(u, 'add_'+attr)
        _add_func(to_add)
        u_attr = getattr(u, attr)
        assert len(u_attr) == len(values)
        assert all(_a_or_reversed_in_b(x, u_attr.indices)
                   for x in values)

    def _check_valid_added_to_populated(self, u, attr, values, to_add):
        assert hasattr(u, attr)
        u_attr = getattr(u, attr)
        original_length = len(u_attr)

        _add_func = getattr(u, 'add_'+attr)
        _add_func(to_add)
        u_attr = getattr(u, attr)
        assert len(u_attr) == len(values) + original_length
        assert all(_a_or_reversed_in_b(x, u_attr.indices)
                   for x in values)

    def _check_invalid_addition(self, u, attr, to_add, err_msg):
        _add_func = getattr(u, 'add_'+attr)
        with pytest.raises(ValueError) as excinfo:
            _add_func(to_add)
        assert err_msg in str(excinfo.value)

    @pytest.mark.parametrize(
        'attr,values', small_atom_indices
    )
    def test_add_indices_to_empty(self, empty, attr, values):
        self._check_valid_added_to_empty(empty, attr, values, values)

    def test_add_reversed_duplicates(self, empty):
        assert not hasattr(empty, 'bonds')
        empty.add_bonds([[0, 1], [1, 0]])
        assert len(empty.bonds) == 1
        assert_array_equal(empty.bonds.indices, np.array([[0, 1]]))

    @pytest.mark.parametrize(
        'attr,values', large_atom_indices
    )
    def test_add_indices_to_populated(self, universe, attr, values):
        self._check_valid_added_to_populated(universe, attr, values, values)

    @pytest.mark.parametrize(
        'attr,values', small_atom_indices
    )
    def test_add_atomgroup_to_empty(self, empty, attr, values):
        ag = [empty.atoms[x] for x in values]
        self._check_valid_added_to_empty(empty, attr, values, ag)

    @pytest.mark.parametrize(
        'attr,values', large_atom_indices
    )
    def test_add_atomgroup_to_populated(self, universe, attr, values):
        ag = [universe.atoms[x] for x in values]
        self._check_valid_added_to_populated(universe, attr, values, ag)

    @pytest.mark.parametrize(
        'attr,values', small_atom_indices
    )
    def test_add_atomgroup_wrong_universe_error(self, universe, empty, attr, values):
        ag = [empty.atoms[x] for x in values]
        self._check_invalid_addition(universe, attr, ag, 'different Universes')

    @pytest.mark.parametrize(
        'attr,values', large_atom_indices
    )
    def test_add_topologyobjects_to_populated(self, universe, attr, values):
        topologyobjects = [getattr(universe.atoms[x], attr[:-1]) for x in values]
        self._check_valid_added_to_populated(universe, attr, values, topologyobjects)

    @pytest.mark.parametrize(
        'attr,values', small_atom_indices
    )
    def test_add_topologyobjects_wrong_universe_error(self, universe, empty, attr, values):
        tobj = [getattr(universe.atoms[x], attr[:-1]) for x in values]
        self._check_invalid_addition(empty, attr, tobj, 'different Universes')

    @pytest.mark.parametrize(
        'attr,values', large_atom_indices
    )
    def test_add_topologygroups_to_populated(self, universe, attr, values):
        topologygroup = mda.core.topologyobjects.TopologyGroup(np.array(values),
                                                               universe)
        self._check_valid_added_to_populated(universe, attr, values, topologygroup)

    @pytest.mark.parametrize(
        'attr,values', small_atom_indices
    )
    def test_add_topologygroup_wrong_universe_error(self, universe, empty, attr, values):
        tg = mda.core.topologyobjects.TopologyGroup(np.array(values),
                                                    universe)
        self._check_invalid_addition(empty, attr, tg, 'different Universes')

    @pytest.mark.parametrize(
        'attr,values', small_atom_indices
    )
    def test_add_topologygroup_different_universe(self, universe, empty, attr, values):
        tg = mda.core.topologyobjects.TopologyGroup(np.array(values),
                                                               universe)
        self._check_valid_added_to_empty(empty, attr, values, tg.to_indices())

    @pytest.mark.parametrize(
        'attr,values', (
            ('impropers', [[0, 111], [22, 3]]),
            ('dihedrals', [[0, 111, 2], [3, 44, 5]]),
            ('angles', [[8, 222, 1, 3], [44, 5, 6, 7], [111, 2, 3, 13]]),
            ('bonds', [[1, 6, 771, 2], [5, 3, 433, 2]]),
        )
    )
    def test_add_wrong_topologygroup_error(self, universe, attr, values):
        arr = np.array(values)
        tg = mda.core.topologyobjects.TopologyGroup(arr, universe)
        self._check_invalid_addition(universe, attr, tg, 'iterable of tuples with')

    @pytest.mark.parametrize(
        'attr,values', (
            ('bonds', [[0, -111], [22, 3]]),
            ('angles', [[0, 11111, 2], [3, 44, 5]]),
            ('dihedrals', [[8, 222, 28888, 3], [44, 5, 6, 7], [111, 2, 3, 13]]),
            ('impropers', [[1, 6, 77133, 2], [5, 3, 433, 2]]),
        )
    )
    def test_add_nonexistent_indices_error(self, universe, attr, values):
        self._check_invalid_addition(universe, attr, values, 'nonexistent atom indices')

    @pytest.mark.parametrize(
        'attr,n', (
            ('bonds', 2),
            ('angles', 3),
            ('dihedrals', 4),
            ('impropers', 4),
        )
    )
    def test_add_wrong_number_of_atoms_error(self, universe, attr, n):
        errmsg = ('{} must be an iterable of '
                  'tuples with {} atom indices').format(attr, n)
        idx = [(0, 1), (0, 1, 2), (8, 22, 1, 3), (5, 3, 4, 2)]
        self._check_invalid_addition(universe, attr, idx, errmsg)

    def test_add_bonds_refresh_fragments(self, empty):
        with pytest.raises(NoDataError):
            getattr(empty.atoms, 'fragments')

        empty.add_bonds([empty.atoms[:2]])
        assert len(empty.atoms.fragments) == len(empty.atoms)-1

        empty.add_bonds([empty.atoms[2:4]])
        assert len(empty.atoms.fragments) == len(empty.atoms)-2

    @pytest.mark.parametrize(
        'attr,values', small_atom_indices
    )
    def test_roundtrip(self, empty, attr, values):
        _add_func = getattr(empty, 'add_'+attr)
        _add_func(values)
        u_attr = getattr(empty, attr)
        assert len(u_attr) == len(values)

        _delete_func = getattr(empty, 'delete_'+attr)
        _delete_func(values)
        u_attr = getattr(empty, attr)
        assert len(u_attr) == 0

class TestDeleteTopologyObjects(object):

    TOP = {'bonds': [(0, 1), (2, 3), (3, 4), (4, 5), (7, 8)],
           'angles': [(0, 1, 2), (3, 4, 5), (8, 2, 4)],
           'dihedrals': [(9, 2, 3, 4), (1, 3, 4, 2), (8, 22, 1, 3), (4, 5, 6, 7), (11, 2, 3, 13)],
           'impropers': [(1, 3, 5, 2), (1, 6, 7, 2), (5, 3, 4, 2)]}

    existing_atom_indices = (
            ('bonds', [[0, 1], [2, 3]]),
            ('angles', [[0, 1, 2], [3, 4, 5]]),
            ('dihedrals', [[8, 22, 1, 3], [4, 5, 6, 7], [11, 2, 3, 13]]),
            ('impropers', [[1, 6, 7, 2], [5, 3, 4, 2]]),
        )
    nonexisting_atom_indices = (
            ('bonds', [[2, 3], [7, 8], [0, 4]]),
            ('angles', [[0, 1, 2], [8, 2, 8], [1, 1, 1]]),
            ('dihedrals', [[0, 0, 0, 0], [1, 1, 1, 1]]),
            ('impropers', [[8, 22, 1, 3],]),
        )

    @pytest.fixture()
    def universe(self):
        u = make_Universe(size=(125, 25, 5))
        for attr, values in self.TOP.items():
            u._add_topology_objects(attr, values)
        return u

    @pytest.fixture()
    def universe2(self):
        u = make_Universe(size=(125, 25, 5))
        for attr, values in self.TOP.items():
            u._add_topology_objects(attr, values)
        return u

    def _check_valid_deleted(self, u, attr, values, to_delete):
        u_attr = getattr(u, attr)
        original_length = len(self.TOP[attr])
        assert len(u_attr) == original_length

        _delete_func = getattr(u, 'delete_'+attr)
        _delete_func(to_delete)
        u_attr = getattr(u, attr)
        assert len(u_attr) == original_length-len(values)

        not_deleted = [x for x in self.TOP[attr] if list(x) not in values]
        assert all([x in u_attr.indices or x[::-1] in u_attr.indices
                   for x in not_deleted])

    def _check_invalid_deleted(self, u, attr, to_delete, err_msg):
        u_attr = getattr(u, attr)
        original_length = len(self.TOP[attr])
        assert len(u_attr) == original_length
        _delete_func = getattr(u, 'delete_'+attr)
        with pytest.raises(ValueError) as excinfo:
            _delete_func(to_delete)
        assert err_msg in str(excinfo.value)

    @pytest.mark.parametrize(
        'attr,values', existing_atom_indices
    )
    def test_delete_valid_indices(self, universe, attr, values):
        self._check_valid_deleted(universe, attr, values, values)

    @pytest.mark.parametrize(
        'attr,values', nonexisting_atom_indices
    )
    def test_delete_missing_indices(self, universe, attr, values):
        self._check_invalid_deleted(universe, attr, values, 'Cannot delete nonexistent')

    @pytest.mark.parametrize(
        'attr,values', existing_atom_indices
    )
    def test_delete_valid_atomgroup(self, universe, attr, values):
        ag = [universe.atoms[x] for x in values]
        self._check_valid_deleted(universe, attr, values, ag)

    @pytest.mark.parametrize(
        'attr,values', existing_atom_indices
    )
    def test_delete_atomgroup_wrong_universe_error(self, universe, universe2, attr, values):
        ag = [universe.atoms[x] for x in values]
        self._check_invalid_deleted(universe2, attr, ag, 'different Universes')

    @pytest.mark.parametrize(
        'attr,values', nonexisting_atom_indices
    )
    def test_delete_missing_atomgroup(self, universe, attr, values):
        ag = [universe.atoms[x] for x in values]
        self._check_invalid_deleted(universe, attr, ag, 'Cannot delete nonexistent')

    @pytest.mark.parametrize(
        'attr,values', existing_atom_indices
    )
    def test_delete_mixed_type(self, universe, attr, values):
        mixed = [universe.atoms[values[0]]] + values[1:]
        self._check_valid_deleted(universe, attr, values, mixed)

    @pytest.mark.parametrize(
        'attr,values', existing_atom_indices
    )
    def test_delete_valid_topologyobjects(self, universe, attr, values):
        to = [getattr(universe.atoms[x], attr[:-1]) for x in values]
        self._check_valid_deleted(universe, attr, values, to)

    @pytest.mark.parametrize(
        'attr,values', existing_atom_indices
    )
    def test_delete_topologyobjects_wrong_universe(self, universe, universe2, attr, values):
        u1 = [getattr(universe.atoms[x], attr[:-1]) for x in values[:-1]]
        u2 = [getattr(universe2.atoms[values[-1]], attr[:-1])]
        self._check_invalid_deleted(universe, attr, u1+u2, 'different Universes')

    @pytest.mark.parametrize(
        'attr,values', existing_atom_indices
    )
    def test_delete_valid_topologygroup(self, universe, attr, values):
        arr = np.array(values)
        tg = mda.core.topologyobjects.TopologyGroup(arr, universe)
        self._check_valid_deleted(universe, attr, values, tg)

    @pytest.mark.parametrize(
        'attr,values', existing_atom_indices
    )
    def test_delete_topologygroup_wrong_universe_error(self, universe, universe2, attr, values):
        arr = np.array(values)
        tg = mda.core.topologyobjects.TopologyGroup(arr, universe2)
        self._check_invalid_deleted(universe, attr, tg, 'different Universes')

    @pytest.mark.parametrize(
        'attr,values', existing_atom_indices
    )
    def test_delete_topologygroup_different_universe(self, universe, universe2, attr, values):
        arr = np.array(values)
        tg = mda.core.topologyobjects.TopologyGroup(arr, universe2)
        self._check_valid_deleted(universe, attr, values, tg.to_indices())

    @pytest.mark.parametrize(
        'attr,n', (
            ('bonds', 2),
            ('angles', 3),
            ('dihedrals', 4),
            ('impropers', 4),
        )
    )
    def test_delete_wrong_number_of_atoms_error(self, universe, attr, n):
        idx = [(0, 1), (0, 1, 2), (8, 22, 1, 3), (5, 3, 4, 2)]
        errmsg = ('{} must be an iterable of '
                  'tuples with {} atom indices').format(attr, n)
        self._check_invalid_deleted(universe, attr, idx, errmsg)

    @pytest.mark.parametrize(
        'attr,values', existing_atom_indices
    )
    def test_delete_missing_attr(self, attr, values):
        u = make_Universe()
        assert not hasattr(u, attr)
        _delete_func = getattr(u, 'delete_'+attr)
        with pytest.raises(ValueError) as excinfo:
            _delete_func(values)
        assert "There are no" in str(excinfo.value)

    def test_delete_bonds_refresh_fragments(self, universe):
        n_fragments = len(universe.atoms.fragments)
        universe.delete_bonds([universe.atoms[[2, 3]]])
        assert len(universe.atoms.fragments) == n_fragments + 1

    @pytest.mark.parametrize(
        'attr,values', existing_atom_indices
    )
    def test_roundtrip(self, universe, attr, values):
        u_attr = getattr(universe, attr)
        original_length = len(self.TOP[attr])
        assert len(u_attr) == original_length

        _delete_func = getattr(universe, 'delete_'+attr)
        _delete_func(values)
        nu_attr = getattr(universe, attr)
        assert len(nu_attr) == original_length-len(values)

        _add_func = getattr(universe, 'add_'+attr)
        _add_func(values)
        nu_attr = getattr(universe, attr)
        assert len(nu_attr) == original_length

        assert_array_equal(u_attr.indices, nu_attr.indices)


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
            u = mda.Universe.empty(n_atoms=10, n_residues=2, n_segments=2,
                                   atom_resindex=res)

    def test_no_trivial_warning(self):
        """
        Make sure that no warning is raised about atom_resindex and
        residue_segindex when n_residues or n_segments is equal to 1.
        """
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            u = mda.Universe.empty(n_atoms=10, n_residues=1, n_segments=1)

    def test_trajectory(self):
        u = mda.Universe.empty(10, trajectory=True)

        assert len(u.atoms) == 10
        assert u.atoms.positions.shape == (10, 3)

    def test_trajectory_multiple_frames(self):
        u = mda.Universe.empty(10, n_frames=5)

        assert_equal(u.trajectory.n_frames, 5)

    def test_trajectory_iteration(self):
        u = mda.Universe.empty(10, trajectory=True)

        assert len(u.trajectory) == 1
        timesteps =[]
        for ts in u.trajectory:
            timesteps.append(ts.frame)
        assert len(timesteps) == 1

    def test_velocities(self):
        u = mda.Universe.empty(10, trajectory=True, velocities=True)

        assert u.atoms.positions.shape == (10, 3)
        assert u.atoms.velocities.shape == (10, 3)

    def test_forces(self):
        u = mda.Universe.empty(10, trajectory=True, forces=True)

        assert u.atoms.positions.shape == (10, 3)
        assert u.atoms.forces.shape == (10, 3)

    def test_empty_no_atoms(self):
        u = mda.Universe.empty(0)
        assert len(u.atoms) == 0
        assert len(u.residues) == 0
        assert len(u.segments) == 0

    def test_empty_creation_raises_error(self):
        with pytest.raises(TypeError) as exc:
            u = mda.Universe()
        assert 'Universe.empty' in str(exc.value)


def test_deprecate_b_tempfactors():
    u = mda.Universe(PDB)
    values = np.arange(len(u.atoms))
    with pytest.warns(DeprecationWarning, match="use the tempfactor"):
        u.add_TopologyAttr("bfactors", values)
    assert_array_equal(u.atoms.tempfactors, values)


class Thingy:
    def __init__(self, val):
        self.v = val


class ThingyParser(TopologyReaderBase):
    format='THINGY'

    @staticmethod
    def _format_hint(thing):
        return isinstance(thing, Thingy)

    def parse(self, **kwargs):
        return mda.core.topology.Topology(n_atoms=10)


class TestOnlyTopology:
    def test_only_top(self):
        # issue 3443
        t = Thingy(20)

        with pytest.warns(UserWarning,
                          match="No coordinate reader found for"):
            u = mda.Universe(t)

        assert len(u.atoms) == 10
