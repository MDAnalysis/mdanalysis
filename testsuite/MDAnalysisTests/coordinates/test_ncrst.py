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
import MDAnalysis as mda
import numpy as np
from scipy.io import netcdf
from MDAnalysis.coordinates.INPCRD import NCRSTReader
import pytest
from numpy.testing import (
    assert_equal,
    assert_almost_equal
)

from MDAnalysisTests.datafiles import (PFncdf_Top, PFncdf_Trj,
                                       GRO, TRR, XYZ_mini)
from MDAnalysisTests.coordinates.reference import (RefVGV, RefTZ2)


class _NCRestartTest(object):

    prec = 5

    @pytest.fixture(scope='class')
    def universe(self):
        return mda.Universe(self.topology, self.filename, mmap=self.mmap)

    def test_wrong_natoms(self):
        with pytest.raises(ValueError):
            NCRSTReader(self.filename, n_atoms=2, mmap=self.mmap)

    def test_mmap_value(self):
        reader = NCRSTReader(self.filename, mmap=self.mmap)
        if not self.mmap:
            assert_equal(reader._mmap, None)
        else:
            assert_equal(reader._mmap, self.mmap)

    def test_remarks(self, universe):
        title = universe.trajectory.remarks.decode('utf-8')
        assert_equal(self.title, title)

    def test_n_atoms(self, universe):
        assert_equal(self.n_atoms, len(universe.atoms))

    def test_numres(self, universe):
        assert_equal(self.n_res, len(universe.residues))

    def test_time(self, universe):
        assert_equal(self.time, universe.trajectory.time)

    def test_frame(self, universe):
        assert_equal(universe.trajectory.frame, 0)

    # Test slices as per base
    def test_full_slice(self, universe):
        trj_iter = universe.trajectory[:]
        frames = [ts.frame for ts in trj_iter]
        assert_equal(frames, np.arange(universe.trajectory.n_frames))

    def test_last_slice(self, universe):
        trj_iter = universe.trajectory[-1:]
        frames = [ts.frame for ts in trj_iter]
        assert_equal(frames, np.arange(universe.trajectory.n_frames))

    def test_num_frames(self, universe):
        assert_equal(universe.trajectory.n_frames, 1)

    def test_dt(self, universe):
        with pytest.warns(UserWarning) as record:
            assert_equal(universe.trajectory.dt, 1.0)

        assert len(record) == 1
        wmsg = "Reader has no dt information, set to 1.0 ps"
        assert str(record[0].message.args[0]) == wmsg

    def test_units(self, universe):
        units = {'time': 'ps', 'length': 'Angstrom',
                 'velocity': 'Angstrom/ps', 'force': 'kcal/(mol*Angstrom)'}
        assert_equal(universe.trajectory.units, units)

    def test_has_velocities(self, universe):
        assert_equal(self.has_velocities, universe.trajectory.has_velocities)

    def test_has_forces(self, universe):
        assert_equal(self.has_forces, universe.trajectory.has_forces)


class Read_ACECFV(_NCRestartTest):
    """ACE in water with Coordinates, Forces and Velocities"""

    # Placeholder values
    # topology = ACE
    # filename = ACE-PARM7
    title = 'Cpptraj Generated Restart'
    n_atoms = 1398
    n_res = 465
    time = 1.0
    has_velocities = True
    has_forces = True

    def test_positions(self, universe, selstr, expected):
        ag = universe.select_atoms(selstr)
        assert_almost_equal(expected, ag.positions, self.prec)

    def test_velocities(self, universe, selstr, expected):
        ag = universe.select_atoms(selstr)
        assert_almost_equal(expected, ag.velocities, self.prec)

    def test_forces(self, universe, selstr, expected):
        ag = universe.select_atoms(selstr)
        # We convert from kj back to the original kcal value
        assert_almost_equal(expected, ag.forces * 4.184, self.prec)


class ExceptionsNCRSTGenerator(object):
    """A basic modular ncrst writer based on :class:`NCDFWriter`"""

    def create_ncrst(self):
        # Create under context manager
        with netcdf.netcdf_file(self.filename, mode='w',
                                version=self.version_byte) as ncrst:
            # Top level attributes
            if self.Conventions:
                setattr(ncrst, 'Conventions', self.Conventions)
            if self.ConventionVersion:
                setattr(ncrst, 'ConventionVersion', self.ConventionVersion)
            if self.program:
                setattr(ncrst, 'program', self.program)
            if self.programVersion:
                setattr(ncrst, 'programVersion', self.programVersion)

            # Dimensions
            if self.n_atoms:
                ncrst.createDimension('atom', self.n_atoms)
            if self.spatial:
                ncrst.createDimension('spatial', self.spatial)

            # Variables
            if self.time:
                time = ncrst.createVariable('time', 'f8')
                time = 1.0
                setattr(time, 'units', self.time)
            # Spatial or atom dependent variables
            if (self.spatial) and (self.n_atoms):
                if self.coords:
                    coords = ncrst.createVariable('coordinates', 'f8',
                                                  ('atom', 'spatial'))
                    setattr(coords, 'units', self.coords)
                    coords[:] = np.asarray(range(self.spatial),
                                           dtype=np.float32)
                spatial = ncrst.createVariable('spatial', 'c', ('spatial',))
                spatial[:] = np.asarray(list('xyz'))
                velocs = ncrst.createVariable('velocities', 'f8',
                                              ('atom', 'spatial'))
                setattr(velocs, 'units', 'angstrom/picosecond')
                velocs[:] = np.asarray(range(self.spatial), dtype=np.float32)
                forces = ncrst.createVariable('forces', 'f8',
                                              ('atom', 'spatial'))
                setattr(forces, 'units', 'kilocalorie/mole/angstrom')
                forces[:] = np.asarray(range(self.spatial), dtype=np.float32)

            # self.scale_factor overrides which variable gets a scale_factor
            if self.scale_factor:
                setattr(ncrst.variables[self.scale_factor],
                        'scale_factor', 2.0)

    # Default definitions
    filename = 'test.ncrst'
    version_byte = 2
    Conventions = 'AMBERRESTART'
    ConventionVersion = '1.0'
    program = 'mda test_writer'
    programVersion = 'V42'
    n_atoms = 1
    spatial = 3
    coords = 'angstrom'
    time = 'picosecond'
    scale_factor = None


def TestNCRSTReaderExceptions(ExceptionsNCRSTGenerator):

   # TO DO: merge these so you get TypeErrors, KeyErrors, etc...

    def test_VersionbyteError(self, tmpdir):
        self.version_byte = 1
        with tmpdir.as_cwd():
            self.create_ncrst()
            with pytest.raises(TypeError):
                NCRSTReader(self.filename)

    def test_WrongConventionsError(self, tmpdir):
        self.Conventions = 'Foo'
        with tmpdir.as_cwd():
            self.create_ncrst()
            with pytest.raises(TypeError):
                NCRSTReader(self.filename)

    def test_NoConventionsError(self, tmpdir):
        self.Conventions = None
        with tmpdir.as_cwd():
            self.create_ncrst()
            with pytest.raises(KeyError):
                NCRSTReader(self.filename)

    def test_WrongConventionVersionError(self, tmpdir):
        self.ConventionVersion = '2.0'
        with tmpdir.as_cwd():
            self.create_ncrst()
            with pytest.raises(TypeError):
                NCRSTReader(self.filename)

    def test_NoConventionVersionError(self, tmpdir):
        self.ConventionVersion = None
        with tmpdir.as_cwd():
            self.create_ncrst()
            with pytest.raises(KeyError):
                NCRSTReader(self.filename)

    def test_WrongSpatialError(self, tmpdir):
        self.spatial = 2
        with tmpdir.as_cwd():
            self.create_ncrst()
            with pytest.raises(TypeError):
                NCRSTReader(self.filename)

    def test_NoSpatialError(self, tmpdir):
        self.spatial = None
        with tmpdir.as_cwd():
            self.create_ncrst()
            with pytest.raises(KeyError):
                NCRSTReader(self.filename)

    def test_NoAtomsError(self, tmpdir):
        self.n_atoms = None
        with tmpdir.as_cwd():
            self.create_ncrst()
            with pytest.raises(KeyError):
                NCRSTReader(self.filename)

    def test_WrongTimeError(self, tmpdir):
        self.time = 'femtosecond'
        with tmpdir.as_cwd():
            self.create_ncrst()
            with pytest.raises(NotImplementedError):
                NCRSTReader(self.filename)

    def test_WrongCoordsError(self, tmpdir):
        self.coords = 'nanometer'
        with tmpdir.as_cwd():
            self.create_ncrst()
            with pytest.raises(NotImplementedError):
                NCRSTReader(self.filename)

    def test_NoCoordsError(self, tmpdir):
        self.coords = None
        with tmpdir.as_cwd():
            self.create_ncrst()
            with pytest.raises(KeyError):
                NCRSTReader(self.filename)

    @pytest.mark.parametrize('variable', [
        'time', 'coordinates', 'spatial',
        'velocities', 'forces'
    ])
    def test_scale_factor(self, tmpdir, variable):
        self.scale_factor = variable
        with tmpdir.as_cwd():
            self.create_ncrst()
            with pytest.raises(NotImplementedError):
                NCRSTReader(self.filename)

# Errors:
# scale_factor
# Warnings:
# Wrong ConventionVersion
# Missing program value
# Missing programVersion value
# No time
