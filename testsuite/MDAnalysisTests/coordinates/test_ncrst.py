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

from MDAnalysisTests.datafiles import (PRMNCRSTbox, NCRSTbox)

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


class TestReadACEBox(_NCRestartTest):
    """ACE in water with Coordinates, Forces and Velocities"""

    # Placeholder values
    topology = PRMNCRSTbox
    filename = NCRSTbox
    title = 'Cpptraj Generated Restart'
    n_atoms = 1398
    n_res = 465
    time = 1.0
    has_velocities = True
    has_forces = True
    mmap = None

    def test_positions_resid1(self, universe):
        ag = universe.select_atoms('resid 1')
        expected = np.array([[16.91830635, 11.8711319, 15.46306515],
                             [16.02449417, 12.34413528, 15.86983871],
                             [16.05173683, 12.2475071, 16.95520592],
                             [15.09445763, 12.00988102, 15.41004944],
                             [16.09916306, 13.7953701, 15.57103825],
                             [16.00987625, 14.24524307, 14.39752007]],
                            dtype=np.float32)
        assert_almost_equal(expected, ag.positions, self.prec)

    def test_positions_resid42(self, universe):
        ag = universe.select_atoms('resid 42')
        expected = np.array([[20.67843246, 21.18344688, 20.79789162],
                             [19.85808182, 21.5282135, 21.15058327],
                             [20.97104073, 20.54719162, 21.45041847]],
                            dtype=np.float32)
        assert_almost_equal(expected, ag.positions, self.prec)

#    def test_velocities(self, universe, selstr, expected):
#        ag = universe.select_atoms(selstr)
#        assert_almost_equal(expected, ag.velocities, self.prec)
#
#    def test_forces(self, universe, selstr, expected):
#        ag = universe.select_atoms(selstr)
#        # We convert from kj back to the original kcal value
#        assert_almost_equal(expected, ag.forces * 4.184, self.prec)


class _NCRSTGenerator(object):
    """A basic modular ncrst writer based on :class:`NCDFWriter`"""

    def create_ncrst(self, params):
        # Create under context manager
        with netcdf.netcdf_file(params['filename'], mode='w',
                                version=params['version_byte']) as ncrst:
            # Top level attributes
            if params['Conventions']:
                setattr(ncrst, 'Conventions', params['Conventions'])
            if params['ConventionVersion']:
                setattr(ncrst, 'ConventionVersion',
                        params['ConventionVersion'])
            if params['program']:
                setattr(ncrst, 'program', params['program'])
            if params['programVersion']:
                setattr(ncrst, 'programVersion', params['programVersion'])

            # Dimensions
            if params['n_atoms']:
                ncrst.createDimension('atom', params['n_atoms'])
            if params['spatial']:
                ncrst.createDimension('spatial', params['spatial'])
            if params['time']:
                ncrst.createDimension('time', 1)

            # Variables
            if params['time']:
                time = ncrst.createVariable('time', 'd', ('time',))
                setattr(time, 'units', params['time'])
                time[:] = 1.0
            # Spatial or atom dependent variables
            if (params['spatial']) and (params['n_atoms']):
                if params['coords']:
                    coords = ncrst.createVariable('coordinates', 'f8',
                                                  ('atom', 'spatial'))
                    setattr(coords, 'units', params['coords'])
                    coords[:] = np.asarray(range(params['spatial']),
                                           dtype=np.float32)
                spatial = ncrst.createVariable('spatial', 'c', ('spatial',))
                spatial[:] = np.asarray(list('xyz')[:params['spatial']])
                velocs = ncrst.createVariable('velocities', 'f8',
                                              ('atom', 'spatial'))
                setattr(velocs, 'units', 'angstrom/picosecond')
                velocs[:] = np.asarray(range(params['spatial']),
                                       dtype=np.float32)
                forces = ncrst.createVariable('forces', 'f8',
                                              ('atom', 'spatial'))
                setattr(forces, 'units', 'kilocalorie/mole/angstrom')
                forces[:] = np.asarray(range(params['spatial']),
                                       dtype=np.float32)

            # self.scale_factor overrides which variable gets a scale_factor
            if params['scale_factor']:
                setattr(ncrst.variables[params['scale_factor']],
                        'scale_factor', 2.0)

    def gen_params(self, key=None, value=None):
        """Generate writer parameters, key and value can be used to overwrite
        a given dictionary entry
        """

        params = {
            'filename': 'test.ncrst',
            'version_byte': 2,
            'Conventions': 'AMBERRESTART',
            'ConventionVersion': '1.0',
            'program': 'mda test_writer',
            'programVersion': 'V42',
            'n_atoms': 1,
            'spatial': 3,
            'coords': 'angstrom',
            'time': 'picosecond',
            'scale_factor': None
        }

        if key:
            params[key] = value

        return params


class TestNCRSTReaderExceptionsWarnings(_NCRSTGenerator):

    @pytest.mark.parametrize('key,value', (
        ('version_byte', 1),
        ('Conventions', 'Foo'),
        ('spatial', 2)
    ))
    def test_read_type_errors(self, tmpdir, key, value):
        params = self.gen_params(key=key, value=value)
        with tmpdir.as_cwd():
            self.create_ncrst(params)
            with pytest.raises(TypeError):
                NCRSTReader(params['filename'])

    @pytest.mark.parametrize('key,value', (
        ('Conventions', None),
        ('ConventionVersion', None)
    ))
    def test_attribute_errors(self, tmpdir, key, value):
        params = self.gen_params(key=key, value=value)
        with tmpdir.as_cwd():
            self.create_ncrst(params)
            with pytest.raises(AttributeError):
                NCRSTReader(params['filename'])

    @pytest.mark.parametrize('key,value', (
        ('spatial', None),
        ('n_atoms', None),
        ('coords', None)
    ))
    def test_key_errors(self, tmpdir, key, value):
        params = self.gen_params(key=key, value=value)
        with tmpdir.as_cwd():
            self.create_ncrst(params)
            with pytest.raises(KeyError):
                NCRSTReader(params['filename'])

    @pytest.mark.parametrize('key,value', (
        ('time', 'femtosecond'),
        ('coords', 'nanometer')
    ))
    def test_not_implemented_errors(self, tmpdir, key, value):
        params = self.gen_params(key=key, value=value)
        with tmpdir.as_cwd():
            self.create_ncrst(params)
            with pytest.raises(NotImplementedError):
                NCRSTReader(params['filename'])

    @pytest.mark.parametrize('value', [
        'time', 'coordinates', 'spatial',
        'velocities', 'forces'
    ])
    def test_scale_factor(self, tmpdir, value):
        params = self.gen_params(key='scale_factor', value=value)
        with tmpdir.as_cwd():
            self.create_ncrst(params)
            with pytest.raises(NotImplementedError):
                NCRSTReader(params['filename'])

    def test_conventionversion_warn(self, tmpdir):
        params = self.gen_params(key='ConventionVersion', value='2.0')
        with tmpdir.as_cwd():
            self.create_ncrst(params)
            with pytest.warns(UserWarning) as record:
                NCRSTReader(params['filename'])

            assert len(record) == 1
            wmsg = "NCRST format is 2.0 but the reader implements format 1.0"
            assert str(record[0].message.args[0]) == wmsg

    @pytest.mark.parametrize('key,value', (
        ('program', None),
        ('programVersion', None)
    ))
    def test_program_warn(self, tmpdir, key, value):
        params = self.gen_params(key=key, value=value)
        with tmpdir.as_cwd():
            self.create_ncrst(params)
            with pytest.warns(UserWarning) as record:
                NCRSTReader(params['filename'])

            assert len(record) == 1
            wmsg = ("This NCRST file may not fully adhere to AMBER "
                    "standards as either the `program` or `programVersion` "
                    "attributes are missing")
            assert str(record[0].message.args[0]) == wmsg

    def test_notime_warn(self, tmpdir):
        params = self.gen_params(key='time', value=None)
        with tmpdir.as_cwd():
            self.create_ncrst(params)
            with pytest.warns(UserWarning) as record:
                NCRSTReader(params['filename'])

            # Lack of time triggers two warnings
            assert len(record) == 2
            wmsg1 = ("NCRestart file {0} does not contain time information. "
                     "This should be expected if the file was not created "
                     "from an MD trajectory (e.g. a minimization)".format(
                      params['filename']))
            assert str(record[0].message.args[0]) == wmsg1
            wmsg2 = ("Reader has no dt information, set to 1.0 ps")
            assert str(record[1].message.args[0]) == wmsg2
