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
import MDAnalysis as mda
import numpy as np
import sys

from scipy.io import netcdf_file

import pytest
from numpy.testing import (
    assert_equal,
    assert_almost_equal
)
from MDAnalysis.coordinates.TRJ import NCDFReader, NCDFWriter

from MDAnalysisTests.datafiles import (PFncdf_Top, PFncdf_Trj,
                                       GRO, TRR, XYZ_mini,
                                       PRM_NCBOX, TRJ_NCBOX, DLP_CONFIG,
                                       CPPTRAJ_TRAJ_TOP, CPPTRAJ_TRAJ)
from MDAnalysisTests.coordinates.test_trj import _TRJReaderTest
from MDAnalysisTests.coordinates.reference import (RefVGV, RefTZ2)
from MDAnalysisTests import make_Universe
from MDAnalysisTests.util import block_import


class _NCDFReaderTest(_TRJReaderTest):
    prec = 3

    @pytest.fixture()
    def universe(self):
        return mda.Universe(self.topology, self.filename)

    def test_slice_iteration(self, universe):
        frames = [ts.frame for ts in universe.trajectory[4:-2:4]]
        assert_equal(frames,
                     np.arange(universe.trajectory.n_frames)[4:-2:4],
                     err_msg="slicing did not produce the expected frames")

    def test_metadata(self, universe):
        data = universe.trajectory.trjfile
        assert_equal(data.Conventions.decode('utf-8'), 'AMBER')
        assert_equal(data.ConventionVersion.decode('utf-8'), '1.0')

    def test_dt(self, universe):
        ref = 0.0
        assert_almost_equal(ref, universe.trajectory.dt, self.prec)
        assert_almost_equal(ref, universe.trajectory.ts.dt, self.prec)

    def test_get_writer(self, universe):
        with universe.trajectory.Writer('out.ncdf') as w:
            assert w.n_atoms == len(universe.atoms)
            assert w.remarks.startswith('AMBER NetCDF format')

    def test_get_writer_custom_n_atoms(self, universe):
        with universe.trajectory.Writer('out.ncdf', n_atoms=42,
                                        remarks='Hi!') as w:
            assert w.n_atoms == 42
            assert w.remarks == 'Hi!'

    def test_wrong_natoms(self):
        with pytest.raises(ValueError):
            mda.coordinates.TRJ.NCDFReader(self.filename, n_atoms=2)

    def test_read_on_closed(self, universe):
        universe.trajectory.close()

        with pytest.raises(IOError):
            universe.trajectory.__getitem__(2)

    def test_mmap_kwarg(self, universe):
        # default is None
        assert universe.trajectory._mmap == None


# Ugly way to create the tests for mmap

class _NCDFReaderTest_mmap_None(_NCDFReaderTest):
    @pytest.fixture()
    def universe(self):
        return mda.Universe(self.topology, self.filename, mmap=None)


class _NCDFReaderTest_mmap_True(_NCDFReaderTest):
    @pytest.fixture()
    def universe(self):
        return mda.Universe(self.topology, self.filename, mmap=True)

    def test_mmap_kwarg(self, universe):
        # default is None
        assert universe.trajectory._mmap == True


class _NCDFReaderTest_mmap_False(_NCDFReaderTest):
    @pytest.fixture()
    def universe(self):
        return mda.Universe(self.topology, self.filename, mmap=False)

    def test_mmap_kwarg(self, universe):
        assert universe.trajectory._mmap == False


class TestNCDFReader(_NCDFReaderTest, RefVGV):
    pass


class TestNCDFReader_mmap_None(_NCDFReaderTest_mmap_None, RefVGV):
    pass


class TestNCDFReader_mmap_True(_NCDFReaderTest_mmap_True, RefVGV):
    pass


class TestNCDFReader_mmap_False(_NCDFReaderTest_mmap_False, RefVGV):
    pass


class TestNCDFReaderTZ2(_NCDFReaderTest, RefTZ2):
    pass


class TestNCDFReader2(object):
    """NCDF Trajectory with positions and forces.

    Contributed by Albert Solernou
    """
    prec = 3

    @pytest.fixture(scope='class')
    def u(self):
        return mda.Universe(PFncdf_Top, PFncdf_Trj)

    def test_positions_1(self, u):
        """Check positions on first frame"""
        u.trajectory[0]
        ref_1 = np.array([[-0.11980818, 18.70524979, 11.6477766
                           ], [-0.44717646, 18.61727142, 12.59919548],
                          [-0.60952115, 19.47885513, 11.22137547]],
                         dtype=np.float32)
        assert_almost_equal(ref_1, u.atoms.positions[:3], self.prec)

    def test_positions_2(self, u):
        """Check positions on second frame"""
        u.trajectory[1]
        ref_2 = np.array([[-0.13042036, 18.6671524, 11.69647026
                           ], [-0.46643803, 18.60186768, 12.646698],
                          [-0.46567637, 19.49173927, 11.21922874]],
                         dtype=np.float32)
        assert_almost_equal(ref_2, u.atoms.positions[:3], self.prec)

    def test_forces_1(self, u):
        """Check forces on first frame"""
        u.trajectory[0]
        ref_1 = np.array([[49.23017883, -97.05565643, -86.09863281
                           ], [2.97547197, 29.84169388, 11.12069607],
                          [-15.93093777, 14.43616867, 30.25889015]],
                         dtype=np.float32)
        assert_almost_equal(ref_1, u.atoms.forces[:3], self.prec)

    def test_forces_2(self, u):
        """Check forces on second frame"""
        u.trajectory[1]
        ref_2 = np.array([[116.39096832, -145.44448853, -151.3155365
                           ], [-18.90058327, 27.20145798, 1.95245135],
                          [-31.08556366, 14.95863628, 41.10367966]],
                         dtype=np.float32)
        assert_almost_equal(ref_2, u.atoms.forces[:3], self.prec)

    def test_time_1(self, u):
        """Check time on first frame"""
        ref = 35.02
        assert_almost_equal(ref, u.trajectory[0].time, self.prec)

    def test_time_2(self, u):
        """Check time on second frame"""
        ref = 35.04
        assert_almost_equal(ref, u.trajectory[1].time, self.prec)

    def test_dt(self, u):
        ref = 0.02
        assert_almost_equal(ref, u.trajectory.dt, self.prec)
        assert_almost_equal(ref, u.trajectory.ts.dt, self.prec)

    def test_box(self, u):
        for ts in u.trajectory:
            assert ts.dimensions is None


class TestNCDFReader3(object):
    """NCDF trajectory with box, positions, forces and velocities

    Added to address Issue #2323
    """
    prec = 3

    # Expected coordinates as stored in Angstrom units
    coord_refs = np.array([
        [[15.249873, 12.578178, 15.191731],
         [14.925511, 13.58888, 14.944009],
         [15.285703, 14.3409605, 15.645962]],
        [[14.799454, 15.214347, 14.714555],
         [15.001984, 15.870884, 13.868363],
         [16.03358, 16.183628, 14.02995]]
    ], dtype=np.float32)

    # Expected forces as stored in kcal/(mol*Angstrom)
    frc_refs = np.array([
        [[8.583388, 1.8023694, -15.0033455],
         [-21.594835, 39.09166, 6.567963],
         [4.363016, -12.135163, 4.4775457]],
        [[-10.106646, -7.870829, -10.385734],
         [7.23599, -12.366022, -9.106191],
         [-4.637955, 11.597565, -6.463743]]
    ], dtype=np.float32)

    # Expected velocities as stored in Angstrom per AKMA time unit
    # These are usually associated with a scale_factor of 20.455
    vel_refs = np.array([
        [[-0.5301689, -0.16311595, -0.31390688],
         [0.00188578, 0.02513031, -0.2687525],
         [0.84072256, 0.09402391, -0.7457009]],
        [[-1.7773226, 1.2307, 0.50276583],
         [-0.13532305, 0.1355039, -0.05567304],
         [-0.6182481, 1.6396415, 0.46686798]]
    ], dtype=np.float32)

    # Expected box values as stored in cell_lengths ([:3]) of Angstrom
    # and cell_angles ([:3]) of degree
    box_refs = np.array([
        [28.81876287, 28.27875261, 27.72616397, 90., 90., 90.],
        [27.06266081, 26.55555665, 26.03664058, 90., 90., 90.]
    ], dtype=np.float32)

    @pytest.fixture(scope='class')
    def universe(self):
        return mda.Universe(PRM_NCBOX, TRJ_NCBOX)

    @pytest.mark.parametrize('index,expected', ((0, 0), (8, 1)))
    def test_positions(self, universe, index, expected):
        universe.trajectory[index]
        assert_almost_equal(self.coord_refs[expected],
                            universe.atoms.positions[:3], self.prec)

    @pytest.mark.parametrize('index,expected', ((0, 0), (8, 1)))
    def test_forces(self, universe, index, expected):
        """Here we multiply the forces by 4.184 to convert from
        kcal to kj in order to verify that MDA has correctly read
        and converted the units from those stored in the NetCDF file.
        """
        universe.trajectory[index]
        assert_almost_equal(self.frc_refs[expected] * 4.184,
                            universe.atoms.forces[:3], self.prec)

    @pytest.mark.parametrize('index,expected', ((0, 0), (8, 1)))
    def test_velocities(self, universe, index, expected):
        """Here we multiply the velocities by 20.455 to match the value of
        `scale_factor` which has been declared in the NetCDF file, which
        should change the values from Angstrom/AKMA time unit to Angstrom/ps.
        """
        universe.trajectory[index]
        assert_almost_equal(self.vel_refs[expected] * 20.455,
                            universe.atoms.velocities[:3], self.prec)

    @pytest.mark.parametrize('index,expected', ((0, 1.0), (8, 9.0)))
    def test_time(self, universe, index, expected):
        assert_almost_equal(expected, universe.trajectory[index].time,
                            self.prec)

    def test_nframes(self, universe):
        assert_equal(10, universe.trajectory.n_frames)

    def test_dt(self, universe):
        ref = 1.0
        assert_almost_equal(ref, universe.trajectory.dt, self.prec)
        assert_almost_equal(ref, universe.trajectory.ts.dt, self.prec)

    @pytest.mark.parametrize('index,expected', ((0, 0), (8, 1)))
    def test_box(self, universe, index, expected):
        universe.trajectory[index]
        assert_almost_equal(self.box_refs[expected], universe.dimensions)


class TestNCDFReader4(object):
    """NCDF Trajectory exported by cpptaj, without `time` variable."""
    prec = 3

    @pytest.fixture(scope='class')
    def u(self):
        return mda.Universe(CPPTRAJ_TRAJ_TOP,
                            [CPPTRAJ_TRAJ, CPPTRAJ_TRAJ])

    def test_chain_times(self, u):
        """Check times entries for a chain of trajectories without
        a defined time variable"""
        ref_times = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
        time_list = [ts.time for ts in u.trajectory]
        assert ref_times == time_list

    def test_dt(self, u):
        ref = 1.0
        assert u.trajectory.dt == pytest.approx(ref)
        assert u.trajectory.ts.dt == pytest.approx(ref)

    def test_warn_user_no_time_information(self, u):
        wmsg = ("NCDF trajectory does not contain `time` information;"
                " `time` will be set as an increasing index")  
        with pytest.warns(UserWarning, match=wmsg):
            u2 = mda.Universe(CPPTRAJ_TRAJ_TOP, CPPTRAJ_TRAJ)


class _NCDFGenerator(object):
    """A class for generating abitrary ncdf files and exhaustively test
    edge cases which might not be found in the wild"""

    def create_ncdf(self, params):
        """A basic modular ncdf writer based on :class:`NCDFWriter`"""
        # Create under context manager
        with netcdf_file(params['filename'], mode='w',
                         version=params['version_byte']) as ncdf:
            # Top level attributes
            if params['Conventions']:
                setattr(ncdf, 'Conventions', params['Conventions'])
            if params['ConventionVersion']:
                setattr(ncdf, 'ConventionVersion',
                        params['ConventionVersion'])
            if params['program']:
                setattr(ncdf, 'program', params['program'])
            if params['programVersion']:
                setattr(ncdf, 'programVersion', params['programVersion'])

            # Dimensions
            if params['frame']:
                ncdf.createDimension('frame', None)
            if params['n_atoms']:
                ncdf.createDimension('atom', params['n_atoms'])
            if params['spatial']:
                ncdf.createDimension('spatial', params['spatial'])
            if params['time']:
                ncdf.createDimension('time', 1)
            ncdf.createDimension('label', 5)
            ncdf.createDimension('cell_spatial', 3)
            ncdf.createDimension('cell_angular', 3)

            # Variables
            if params['time']:
                time = ncdf.createVariable('time', 'd', ('time',))
                setattr(time, 'units', params['time'])
                time[:] = 1.0
            cell_spatial = ncdf.createVariable('cell_spatial', 'c',
                                               ('cell_spatial', ))
            cell_spatial[:] = np.asarray(list('abc'))
            cell_angular = ncdf.createVariable('cell_angular', 'c',
                                               ('cell_angular', 'label'))
            cell_angular[:] = np.asarray([list('alpha'), list('beta '),
                                          list('gamma')])

            # Spatial or atom dependent variables
            if (params['spatial']) and (params['n_atoms']):
                spatial = ncdf.createVariable('spatial', 'c', ('spatial',))
                spatial[:] = np.asarray(list('xyz')[:params['spatial']])
                if params['frame']:
                    if params['coordinates']:
                        coords = ncdf.createVariable('coordinates', 'f4',
                                                     ('frame', 'atom',
                                                      'spatial'))
                    velocs = ncdf.createVariable('velocities', 'f4',
                                                 ('frame', 'atom', 'spatial'))
                    forces = ncdf.createVariable('forces', 'f4',
                                                 ('frame', 'atom', 'spatial'))
                    cell_lengths = ncdf.createVariable('cell_lengths', 'f8',
                                                       ('frame',
                                                        'cell_spatial'))
                    cell_angles = ncdf.createVariable('cell_angles', 'f8',
                                                      ('frame',
                                                       'cell_angular'))
                else:
                    if params['coordinates']:
                        coords = ncdf.createVariable('coordinates', 'f8',
                                                     ('atom', 'spatial'))
                    cell_lengths = ncdf.createVariable('cell_lengths', 'f8',
                                                       ('cell_spatial',))
                    cell_angles = ncdf.createVariable('cell_angles', 'f8',
                                                      ('cell_angular',))
                    velocs = ncdf.createVariable('velocities', 'f8',
                                                 ('atom', 'spatial'))
                    forces = ncdf.createVariable('forces', 'f8',
                                                 ('atom', 'spatial'))

                # Set units
                if params['coordinates']:
                    setattr(coords, 'units', params['coordinates'])
                setattr(velocs, 'units', params['velocities'])
                setattr(forces, 'units', params['forces'])
                setattr(cell_lengths, 'units', params['cell_lengths'])
                setattr(cell_angles, 'units', params['cell_angles'])

                # Assign value
                if params['frame']:
                    for index in range(params['frame']):
                        if params['coordinates']:
                            coords[index, :] = np.asarray(
                             range(params['spatial']), dtype=np.float32)
                        cell_lengths[index, :] = np.array([20., 20., 20.],
                                                          dtype=np.float32)
                        cell_angles[index, :] = np.array([90., 90., 90.],
                                                         dtype=np.float32)
                        velocs[index, :] = np.asarray(range(params['spatial']),
                                                      dtype=np.float32)
                        forces[index, :] = np.asarray(range(params['spatial']),
                                                      dtype=np.float32)
                else:
                    if params['coordinates']:
                        coords[:] = np.asarray(range(params['spatial']),
                                               dtype=np.float32)
                    cell_lengths[:] = np.array([20., 20., 20.],
                                               dtype=np.float32)
                    cell_angles[:] = np.array([90., 90., 90.],
                                              dtype=np.float32)
                    velocs[:] = np.asarray(range(params['spatial']),
                                           dtype=np.float32)
                    forces[:] = np.asarray(range(params['spatial']),
                                           dtype=np.float32)

            # self.scale_factor overrides which variable gets a scale_factor
            if params['scale_factor']:
                setattr(ncdf.variables[params['scale_factor']],
                        'scale_factor', params['scale_factor_value'])

    def gen_params(self, keypair=None, restart=False):
        """Generate writer parameters, keypair can be used to overwrite
        given dictonary entries (expects dictionary)
        """

        params = {
            'filename': 'test.nc',
            'version_byte': 2,
            'Conventions': 'AMBER',
            'ConventionVersion': '1.0',
            'program': 'mda test_writer',
            'programVersion': 'V42',
            'n_atoms': 1,
            'spatial': 3,
            'coordinates': 'angstrom',
            'velocities': 'angstrom/picosecond',
            'forces': 'kilocalorie/mole/angstrom',
            'cell_lengths': 'angstrom',
            'cell_angles': 'degree',
            'time': 'picosecond',
            'scale_factor': None,
            'scale_factor_value': 2.0,
            'frame': 2
        }

        if restart:
            params['filename'] = 'test.ncrst'
            params['frame'] = None

        if keypair:
            for entry in keypair:
                params[entry] = keypair[entry]

        return params


class TestScaleFactorImplementation(_NCDFGenerator):

    prec = 5

    def test_scale_factor_coordinates(self, tmpdir):
        mutation = {'scale_factor': 'coordinates'}
        params = self.gen_params(keypair=mutation, restart=False)
        expected = np.asarray(range(3), dtype=np.float32) * 2.0
        with tmpdir.as_cwd():
            self.create_ncdf(params)
            u = mda.Universe(params['filename'], to_guess=())
            for ts in u.trajectory:
                assert_almost_equal(ts.positions[0], expected, self.prec)

    def test_scale_factor_velocities(self, tmpdir):
        mutation = {'scale_factor': 'velocities', 'scale_factor_value': 3.0}
        params = self.gen_params(keypair=mutation, restart=False)
        expected = np.asarray(range(3), dtype=np.float32) * 3.0
        with tmpdir.as_cwd():
            self.create_ncdf(params)
            u = mda.Universe(params['filename'], to_guess=())
            for ts in u.trajectory:
                assert_almost_equal(ts.velocities[0], expected, self.prec)

    def test_scale_factor_forces(self, tmpdir):
        mutation = {'scale_factor': 'forces', 'scale_factor_value': 10.0}
        params = self.gen_params(keypair=mutation, restart=False)
        expected = np.asarray(range(3), dtype=np.float32) * 10.0 * 4.184
        with tmpdir.as_cwd():
            self.create_ncdf(params)
            u = mda.Universe(params['filename'], to_guess=())
            for ts in u.trajectory:
                assert_almost_equal(ts.forces[0], expected, self.prec)

    @pytest.mark.parametrize('mutation,expected', (
        ({'scale_factor': 'cell_lengths', 'scale_factor_value': 0.75},
         np.array([15., 15., 15., 90., 90., 90.])),
        ({'scale_factor': 'cell_angles', 'scale_factor_value': 0.5},
         np.array([20., 20., 20., 45., 45., 45.]))
    ))
    def test_scale_factor_box(self, tmpdir, mutation, expected):
        params = self.gen_params(keypair=mutation, restart=False)
        with tmpdir.as_cwd():
            self.create_ncdf(params)
            u = mda.Universe(params['filename'], to_guess=())
            for ts in u.trajectory:
                assert_almost_equal(ts.dimensions, expected, self.prec)

    def test_scale_factor_not_float(self, tmpdir):
        mutation = {'scale_factor': 'coordinates',
                    'scale_factor_value': 'parsnips'}
        params = self.gen_params(keypair=mutation, restart=False)
        with tmpdir.as_cwd():
            self.create_ncdf(params)
            errmsg = "b'parsnips' is not a float"
            with pytest.raises(TypeError, match=errmsg):
                u = mda.Universe(params['filename'])


class TestNCDFReaderExceptionsWarnings(_NCDFGenerator):

    @pytest.mark.parametrize('mutation', [
        {'Conventions': 'Foo'},
        {'version_byte': 1},
        {'spatial': 2}
    ])
    def test_type_errors(self, tmpdir, mutation):
        params = self.gen_params(keypair=mutation, restart=False)
        with tmpdir.as_cwd():
            self.create_ncdf(params)
            with pytest.raises(TypeError):
                NCDFReader(params['filename'])

    @pytest.mark.parametrize('mutation', [
        {'Conventions': None},
        {'ConventionVersion': None},
        {'spatial': None},
        {'n_atoms': None},
        {'frame': None}
    ])
    def test_value_errors(self, tmpdir, mutation):
        params = self.gen_params(keypair=mutation, restart=False)
        with tmpdir.as_cwd():
            self.create_ncdf(params)
            with pytest.raises(ValueError):
                NCDFReader(params['filename'])

    @pytest.mark.parametrize('mutation', [
        {'scale_factor': 'cell_spatial'},
        {'time': 'femtosecond'},
        {'coordinates': 'nanometer'},
        {'velocities': 'angstrom/akma'},
        {'forces': 'kilojoule/mole/angstrom'},
        {'cell_lengths': 'nanometer'},
        {'cell_angles': 'radians'}
    ])
    def test_notimplemented_errors(self, tmpdir, mutation):
        params = self.gen_params(keypair=mutation, restart=False)
        with tmpdir.as_cwd():
            self.create_ncdf(params)
            with pytest.raises(NotImplementedError):
                NCDFReader(params['filename'])

    @pytest.mark.parametrize('evaluate,expected', (
        ('yard', 'foot'),
        ('second', 'minute')
    ))
    def test_verify_units_errors(self, evaluate, expected):
        """Directly tests expected failures of _verify_units"""
        with pytest.raises(NotImplementedError):
            NCDFReader._verify_units(evaluate.encode('utf-8'), expected)

    def test_ioerror(self, tmpdir):
        params = self.gen_params(restart=False)
        with tmpdir.as_cwd():
            self.create_ncdf(params)
            with pytest.raises(IOError):
                u = mda.Universe(params['filename'], to_guess=())
                u.trajectory.close()
                u.trajectory[-1]

    def test_conventionversion_warn(self, tmpdir):
        mutation = {'ConventionVersion': '2.0'}
        params = self.gen_params(keypair=mutation, restart=False)
        with tmpdir.as_cwd():
            self.create_ncdf(params)
            with pytest.warns(UserWarning) as record:
                NCDFReader(params['filename'])

            assert len(record) == 1
            wmsg = ("NCDF trajectory format is 2.0 but the reader "
                    "implements format 1.0")
            assert str(record[0].message.args[0]) == wmsg

    @pytest.mark.parametrize('mutation', [
        {'program': None},
        {'programVersion': None}
    ])
    def test_program_warn(self, tmpdir, mutation):
        params = self.gen_params(keypair=mutation, restart=False)
        with tmpdir.as_cwd():
            self.create_ncdf(params)
            with pytest.warns(UserWarning) as record:
                NCDFReader(params['filename'])

            assert len(record) == 1
            wmsg = ("NCDF trajectory test.nc may not fully adhere to AMBER "
                    "standards as either the `program` or `programVersion` "
                    "attributes are missing")
            assert str(record[0].message.args[0]) == wmsg

    def test_no_dt_warning(self, tmpdir):
        """Issue 3166 - not being able to call dt should throw a warning.
        Also at the same time checks that single frame chain reading works"""
        u = mda.Universe(PFncdf_Top, PFncdf_Trj)
        with tmpdir.as_cwd():
            with NCDFWriter('single_frame.nc', u.atoms.n_atoms) as W:
                W.write(u)

            # Using the ChainReader implicitly calls dt() and thus _get_dt()
            wmsg = "Reader has no dt information, set to 1.0 ps"
            with pytest.warns(UserWarning, match=wmsg):
                u2 = mda.Universe(PFncdf_Top, [PFncdf_Trj, 'single_frame.nc'])


class _NCDFWriterTest(object):
    prec = 5

    @pytest.fixture()
    def universe(self):
        return mda.Universe(self.topology, self.filename)

    @pytest.fixture(params=['nc', 'ncdf'])
    def outfile_extensions(self, tmpdir, request):
        # Issue 3030, test all extensions of NCDFWriter
        ext = request.param
        return str(tmpdir) + f'ncdf-writer-1.{ext}'

    @pytest.fixture()
    def outfile(self, tmpdir):
        return str(tmpdir) + 'ncdf-writer-1.ncdf'

    @pytest.fixture()
    def outtop(self, tmpdir):
        return str(tmpdir) + 'ncdf-writer-top.pdb'

    def _test_write_trajectory(self, universe, outfile):
        # explicit import so that we can artifically remove netCDF4
        # before calling
        from MDAnalysis.coordinates import TRJ

        t = universe.trajectory
        with TRJ.NCDFWriter(outfile, t.n_atoms, dt=t.dt) as W:
            self._copy_traj(W, universe)
        self._check_new_traj(universe, outfile)
        # for issue #518 -- preserve float32 data in ncdf output
        # NOTE: This originally failed with the dtype('>f4') instead
        #       of dtype('<f4') == dtype('f') == np.float32, i.e. then
        #       endianness is different. The current hack-ish solution
        #       ignores endianness by comparing the name of the types,
        #       which should be "float32".
        #       See http://docs.scipy.org/doc/numpy-1.10.0/reference/arrays.dtypes.html
        #       and https://github.com/MDAnalysis/mdanalysis/pull/503
        dataset = netcdf_file(outfile, 'r')
        coords = dataset.variables['coordinates']
        time = dataset.variables['time']
        assert_equal(coords[:].dtype.name, np.dtype(np.float32).name,
                     err_msg='ncdf coord output not float32 '
                             'but {}'.format(coords[:].dtype))
        assert_equal(time[:].dtype.name, np.dtype(np.float32).name,
                     err_msg='ncdf time output not float32 '
                             'but {}'.format(time[:].dtype))

    def test_write_trajectory_netCDF4(self, universe, outfile):
        pytest.importorskip("netCDF4")
        return self._test_write_trajectory(universe, outfile)

    def test_write_trajectory_netcdf(self, universe, outfile):
        import MDAnalysis.coordinates.TRJ
        loaded_netCDF4 = sys.modules['MDAnalysis.coordinates.TRJ'].netCDF4
        try:
            # cannot use @block_import('netCDF4') because TRJ was already imported
            # during setup() and already sits in the global module list so we just
            # set it to None because that is what TRJ does if it cannot find netCDF4
            sys.modules['MDAnalysis.coordinates.TRJ'].netCDF4 = None
            assert MDAnalysis.coordinates.TRJ.netCDF4 is None  # should happen if netCDF4 not found
            return self._test_write_trajectory(universe, outfile)
        finally:
            sys.modules['MDAnalysis.coordinates.TRJ'].netCDF4 = loaded_netCDF4

    def test_OtherWriter(self, universe, outfile_extensions):
        t = universe.trajectory
        with t.OtherWriter(outfile_extensions) as W:
            self._copy_traj(W, universe)
        self._check_new_traj(universe, outfile_extensions)

    def _copy_traj(self, writer, universe):
        for ts in universe.trajectory:
            writer.write(universe)

    def _check_new_traj(self, universe, outfile):
        uw = mda.Universe(self.topology, outfile)

        # check that the trajectories are identical for each time step
        for orig_ts, written_ts in zip(universe.trajectory,
                                       uw.trajectory):
            assert_almost_equal(written_ts._pos, orig_ts._pos, self.prec,
                                      err_msg="coordinate mismatch between "
                                              "original and written trajectory at "
                                              "frame %d (orig) vs %d (written)" % (
                                                  orig_ts.frame,
                                                  written_ts.frame))
            # not a good test because in the example trajectory all times are 0
            assert_almost_equal(orig_ts.time, written_ts.time, self.prec,
                                err_msg="Time for step {0} are not the "
                                        "same.".format(orig_ts.frame))
            assert_almost_equal(written_ts.dimensions,
                                      orig_ts.dimensions,
                                      self.prec,
                                      err_msg="unitcells are not identical")
        # check that the NCDF data structures are the same
        nc_orig = universe.trajectory.trjfile
        nc_copy = uw.trajectory.trjfile

        # note that here 'dimensions' is a specific netcdf data structure and
        # not the unit cell dimensions in MDAnalysis
        for k, dim in nc_orig.dimensions.items():
            try:
                dim_new = nc_copy.dimensions[k]
            except KeyError:
                raise AssertionError("NCDFWriter did not write "
                                     "dimension '{0}'".format(k))
            else:
                assert_equal(dim, dim_new,
                             err_msg="Dimension '{0}' size mismatch".format(k))

        for k, v in nc_orig.variables.items():
            try:
                v_new = nc_copy.variables[k]
            except KeyError:
                raise AssertionError("NCDFWriter did not write "
                                     "variable '{0}'".format(k))
            else:
                try:
                    assert_almost_equal(v[:], v_new[:], self.prec,
                                              err_msg="Variable '{0}' not "
                                                      "written correctly".format(
                                                  k))
                except TypeError:
                    assert_equal(v[:], v_new[:],
                                       err_msg="Variable {0} not written "
                                               "correctly".format(k))

    def test_TRR2NCDF(self, outfile):
        trr = mda.Universe(GRO, TRR)
        with mda.Writer(outfile, trr.trajectory.n_atoms,
                        velocities=True, format="ncdf") as W:
            for ts in trr.trajectory:
                W.write(trr)

        uw = mda.Universe(GRO, outfile)

        for orig_ts, written_ts in zip(trr.trajectory,
                                       uw.trajectory):
            assert_almost_equal(written_ts._pos, orig_ts._pos, self.prec,
                                      err_msg="coordinate mismatch between "
                                              "original and written trajectory at "
                                              "frame {0} (orig) vs {1} (written)".format(
                                          orig_ts.frame, written_ts.frame))
            assert_almost_equal(written_ts._velocities,
                                      orig_ts._velocities, self.prec,
                                      err_msg="velocity mismatch between "
                                              "original and written trajectory at "
                                              "frame {0} (orig) vs {1} (written)".format(
                                          orig_ts.frame, written_ts.frame))
            assert_almost_equal(orig_ts.time, written_ts.time, self.prec,
                                err_msg="Time for step {0} are not the "
                                        "same.".format(orig_ts.frame))
            assert_almost_equal(written_ts.dimensions,
                                      orig_ts.dimensions,
                                      self.prec,
                                      err_msg="unitcells are not identical")
        del trr

    def test_write_AtomGroup(self, universe, outfile, outtop):
        """test to write NCDF from AtomGroup (Issue 116)"""
        p = universe.select_atoms("not resname WAT")
        p.write(outtop)
        with mda.Writer(outfile, n_atoms=p.n_atoms, format="ncdf") as W:
            for ts in universe.trajectory:
                W.write(p)

        uw = mda.Universe(outtop, outfile)
        pw = uw.atoms

        for orig_ts, written_ts in zip(universe.trajectory,
                                       uw.trajectory):
            assert_almost_equal(p.positions, pw.positions, self.prec,
                                      err_msg="coordinate mismatch between "
                                              "original and written trajectory at "
                                              "frame %d (orig) vs %d (written)" % (
                                                  orig_ts.frame,
                                                  written_ts.frame))
            assert_almost_equal(orig_ts.time, written_ts.time, self.prec,
                                err_msg="Time for step {0} are not the "
                                        "same.".format(orig_ts.frame))
            assert_almost_equal(written_ts.dimensions,
                                      orig_ts.dimensions,
                                      self.prec,
                                      err_msg="unitcells are not identical")


class TestNCDFWriter(_NCDFWriterTest, RefVGV):
    pass


class TestNCDFWriterTZ2(_NCDFWriterTest, RefTZ2):
    pass


class TestNCDFWriterVelsForces(object):
    """Test writing NCDF trajectories with a mixture of options"""
    prec = 3
    top = XYZ_mini
    n_atoms = 3

    @pytest.fixture()
    def u1(self):
        u = make_Universe(size=(self.n_atoms, 1, 1), trajectory=True,
                          velocities=True, forces=True)
        # Memory reader so changes should be in-place
        u.atoms.velocities += 100
        u.atoms.forces += 200
        return u

    @pytest.fixture()
    def u2(self):
        u = make_Universe(size=(self.n_atoms, 1, 1), trajectory=True,
                          velocities=True, forces=True)
        # Memory reader so changes should be in-place
        u.atoms.positions += 300
        u.atoms.velocities += 400
        u.atoms.forces += 500
        return u

    @pytest.mark.parametrize('pos, vel, force', (
            (True, False, False),
            (True, True, False),
            (True, False, True),
            (True, True, True),
    ))
    def test_write_u(self, pos, vel, force, tmpdir, u1, u2):
        """Write the two reference universes, then open them up and check values

        pos vel and force are bools which define whether these properties
        should be in universe
        """
        outfile = str(tmpdir) + 'ncdf-write-vels-force.ncdf'
        with NCDFWriter(outfile,
                        n_atoms=self.n_atoms,
                        velocities=vel,
                        forces=force) as w:
            w.write(u1)
            w.write(u2)

        # test that the two reference states differ
        for ts1, ts2 in zip(u1.trajectory, u2.trajectory):
            assert_almost_equal(ts1._pos + 300, ts2._pos)
            assert_almost_equal(ts1._velocities + 300, ts2._velocities)
            assert_almost_equal(ts1._forces + 300, ts2._forces)

        u = mda.Universe(self.top, outfile)
        # check the trajectory contents match reference universes
        for ts, ref_ts in zip(u.trajectory, [u1.trajectory.ts, u2.trajectory.ts]):
            if pos:
                assert_almost_equal(ts._pos, ref_ts._pos, self.prec)
            else:
                with pytest.raises(mda.NoDataError):
                    getattr(ts, 'positions')
            if vel:
                assert_almost_equal(ts._velocities, ref_ts._velocities,
                                    self.prec)
            else:
                with pytest.raises(mda.NoDataError):
                    getattr(ts, 'velocities')
            if force:
                assert_almost_equal(ts._forces, ref_ts._forces, self.prec)
            else:
                with pytest.raises(mda.NoDataError):
                    getattr(ts, 'forces')

        u.trajectory.close()


class TestNCDFWriterScaleFactors:
    """Class to check that the NCDF Writer properly handles the
    application of scale factors"""

    @pytest.fixture()
    def outfile(self, tmpdir):
        return str(tmpdir) + 'ncdf-write-scale.ncdf'

    @pytest.fixture()
    def outfile2(self, tmpdir):
        return str(tmpdir) + 'ncdf-write-scale2.ncdf'

    @pytest.fixture()
    def universe(self):
        return mda.Universe(PRM_NCBOX, TRJ_NCBOX)

    def get_scale_factors(self, ncdfile):
        """Get a dictionary of scale factors stored in netcdf file"""
        sfactors = {}
        # being overly cautious by setting mmap to False, probably would
        # be faster & ok to set it to True
        with netcdf_file(ncdfile, mmap=False) as f:
            for var in f.variables:
                if hasattr(f.variables[var], 'scale_factor'):
                    sfactors[var] = f.variables[var].scale_factor

        return sfactors

    def get_variable(self, ncdfile, variable, frame):
        """Return a variable array from netcdf file"""
        with netcdf_file(ncdfile, mmap=False) as f:
            return f.variables[variable][frame]

    def test_write_read_factors_default(self, outfile, universe):
        with universe.trajectory.Writer(outfile) as W:
            W.write(universe.atoms)

        # check scale_factors
        sfactors = self.get_scale_factors(outfile)
        assert len(sfactors) == 1
        assert sfactors['velocities'] == 20.455

    def test_write_bad_scale_factor(self, outfile, universe):
        errmsg = "scale_factor parsnips is not a float"
        with pytest.raises(TypeError, match=errmsg):
            NCDFWriter(outfile, n_atoms=len(universe.atoms),
                       scale_velocities="parsnips")

    @pytest.mark.parametrize('stime, slengths, sangles, scoords, svels, sfrcs', (
            (-2.0, -2.0, -2.0, -2.0, -2.0, -2.0),
            (1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            (2.0, 4.0, 8.0, 16.0, 32.0, 64.0)
    ))
    def test_write_read_write(self, outfile, outfile2, universe, stime,
                              slengths, sangles, scoords, svels, sfrcs):
        """Write out a file with assorted scale_factors, then
        read it back in, then write it out to make sure that the
        new assorted scale_factors have been retained by Write"""

        with NCDFWriter(outfile, n_atoms=len(universe.atoms), velocities=True,
                        forces=True, scale_time=stime,
                        scale_cell_lengths=slengths, scale_cell_angles=sangles,
                        scale_coordinates=scoords, scale_velocities=svels,
                        scale_forces=sfrcs) as W:
            for ts in universe.trajectory:
                W.write(universe.atoms)

        universe2 = mda.Universe(PRM_NCBOX, outfile)

        # Write back checking that the same scale_factors as used
        with universe2.trajectory.Writer(outfile2) as OWriter:
            OWriter.write(universe2.atoms)

        universe3 = mda.Universe(PRM_NCBOX, outfile2)

        # check stored scale_factors
        sfactors1 = self.get_scale_factors(outfile)
        sfactors2 = self.get_scale_factors(outfile2)

        assert sfactors1 == sfactors2
        assert len(sfactors1) == 6
        assert sfactors1['time'] == stime
        assert sfactors1['cell_lengths'] == slengths
        assert sfactors1['cell_angles'] == sangles
        assert sfactors1['coordinates'] == scoords
        assert sfactors1['velocities'] == svels
        assert sfactors1['forces'] == sfrcs

        # check that the stored values are indeed scaled
        assert_almost_equal(universe.trajectory.time / stime,
                            self.get_variable(outfile, 'time', 0), 4)
        assert_almost_equal(universe.dimensions[:3] / slengths,
                            self.get_variable(outfile, 'cell_lengths', 0), 4)
        assert_almost_equal(universe.dimensions[3:] / sangles,
                            self.get_variable(outfile, 'cell_angles', 0), 4)
        assert_almost_equal(universe.atoms.positions / scoords,
                            self.get_variable(outfile, 'coordinates', 0), 4)
        assert_almost_equal(universe.atoms.velocities / svels,
                            self.get_variable(outfile, 'velocities', 0), 4)
        # note: kJ/mol -> kcal/mol = 4.184 conversion
        assert_almost_equal(universe.atoms.forces / (sfrcs * 4.184),
                            self.get_variable(outfile, 'forces', 0), 4)

        # check that the individual components were saved/read properly
        for ts1, ts3 in zip(universe.trajectory, universe3.trajectory):
            assert_almost_equal(ts1.time, ts3.time)
            assert_almost_equal(ts1.dimensions, ts3.dimensions)
            assert_almost_equal(universe.atoms.positions,
                                universe3.atoms.positions, 4)
            assert_almost_equal(universe.atoms.velocities,
                                universe3.atoms.velocities, 4)
            assert_almost_equal(universe.atoms.forces,
                                universe3.atoms.forces, 4)


class TestScipyScaleFactors(TestNCDFWriterScaleFactors):
    """As above, but netCDF4 is disabled since scaleandmask is different
    between the two libraries"""
    @pytest.fixture(autouse=True)
    def block_netcdf4(self, monkeypatch):
        monkeypatch.setattr(sys.modules['MDAnalysis.coordinates.TRJ'], 'netCDF4', None)

    def test_ncdf4_not_present(self, outfile, universe):
        # whilst we're here, let's also test this warning
        with pytest.warns(UserWarning, match="Could not find netCDF4 module"):
            with universe.trajectory.Writer(outfile) as W:
                W.write(universe.atoms)


class TestNCDFWriterUnits(object):
    """Tests that the writer adheres to AMBER convention units"""
    @pytest.fixture()
    def outfile(self, tmpdir):
        return str(tmpdir) + 'ncdf-writer-1.ncdf'

    @pytest.mark.parametrize('var, expected', (
        ('coordinates', 'angstrom'),
        ('time', 'picosecond'),
        ('cell_lengths', 'angstrom'),
        ('cell_angles', 'degree'),
        ('velocities', 'angstrom/picosecond'),
        ('forces', 'kilocalorie/mole/angstrom')
    ))
    def test_writer_units(self, outfile, var, expected):
        trr = mda.Universe(DLP_CONFIG, format='CONFIG')

        with mda.Writer(outfile, trr.trajectory.n_atoms, velocities=True,
                        forces=True, format='ncdf') as W:
            for ts in trr.trajectory:
                W.write(trr)

        with netcdf_file(outfile, mode='r') as ncdf:
            unit = ncdf.variables[var].units.decode('utf-8')
            assert_equal(unit, expected)


class TestNCDFWriterErrorsWarnings(object):
    @pytest.fixture()
    def outfile(self, tmpdir):
        return str(tmpdir) + 'out.ncdf'

    def test_zero_atoms_VE(self, outfile):
        with pytest.raises(ValueError):
            NCDFWriter(outfile, 0)

    def test_wrong_n_atoms(self, outfile):
        with NCDFWriter(outfile, 100) as w:
            u = make_Universe(trajectory=True)
            with pytest.raises(IOError):
                w.write(u)
