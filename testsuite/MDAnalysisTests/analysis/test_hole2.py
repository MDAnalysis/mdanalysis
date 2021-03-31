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
import pytest
import glob
import os
import sys
import textwrap

import numpy as np
import matplotlib
import mpl_toolkits.mplot3d
import errno
from numpy.testing import (
    assert_almost_equal,
    assert_equal,
)

import MDAnalysis as mda
from MDAnalysis.analysis import hole2
from MDAnalysis.analysis.hole2.utils import check_and_fix_long_filename
from MDAnalysis.exceptions import ApplicationError
from MDAnalysisTests.datafiles import PDB_HOLE, MULTIPDB_HOLE, DCD
from MDAnalysisTests import executable_not_found


def rlimits_missing():
    # return True if resources module not accesible (ie setting of rlimits)
    try:
        # on Unix we can manipulate our limits: http://docs.python.org/2/library/resource.html
        import resource

        soft_max_open_files, hard_max_open_files = resource.getrlimit(
            resource.RLIMIT_NOFILE)
    except ImportError:
        return True
    return False


class TestCheckAndFixLongFilename(object):

    max_length = 70
    filename = 'a'*(max_length-4) + '.pdb'

    def test_short(self):
        fixed = check_and_fix_long_filename(self.filename)
        assert self.filename == fixed

    def test_relative(self):
        abspath = os.path.abspath(self.filename)
        if len(abspath) > self.max_length:
            fixed = check_and_fix_long_filename(abspath)
            assert fixed == self.filename

    @pytest.mark.skipif(os.name == 'nt' and sys.maxsize <= 2**32,
                        reason="FileNotFoundError on Win 32-bit")
    def test_symlink_dir(self, tmpdir):
        dirname = 'really_'*20 +'long_name'
        short_name = self.filename[-20:]
        path = os.path.join(dirname, short_name)
        with tmpdir.as_cwd():
            os.makedirs(dirname)
            u = mda.Universe(PDB_HOLE)
            u.atoms.write(path)

            fixed = check_and_fix_long_filename(path)
            assert os.path.islink(fixed)
            assert fixed.endswith(short_name)

    @pytest.mark.skipif(os.name == 'nt' and sys.maxsize <= 2**32,
                        reason="OSError: symbolic link privilege not held")
    def test_symlink_file(self, tmpdir):
        long_name = 'a'*10 + self.filename

        with tmpdir.as_cwd():
            fixed = check_and_fix_long_filename(long_name)
            assert os.path.islink(fixed)
            assert not fixed.endswith(long_name)


@pytest.mark.skipif(executable_not_found("hole"),
                    reason="Test skipped because HOLE not found")
class TestHole(object):
    filename = PDB_HOLE
    random_seed = 31415
    profile_length = 425
    rxn_coord_mean = -1.41225
    radius_min = 1.19707

    def test_correct_input(self, tmpdir):
        with tmpdir.as_cwd():
            hole2.hole(self.filename, random_seed=self.random_seed,
                       infile='hole.inp')

        infile = str(tmpdir.join('hole.inp'))
        with open(infile, 'r') as f:
            contents = f.read()

        hole_input = textwrap.dedent("""
            RADIUS simple2.rad
            SPHPDB hole.sph
            SAMPLE 0.200000
            ENDRAD 22.000000
            IGNORE SOL WAT TIP HOH K   NA  CL 
            SHORTO 0
            RASEED 31415
            """)

        # don't check filenames
        assert contents.endswith(hole_input)

    def test_input_options(self, tmpdir):
        u = mda.Universe(PDB_HOLE)
        cog = u.select_atoms('protein').center_of_geometry()

        with tmpdir.as_cwd():
            hole2.hole(self.filename, random_seed=self.random_seed,
                       infile='hole.inp', cpoint=cog,
                       ignore_residues=[])

        infile = str(tmpdir.join('hole.inp'))
        with open(infile, 'r') as f:
            contents = f.read()

        hole_input = textwrap.dedent("""
            RADIUS simple2.rad
            SPHPDB hole.sph
            SAMPLE 0.200000
            ENDRAD 22.000000
            IGNORE 
            SHORTO 0
            RASEED 31415
            CPOINT -0.0180961507 -0.0122730583 4.1497999943
            """)

        # don't check filenames
        assert contents.endswith(hole_input)

    def test_correct_profile_values(self, tmpdir):
        with tmpdir.as_cwd():
            profiles = hole2.hole(self.filename, random_seed=self.random_seed)

        values = list(profiles.values())
        assert_equal(len(values), 1,
                     err_msg='Should only have 1 HOLE profile')
        profile = values[0]
        assert_equal(len(profile), self.profile_length,
                     err_msg='Wrong number of points in HOLE profile')
        assert_almost_equal(profile.rxn_coord.mean(),
                            self.rxn_coord_mean,
                            err_msg='Wrong mean HOLE rxn_coord')
        assert_almost_equal(profile.radius.min(),
                            self.radius_min,
                            err_msg='Wrong minimum HOLE radius')

    # HOLE cannot read NAMD CHARMM files as is written by MDAnalysis
    # fails with Linux HOLE?
    # def test_dcd(self, tmpdir):
    #     with tmpdir.as_cwd():
    #         u = mda.Universe(PSF, DCD)
    #         n_frames = len(u.trajectory)
    #         keep_frames = 12
    #         step = 3
    #         skip = n_frames-keep_frames
    #         filename = 'temp.pdb'
    #         u.atoms.write(filename)

    #         profiles = hole2.hole(filename, random_seed=self.random_seed, dcd=DCD,
    #                               dcd_iniskip=skip,
    #                               dcd_step=step, infile='hole.inp')
    #         with open('hole.inp', 'r') as f:
    #             contents = f.read()
    #         assert contents.endswith('CHARMS {} {}\n'.format(skip, step-1))

    #     assert_equal(len(profiles), int(keep_frames/step))

    def test_application_error(self, tmpdir):
        with tmpdir.as_cwd():
            with pytest.raises(ApplicationError):
                hole2.hole(self.filename, dcd=DCD)

    def test_output_level(self, tmpdir):
        with tmpdir.as_cwd():
            with pytest.warns(UserWarning, match="needs to be < 3"):
                profiles = hole2.hole(self.filename,
                                      random_seed=self.random_seed,
                                      output_level=100)
            assert len(profiles) == 0

    def test_keep_files(self, tmpdir):
        with tmpdir.as_cwd():
            hole2.hole(self.filename, random_seed=self.random_seed,
                       infile='hole.inp',
                       keep_files=False)
            sphpdbs = tmpdir.join('*.sph')
            assert len(glob.glob(str(sphpdbs))) == 0
            outfiles = tmpdir.join('*.out')
            assert len(glob.glob(str(outfiles))) == 0
            vdwradii = tmpdir.join('simple2.rad')
            assert len(glob.glob(str(vdwradii))) == 0
            pdbfiles = tmpdir.join('*.pdb')
            assert len(glob.glob(str(pdbfiles))) == 0
            oldfiles = tmpdir.join('*.old')
            assert len(glob.glob(str(oldfiles))) == 0


@pytest.mark.skipif(executable_not_found("hole"),
                    reason="Test skipped because HOLE not found")
class BaseTestHole(object):
    filename = MULTIPDB_HOLE
    start = 5
    stop = 7
    random_seed = 31415

    # HOLE is so slow so we only run it once and keep it in
    # the class; note that you may not change universe.trajectory
    # (eg iteration) because this is not safe in parallel
    @pytest.fixture()
    def universe(self):
        return mda.Universe(self.filename)

    @pytest.fixture()
    def hole(self, universe, tmpdir):
        with tmpdir.as_cwd():
            h = hole2.HoleAnalysis(universe)
            h.run(start=self.start, stop=self.stop,
                  random_seed=self.random_seed)
        return h

    @pytest.fixture()
    def frames(self, universe):
        return [ts.frame for ts in universe.trajectory[self.start:self.stop]]

    @pytest.fixture()
    def profiles(self, hole, frames):
        return [hole.profiles[f] for f in frames]


class TestHoleAnalysis(BaseTestHole):

    def test_correct_profile_values(self, hole, frames):
        assert_equal(sorted(hole.profiles.keys()), frames,
                     err_msg="hole.profiles.keys() should contain the frame numbers")
        assert_equal(list(hole.frames), frames,
                     err_msg="hole.frames should contain the frame numbers")
        data = np.transpose([(len(p), p.rxn_coord.mean(), p.radius.min())
                             for p in hole.profiles.values()])
        assert_equal(data[0], [401, 399], err_msg="incorrect profile lengths")
        assert_almost_equal(data[1], [1.98767,  0.0878],
                            err_msg="wrong mean HOLE rxn_coord")
        assert_almost_equal(data[2], [1.19819, 1.29628],
                            err_msg="wrong minimum radius")

    def test_min_radius(self, hole):
        values = np.array([[5., 1.19819],
                           [6., 1.29628]])
        assert_almost_equal(hole.min_radius(), values,
                            err_msg="min_radius() array not correct")

    def test_context_manager(self, universe, tmpdir):
        with tmpdir.as_cwd():
            with hole2.HoleAnalysis(universe) as h:
                h.run()
                h.run()  # run again for *.old files
                h.create_vmd_surface(filename='hole.vmd')

        sphpdbs = tmpdir.join('*.sph')
        assert len(glob.glob(str(sphpdbs))) == 0
        outfiles = tmpdir.join('*.out')
        assert len(glob.glob(str(outfiles))) == 0
        vdwradii = tmpdir.join('simple2.rad')
        assert len(glob.glob(str(vdwradii))) == 0
        pdbfiles = tmpdir.join('*.pdb')
        assert len(glob.glob(str(pdbfiles))) == 0
        oldfiles = tmpdir.join('*.old')
        assert len(glob.glob(str(oldfiles))) == 0
        vmd_file = tmpdir.join('hole.vmd')
        assert len(glob.glob(str(vmd_file))) == 1

    def test_output_level(self, tmpdir, universe):
        with tmpdir.as_cwd():
            with pytest.warns(UserWarning, match='needs to be < 3'):
                h = hole2.HoleAnalysis(universe,
                                       output_level=100)
                h.run(start=self.start,
                      stop=self.stop, random_seed=self.random_seed)

            # no profiles
            assert len(h.profiles) == 0

    def test_cpoint_geometry(self, tmpdir, universe):
        protein = universe.select_atoms('protein')
        cogs = [protein.center_of_geometry() for ts in universe.trajectory]
        with tmpdir.as_cwd():
            h = hole2.HoleAnalysis(universe,
                                   select='protein',
                                   cpoint='center_of_geometry',
                                   write_input_files=True)
            h.run(start=self.start,
                  stop=self.stop, random_seed=self.random_seed)

            infiles = sorted(glob.glob(str(tmpdir.join('hole*.inp'))))
            for file, cog in zip(infiles, cogs[self.start:self.stop]):
                with open(file, 'r') as f:
                    line = f.read().split('CPOINT')[1].split('\n')[0]
                arr = np.array(list(map(float, line.split())))
                assert_almost_equal(arr, cog)

    # plotting
    def test_plot(self, hole, frames, profiles):
        ax = hole.plot(label=True, frames=None, y_shift=1)
        err_msg = "HoleAnalysis.plot() did not produce an Axes instance"
        assert isinstance(ax, matplotlib.axes.Axes), err_msg
        lines = ax.get_lines()[:]
        assert len(lines) == hole.n_frames
        for i, (line, frame, profile) in enumerate(zip(lines, frames, profiles)):
            x, y = line.get_data()
            assert_almost_equal(x, profile.rxn_coord)
            assert_almost_equal(y, profile.radius + i)
            assert line.get_label() == str(frame)

    def test_plot_mean_profile(self, hole, frames, profiles):
        binned, bins = hole.bin_radii(bins=100)
        mean = np.array(list(map(np.mean, binned)))
        stds = np.array(list(map(np.std, binned)))
        midpoints = 0.5 * bins[1:] + bins[:-1]
        ylow = list(mean-(2*stds))
        yhigh = list(mean+(2*stds))

        ax = hole.plot_mean_profile(bins=100, n_std=2)

        # test fillbetween standard deviation
        children = ax.get_children()
        poly = []
        for x in children:
            if isinstance(x, matplotlib.collections.PolyCollection):
                poly.append(x)
        assert len(poly) == 1
        xp, yp = poly[0].get_paths()[0].vertices.T
        assert_almost_equal(np.unique(xp), np.unique(midpoints))
        assert_almost_equal(np.unique(yp), np.unique(ylow+yhigh))

        # test mean line
        lines = ax.get_lines()
        assert len(lines) == 1
        xl, yl = lines[0].get_data()
        assert_almost_equal(xl, midpoints)
        assert_almost_equal(yl, mean)

    @pytest.mark.skipif(sys.version_info < (3, 1),
                        reason="get_data_3d requires 3.1 or higher")
    def test_plot3D(self, hole, frames, profiles):
        ax = hole.plot3D(frames=None, r_max=None)
        err_msg = "HoleAnalysis.plot3D() did not produce an Axes3D instance"
        assert isinstance(ax, mpl_toolkits.mplot3d.Axes3D), err_msg
        lines = ax.get_lines()[:]
        assert len(lines) == hole.n_frames

        for line, frame, profile in zip(lines, frames, profiles):
            x, y, z = line.get_data_3d()
            assert_almost_equal(x, profile.rxn_coord)
            assert_almost_equal(np.unique(y), [frame])
            assert_almost_equal(z, profile.radius)
            assert line.get_label() == str(frame)
        
    @pytest.mark.skipif(sys.version_info < (3, 1),
                        reason="get_data_3d requires 3.1 or higher")
    def test_plot3D_rmax(self, hole, frames, profiles):
        ax = hole.plot3D(r_max=2.5)
        err_msg = "HoleAnalysis.plot3D(rmax=float) did not produce an Axes3D instance"
        assert isinstance(ax, mpl_toolkits.mplot3d.Axes3D), err_msg

        lines = ax.get_lines()[:]

        for line, frame, profile in zip(lines, frames, profiles):
            x, y, z = line.get_data_3d()
            assert_almost_equal(x, profile.rxn_coord)
            assert_almost_equal(np.unique(y), [frame])
            radius = np.where(profile.radius > 2.5, np.nan, profile.radius)
            assert_almost_equal(z, radius)
            assert line.get_label() == str(frame)

    @pytest.mark.skipif(sys.version_info > (3, 1),
                        reason="get_data_3d requires 3.1 or higher")
    def test_plot3D(self, hole, frames, profiles):
        ax = hole.plot3D(frames=None, r_max=None)
        err_msg = "HoleAnalysis.plot3D() did not produce an Axes3D instance"
        assert isinstance(ax, mpl_toolkits.mplot3d.Axes3D), err_msg
        lines = ax.get_lines()[:]
        assert len(lines) == hole.n_frames

        for line, frame, profile in zip(lines, frames, profiles):
            x, y = line.get_data()
            assert_almost_equal(x, profile.rxn_coord)
            assert_almost_equal(np.unique(y), [frame])
            assert line.get_label() == str(frame)

    @pytest.mark.skipif(sys.version_info > (3, 1),
                        reason="get_data_3d requires 3.1 or higher")
    def test_plot3D_rmax(self, hole, frames, profiles):
        ax = hole.plot3D(r_max=2.5)
        err_msg = "HoleAnalysis.plot3D(rmax=float) did not produce an Axes3D instance"
        assert isinstance(ax, mpl_toolkits.mplot3d.Axes3D), err_msg

        lines = ax.get_lines()[:]

        for line, frame, profile in zip(lines, frames, profiles):
            x, y = line.get_data()
            assert_almost_equal(x, profile.rxn_coord)
            assert_almost_equal(np.unique(y), [frame])
            assert line.get_label() == str(frame)


class TestHoleAnalysisLong(BaseTestHole):

    start = 0
    stop = 11

    rmsd = np.array([6.10501252e+00, 4.88398472e+00, 3.66303524e+00, 2.44202454e+00,
                     1.22100521e+00, 1.67285541e-07, 1.22100162e+00, 2.44202456e+00,
                     3.66303410e+00, 4.88398478e+00, 6.10502262e+00])

    @pytest.fixture
    def order_parameter_keys_values(self, hole):
        op = hole.over_order_parameters(self.rmsd, frames=None)
        return op.keys(), op.values()

    def test_gather(self, hole):
        gd = hole.gather(flat=False)
        for i, p in enumerate(hole.profiles.values()):
            assert_almost_equal(p.rxn_coord, gd['rxn_coord'][i])
            assert_almost_equal(p.radius, gd['radius'][i])
            assert_almost_equal(p.cen_line_D, gd['cen_line_D'][i])

    def test_gather_flat(self, hole):
        gd = hole.gather(flat=True)
        i = 0
        for p in hole.profiles.values():
            j = i+len(p.rxn_coord)
            assert_almost_equal(p.rxn_coord, gd['rxn_coord'][i:j])
            assert_almost_equal(p.radius, gd['radius'][i:j])
            assert_almost_equal(p.cen_line_D, gd['cen_line_D'][i:j])
            i = j
        assert_equal(i, len(gd['rxn_coord']))

    def test_min_radius(self, hole):
        rad = hole.min_radius()
        for (f1, p), (f2, r) in zip(hole.profiles.items(), rad):
            assert_equal(f1, f2)
            assert_almost_equal(min(p.radius), r)

    def test_over_order_parameters(self, hole):
        op = self.rmsd
        profiles = hole.over_order_parameters(op, frames=None)
        assert len(op) == len(profiles)

        for key, rmsd in zip(profiles.keys(), np.sort(op)):
            assert key == rmsd

        idx = np.argsort(op)
        arr = np.array(list(hole.profiles.values()), dtype=object)
        for op_prof, arr_prof in zip(profiles.values(), arr[idx]):
            assert op_prof is arr_prof

    def test_over_order_parameters_file(self, hole, tmpdir):
        op = self.rmsd
        with tmpdir.as_cwd():
            np.savetxt('rmsd.dat', self.rmsd)
            profiles = hole.over_order_parameters('rmsd.dat', frames=None)

        assert len(op) == len(profiles)

        for key, rmsd in zip(profiles.keys(), np.sort(op)):
            assert key == rmsd

        idx = np.argsort(op)
        arr = np.array(list(hole.profiles.values()), dtype=object)
        for op_prof, arr_prof in zip(profiles.values(), arr[idx]):
            assert op_prof is arr_prof

    def test_over_order_parameters_missing_file(self, hole):
        with pytest.raises(ValueError) as exc:
            hole.over_order_parameters('missing.dat')
        assert 'not found' in str(exc.value)

    def test_over_order_parameters_invalid_file(self, hole):
        with pytest.raises(ValueError) as exc:
            hole.over_order_parameters(PDB_HOLE)
        assert 'Could not parse' in str(exc.value)

    def test_over_order_parameters_frames(self, hole):
        op = self.rmsd
        n_frames = 7
        profiles = hole.over_order_parameters(op, frames=np.arange(n_frames))
        assert len(profiles) == n_frames
        for key, rmsd in zip(profiles.keys(), np.sort(op[:n_frames])):
            assert key == rmsd

        idx = np.argsort(op[:n_frames])
        values = list(hole.profiles.values())[:n_frames]
        arr = np.array(values, dtype=object)
        for op_prof, arr_prof in zip(profiles.values(), arr[idx]):
            assert op_prof is arr_prof

    def test_bin_radii(self, hole):
        radii, bins = hole.bin_radii(bins=100)
        dct = hole.gather(flat=True)
        coords = dct['rxn_coord']

        assert len(bins) == 101
        assert_almost_equal(bins[0], coords.min())
        assert_almost_equal(bins[-1], coords.max())
        assert len(radii) == (len(bins)-1)

        # check first frame profile
        first = hole.profiles[0]
        for row in first:
            coord = row.rxn_coord
            rad = row.radius
            for i, (lower, upper) in enumerate(zip(bins[:-1], bins[1:])):
                if coord > lower and coord <= upper:
                    assert rad in radii[i]
                    break
            else:
                raise AssertionError('Radius not in binned radii')

    @pytest.mark.parametrize('midpoint', [1.5, 1.8, 2.0, 2.5])
    def test_bin_radii_range(self, hole, midpoint):
        radii, bins = hole.bin_radii(bins=100, 
                                     range=(midpoint, midpoint))
        dct = hole.gather(flat=True)
        coords = dct['rxn_coord']

        assert len(bins) == 101
        low = midpoint - 0.5
        high = midpoint + 0.5
        assert_almost_equal(bins[0], low)
        assert_almost_equal(bins[-1], high)
        assert len(radii) == (len(bins)-1)

        # check first frame profile
        first = hole.profiles[0]
        for row in first:
            coord = row.rxn_coord
            rad = row.radius
            if coord > low and coord <= high:
                for i, (lower, upper) in enumerate(zip(bins[:-1], bins[1:])):
                    if coord > lower and coord <= upper:
                        assert rad in radii[i]
                        break
                else:
                    raise AssertionError('Radius not in binned radii')
            else:
                assert not any([rad in x for x in radii])

    def test_bin_radii_edges(self, hole):
        brange = list(np.linspace(1.0, 2.0, num=101, endpoint=True))
        moved = brange[30:] + brange[10:30] + brange[:10]
        e_radii, e_bins = hole.bin_radii(bins=moved, range=(0.0, 0.0))
        r_radii, r_bins = hole.bin_radii(bins=100, range=(1.5, 1.5))
        assert_almost_equal(e_bins, r_bins)
        for e, r in zip(e_radii, r_radii):
            assert_almost_equal(e, r)
        
    def test_histogram_radii(self, hole):
        means, _ = hole.histogram_radii(aggregator=np.mean,
                                        bins=100)
        radii, _ = hole.bin_radii(bins=100)
        assert means.shape == (100,)
        for r, m in zip(radii, means):
            assert_almost_equal(r.mean(), m)

    # plotting

    def test_plot_select_frames(self, hole, frames, profiles):
        ax = hole.plot(label=True, frames=[2, 3], y_shift=1)
        err_msg = "HoleAnalysis.plot() did not produce an Axes instance"
        assert isinstance(ax, matplotlib.axes.Axes), err_msg
        lines = ax.get_lines()[:]
        assert len(lines) == 2
        for i, (line, frame, profile) in enumerate(zip(lines, frames[2:4], profiles[2:4])):
            x, y = line.get_data()
            assert_almost_equal(x, profile.rxn_coord)
            assert_almost_equal(y, profile.radius + i)
            assert line.get_label() == str(frame)

    @pytest.mark.parametrize('agg', [np.max, np.mean, np.std, np.min])
    def test_plot_order_parameters(self, hole, order_parameter_keys_values,
                                   agg):
        opx = np.array(list(order_parameter_keys_values[0]))
        opy = np.array([agg(p.radius) for p in order_parameter_keys_values[1]])
        ax = hole.plot_order_parameters(self.rmsd, aggregator=agg, frames=None)
        err_msg = ("HoleAnalysis.plot_order_parameters()"
                   "did not produce an Axes instance")
        assert isinstance(ax, matplotlib.axes.Axes), err_msg

        lines = ax.get_lines()
        assert len(lines) == 1
        x, y = lines[0].get_data()
        assert_almost_equal(x, opx)
        assert_almost_equal(y, opy)

    @pytest.mark.skipif(sys.version_info < (3, 1), 
                        reason="get_data_3d requires 3.1 or higher")
    def test_plot3D_order_parameters(self, hole, order_parameter_keys_values):
        opx = np.array(list(order_parameter_keys_values[0]))
        profiles = np.array(list(order_parameter_keys_values[1]))
        ax = hole.plot3D_order_parameters(self.rmsd, frames=None)
        err_msg = ("HoleAnalysis.plot3D_order_parameters() "
                   "did not produce an Axes3D instance")
        assert isinstance(ax, mpl_toolkits.mplot3d.Axes3D), err_msg

        lines = ax.get_lines()
        assert len(lines) == hole.n_frames
        for line, opx_, profile in zip(lines, opx, profiles):
            x, y, z = line.get_data_3d()
            assert_almost_equal(x, profile.rxn_coord)
            assert_almost_equal(np.unique(y), np.array([opx_]))
            assert_almost_equal(z, profile.radius)

    @pytest.mark.skipif(sys.version_info > (3, 1), 
                        reason="get_data_3d requires 3.1 or higher")
    def test_plot3D_order_parameters(self, hole, order_parameter_keys_values):
        opx = np.array(list(order_parameter_keys_values[0]))
        profiles = np.array(list(order_parameter_keys_values[1]))
        ax = hole.plot3D_order_parameters(self.rmsd, frames=None)
        err_msg = ("HoleAnalysis.plot3D_order_parameters() "
                   "did not produce an Axes3D instance")
        assert isinstance(ax, mpl_toolkits.mplot3d.Axes3D), err_msg

        lines = ax.get_lines()
        assert len(lines) == hole.n_frames
        for line, opx_, profile in zip(lines, opx, profiles):
            x, y = line.get_data()
            assert_almost_equal(x, profile.rxn_coord)
            assert_almost_equal(np.unique(y), np.array([opx_]))

@pytest.mark.skipif(executable_not_found("hole"),
                    reason="Test skipped because HOLE not found")
class TestHoleModule(object):
    try:
        # on Unix we can manipulate our limits: http://docs.python.org/2/library/resource.html
        import resource
        soft_max_open_files, hard_max_open_files = resource.getrlimit(
            resource.RLIMIT_NOFILE)
    except ImportError:
        pass

    @staticmethod
    @pytest.fixture()
    def universe():
        return mda.Universe(MULTIPDB_HOLE)

    @pytest.mark.skipif(rlimits_missing,
                        reason="Test skipped because platform does not allow setting rlimits")
    def test_hole_module_fd_closure(self, universe, tmpdir):
        """test open file descriptors are closed (MDAnalysisTests.analysis.test_hole.TestHoleModule): Issue 129"""
        # If Issue 129 isn't resolved, this function will produce an OSError on
        # the system, and cause many other tests to fail as well.
        #
        # Successful test takes ~10 s, failure ~2 s.

        # Hasten failure by setting "ulimit -n 64" (can't go too low because of open modules etc...)
        import resource

        # ----- temporary hack -----
        # on Mac OS X (on Travis) we run out of open file descriptors
        # before even starting this test (see
        # https://github.com/MDAnalysis/mdanalysis/pull/901#issuecomment-231938093);
        # if this issue is solved by #363 then revert the following
        # hack:
        #
        import platform
        if platform.platform() == "Darwin":
            max_open_files = 512
        else:
            max_open_files = 64
        #
        # --------------------------

        resource.setrlimit(resource.RLIMIT_NOFILE,
                           (max_open_files, self.hard_max_open_files))

        with tmpdir.as_cwd():
            try:
                H = hole2.HoleAnalysis(universe, cvect=[0, 1, 0], sample=20.0)
            finally:
                self._restore_rlimits()

            # pretty unlikely that the code will get through 2 rounds if the MDA
            # issue 129 isn't fixed, although this depends on the file descriptor
            # open limit for the machine in question
            try:
                for i in range(2):
                    # will typically get an OSError for too many files being open after
                    # about 2 seconds if issue 129 isn't resolved
                    H.run()
            except OSError as err:
                if err.errno == errno.EMFILE:
                    raise pytest.fail(
                        "hole2.HoleAnalysis does not close file descriptors (Issue 129)")
                raise
            finally:
                # make sure to restore open file limit !!
                self._restore_rlimits()

    def _restore_rlimits(self):
        try:
            import resource
            resource.setrlimit(resource.RLIMIT_NOFILE,
                               (self.soft_max_open_files, self.hard_max_open_files))
        except ImportError:
            pass
