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
from __future__ import print_function, division, absolute_import


from six.moves import zip
import numpy as np
import pytest
import sys
import warnings

from numpy.testing import assert_equal, assert_almost_equal

import gridData.OpenDX

import MDAnalysis as mda
from MDAnalysis.analysis import density

from MDAnalysisTests.datafiles import TPR, XTC, GRO, PDB_full
from mock import Mock, patch
from MDAnalysisTests.util import block_import


class TestDensity(object):
    nbins = 3, 4, 5
    counts = 100
    Lmax = 10.

    @pytest.fixture(scope='class')
    def bins(self):
        return [np.linspace(0, self.Lmax, n + 1) for n in self.nbins]

    @pytest.fixture()
    def h_and_edges(self, bins):
        return np.histogramdd(
            self.Lmax * np.sin(
                np.linspace(0, 1, self.counts * 3)).reshape(self.counts, 3),
            bins=bins)

    @pytest.fixture()
    def D(self, h_and_edges):
        h, edges = h_and_edges
        d = density.Density(h, edges, parameters={'isDensity': False},
                            units={'length': 'A'})
        d.make_density()
        return d

    def test_shape(self, D):
        assert D.grid.shape == self.nbins

    def test_edges(self, bins, D):
        for dim, (edges, fixture) in enumerate(zip(D.edges, bins)):
            assert_almost_equal(edges, fixture,
                                err_msg="edges[{0}] mismatch".format(dim))

    def test_midpoints(self, bins, D):
        midpoints = [0.5*(b[:-1] + b[1:]) for b in bins]
        for dim, (mp, fixture) in enumerate(zip(D.midpoints, midpoints)):
            assert_almost_equal(mp, fixture,
                                err_msg="midpoints[{0}] mismatch".format(dim))

    def test_delta(self, D):
        deltas = np.array([self.Lmax])/np.array(self.nbins)
        assert_almost_equal(D.delta, deltas)

    def test_grid(self, D):
        dV = D.delta.prod()  # orthorhombic grids only!
        # counts = (rho[0] * dV[0] + rho[1] * dV[1] ...) = sum_i rho[i] * dV
        assert_almost_equal(D.grid.sum() * dV, self.counts)

    def test_origin(self, bins, D):
        midpoints = [0.5*(b[:-1] + b[1:]) for b in bins]
        origin = [m[0] for m in midpoints]
        assert_almost_equal(D.origin, origin)

    def test_check_set_unit_keyerror(self, D):
        units = {'weight': 'A'}
        with pytest.raises(ValueError):
            D._check_set_unit(units)

    def test_check_set_unit_attributeError(self, D):
        units = []
        with pytest.raises(ValueError):
            D._check_set_unit(units)

    def test_check_set_unit_nolength(self, D):
        del D.units['length']
        units = {'density': 'A^{-3}'}
        with pytest.raises(ValueError):
            D._check_set_unit(units)

    @pytest.mark.parametrize('dxtype',
                             ("float", "double", "int", "byte"))
    def test_export_types(self, D, dxtype, tmpdir, outfile="density.dx"):
        with tmpdir.as_cwd():
            D.export(outfile, type=dxtype)

            dx = gridData.OpenDX.field(0)
            dx.read(outfile)
            data = dx.components['data']
        assert data.type == dxtype

class DensityParameters(object):
    topology = TPR
    trajectory = XTC
    delta = 2.0
    selections = {'static': "name OW",
                  'dynamic': "name OW and around 4 (protein and resid 1-10)",
                  'solute': "protein and not name H*",
                  }
    references = {'static':
                      {'meandensity': 0.016764271713091212, },
                  'static_sliced':
                      {'meandensity': 0.016764270747693617, },
                  'static_defined':
                      {'meandensity': 0.0025000000000000005, },
                  'static_defined_unequal':
                      {'meandensity': 0.006125, },
                  'dynamic':
                      {'meandensity': 0.0012063418843728784, },
                  'notwithin':
                      {'meandensity': 0.015535385132107926, },
                  }
    cutoffs = {'notwithin': 4.0, }
    gridcenters = {'static_defined': np.array([56.0, 45.0, 35.0]),
                   'error1': np.array([56.0, 45.0]),
                   'error2': [56.0, 45.0, "MDAnalysis"],
                   }
    precision = 5
    outfile = 'density.dx'

    @pytest.fixture()
    def universe(self):
        return mda.Universe(self.topology, self.trajectory)

class TestDensityAnalysis(DensityParameters):
    def check_DensityAnalysis(self, ag, ref_meandensity,
                                    tmpdir, runargs=None, **kwargs):
        runargs = runargs if runargs else {}
        with tmpdir.as_cwd():
            D = density.DensityAnalysis(ag, delta=self.delta, **kwargs).run(**runargs)
            assert_almost_equal(D.density.grid.mean(), ref_meandensity,
                                err_msg="mean density does not match")
            D.density.export(self.outfile)

            D2 = density.Density(self.outfile)
            assert_almost_equal(D.density.grid, D2.grid, decimal=self.precision,
                                err_msg="DX export failed: different grid sizes")

    @pytest.mark.parametrize("mode", ("static", "dynamic"))
    def test_run(self, mode, universe, tmpdir):
        updating = (mode == "dynamic")
        self.check_DensityAnalysis(
            universe.select_atoms(self.selections[mode], updating=updating),
            self.references[mode]['meandensity'],
            tmpdir=tmpdir
        )

    def test_sliced(self, universe, tmpdir):
        self.check_DensityAnalysis(
            universe.select_atoms(self.selections['static']),
            self.references['static_sliced']['meandensity'],
            tmpdir=tmpdir,
            runargs=dict(start=1, stop=-1, step=2),
        )

    def test_userdefn_eqbox(self, universe, tmpdir):
        with warnings.catch_warnings():
            # Do not need to see UserWarning that box is too small
            warnings.simplefilter("ignore")
            self.check_DensityAnalysis(
                universe.select_atoms(self.selections['static']),
                self.references['static_defined']['meandensity'],
                tmpdir=tmpdir,
                gridcenter=self.gridcenters['static_defined'],
                xdim=10.0,
                ydim=10.0,
                zdim=10.0,
            )

    def test_userdefn_neqbox(self, universe, tmpdir):
        self.check_DensityAnalysis(
            universe.select_atoms(self.selections['static']),
            self.references['static_defined_unequal']['meandensity'],
            tmpdir=tmpdir,
            gridcenter=self.gridcenters['static_defined'],
            xdim=10.0,
            ydim=15.0,
            zdim=20.0,
        )

    def test_userdefn_boxshape(self, universe):
        D = density.DensityAnalysis(
            universe.select_atoms(self.selections['static']),
            delta=1.0, xdim=8.0, ydim=12.0, zdim=17.0,
            gridcenter=self.gridcenters['static_defined']).run()
        assert D.density.grid.shape == (8, 12, 17)

    def test_warn_userdefn_padding(self, universe):
        regex = ("Box padding \(currently set at 1\.0\) is not used "
                 "in user defined grids\.")
        with pytest.warns(UserWarning, match=regex):
            D = density.DensityAnalysis(
                universe.select_atoms(self.selections['static']),
                delta=self.delta, xdim=100.0, ydim=100.0, zdim=100.0, padding=1.0,
                gridcenter=self.gridcenters['static_defined']).run(step=5)

    def test_warn_userdefn_smallgrid(self, universe):
        regex = ("Atom selection does not fit grid --- "
                 "you may want to define a larger box")
        with pytest.warns(UserWarning, match=regex):
            D = density.DensityAnalysis(
                universe.select_atoms(self.selections['static']),
                delta=self.delta, xdim=1.0, ydim=2.0, zdim=2.0, padding=0.0,
                gridcenter=self.gridcenters['static_defined']).run(step=5)

    def test_ValueError_userdefn_gridcenter_shape(self, universe):
        # Test len(gridcenter) != 3
        with pytest.raises(ValueError, match="gridcenter must be a 3D coordinate"):
            D = density.DensityAnalysis(
                universe.select_atoms(self.selections['static']),
                delta=self.delta, xdim=10.0, ydim=10.0, zdim=10.0,
                gridcenter=self.gridcenters['error1']).run(step=5)

    def test_ValueError_userdefn_gridcenter_type(self, universe):
        # Test gridcenter includes non-numeric strings
        with pytest.raises(ValueError, match="Non-number values assigned to gridcenter"):
            D = density.DensityAnalysis(
                universe.select_atoms(self.selections['static']),
                delta=self.delta, xdim=10.0, ydim=10.0, zdim=10.0,
                gridcenter=self.gridcenters['error2']).run(step=5)

    def test_ValueError_userdefn_xdim_type(self, universe):
        # Test xdim != int or float
        with pytest.raises(ValueError, match="xdim, ydim, and zdim must be numbers"):
            D = density.DensityAnalysis(
                universe.select_atoms(self.selections['static']),
                delta=self.delta, xdim="MDAnalysis", ydim=10.0, zdim=10.0,
                gridcenter=self.gridcenters['static_defined']).run(step=5)


# remove in 2.0.0
class Test_density_from_Universe(DensityParameters):
    def check_density_from_Universe(self, atomselection, ref_meandensity,
                                    universe, tmpdir, **kwargs):
        with tmpdir.as_cwd():
            D = density.density_from_Universe(
                universe, select=atomselection,
                delta=self.delta, **kwargs)
            assert_almost_equal(D.grid.mean(), ref_meandensity,
                                err_msg="mean density does not match")

            D.export(self.outfile)

            D2 = density.Density(self.outfile)
            assert_almost_equal(D.grid, D2.grid, decimal=self.precision,
                                err_msg="DX export failed: different grid sizes")

    def test_density_from_Universe(self, universe, tmpdir):
        self.check_density_from_Universe(
            self.selections['static'],
            self.references['static']['meandensity'],
            universe=universe,
            tmpdir=tmpdir
        )

    def test_sliced(self, universe, tmpdir):
        self.check_density_from_Universe(
            self.selections['static'],
            self.references['static_sliced']['meandensity'],
            start=1, stop=-1, step=2,
            universe=universe,
            tmpdir=tmpdir
            )

    def test_update_selection(self, universe, tmpdir):
        self.check_density_from_Universe(
            self.selections['dynamic'],
            self.references['dynamic']['meandensity'],
            update_selection=True,
            universe=universe,
            tmpdir=tmpdir
        )

    def test_notwithin(self, universe, tmpdir):
        self.check_density_from_Universe(
            self.selections['static'],
            self.references['notwithin']['meandensity'],
            soluteselection=self.selections['solute'],
            cutoff=self.cutoffs['notwithin'],
            universe=universe,
            tmpdir=tmpdir
        )

    def test_userdefn_eqbox(self, universe, tmpdir):
        self.check_density_from_Universe(
            self.selections['static'],
            self.references['static_defined']['meandensity'],
            gridcenter=self.gridcenters['static_defined'],
            xdim=10.0,
            ydim=10.0,
            zdim=10.0,
            universe=universe,
            tmpdir=tmpdir
        )

    def test_userdefn_neqbox(self, universe, tmpdir):
        self.check_density_from_Universe(
            self.selections['static'],
            self.references['static_defined_unequal']['meandensity'],
            gridcenter=self.gridcenters['static_defined'],
            xdim=10.0,
            ydim=15.0,
            zdim=20.0,
            universe=universe,
            tmpdir=tmpdir
        )

    def test_userdefn_boxshape(self, universe):
        D = density.density_from_Universe(
            universe, select=self.selections['static'],
            delta=1.0, xdim=8.0, ydim=12.0, zdim=17.0,
            gridcenter=self.gridcenters['static_defined'])
        assert D.grid.shape == (8, 12, 17)

    def test_userdefn_padding(self, universe):
        regex = ("Box padding \(currently set at 1\.0\) is not used "
                 "in user defined grids\.")
        with pytest.warns(UserWarning, match=regex):
            D = density.density_from_Universe(
                universe, select=self.selections['static'],
                delta=self.delta, xdim=100.0, ydim=100.0, zdim=100.0, padding=1.0,
                gridcenter=self.gridcenters['static_defined'])

    def test_userdefn_selwarning(self, universe):
        regex = ("Atom selection does not fit grid --- "
                "you may want to define a larger box")
        with pytest.warns(UserWarning, match=regex):
            D = density.density_from_Universe(
                universe, select=self.selections['static'],
                delta=self.delta, xdim=1.0, ydim=2.0, zdim=2.0, padding=0.0,
                gridcenter=self.gridcenters['static_defined'])

    def test_userdefn_ValueErrors(self, universe):
        # Test len(gridcenter) != 3
        with pytest.raises(ValueError):
            D = density.density_from_Universe(
                universe, select=self.selections['static'],
                delta=self.delta, xdim=10.0, ydim=10.0, zdim=10.0,
                gridcenter=self.gridcenters['error1'])
        # Test gridcenter includes non-numeric strings
        with pytest.raises(ValueError):
            D = density.density_from_Universe(
                universe, select=self.selections['static'],
                delta=self.delta, xdim=10.0, ydim=10.0, zdim=10.0,
                gridcenter=self.gridcenters['error2'])
        # Test xdim != int or float
        with pytest.raises(ValueError):
            D = density.density_from_Universe(
                universe, select=self.selections['static'],
                delta=self.delta, xdim="MDAnalysis", ydim=10.0, zdim=10.0,
                gridcenter=self.gridcenters['static_defined'])

    def test_has_DeprecationWarning(self, universe):
        with pytest.warns(DeprecationWarning,
                          match="will be removed in release 2.0.0"):
            D = density.density_from_Universe(
                universe, select=self.selections['static'],
                delta=self.delta)


class TestNotWithin(object):
    # tests notwithin_coordinates_factory
    # only checks that KDTree and distance_array give same results

    @staticmethod
    @pytest.fixture()
    def u():
        return mda.Universe(GRO)

    @pytest.mark.skipif(sys.maxsize <= 2**32,
                        reason="non-kdtree density too large for 32-bit")
    def test_within(self, u):
        vers1 = density.notwithin_coordinates_factory(u, 'resname SOL',
                                                      'protein', 2,
                                                      not_within=False,
                                                      use_kdtree=True)()
        vers2 = density.notwithin_coordinates_factory(u, 'resname SOL',
                                                      'protein', 2,
                                                      not_within=False,
                                                      use_kdtree=False)()

        assert_equal(vers1, vers2)

    @pytest.mark.skipif(sys.maxsize <= 2**32,
                        reason="non-kdtree density too large for 32-bit")
    def test_not_within(self, u):
        vers1 = density.notwithin_coordinates_factory(u, 'resname SOL',
                                                     'protein', 2,
                                                     not_within=True,
                                                     use_kdtree=True)()
        vers2 = density.notwithin_coordinates_factory(u, 'resname SOL',
                                                     'protein', 2,
                                                     not_within=True,
                                                     use_kdtree=False)()
        assert_equal(vers1, vers2)

    def test_has_DeprecationWarning(self, u):
        with pytest.warns(DeprecationWarning,
                          match="will be removed in release 2.0.0"):
            density.notwithin_coordinates_factory(u, 'resname SOL',
                                                  'protein', 2,
                                                  not_within=False,
                                                  use_kdtree=True)()


class TestBfactor2RMSF(object):
    @pytest.mark.parametrize('B', [20, 20.0, 20.123,
                                   np.array([45.678, 67.1, 0.0, 100.2])])
    def test_Bfactor2RMSF(self, B):
        values = density.Bfactor2RMSF(B)
        ref = np.sqrt(3. * B / 8.) / np.pi  # original equation
        assert_almost_equal(values, ref)

    def test_has_DeprecationWarning(self):
        with pytest.warns(DeprecationWarning,
                          match="will be removed in release 2.0.0"):
            density.Bfactor2RMSF(100.0)

class Test_density_from_PDB(object):

    @pytest.mark.parametrize('sigma,ref_shape,ref_gridsum',
                             [(None, (21, 26, 29), 108.710149),
                              (1.5,  (21, 26, 29), 160.050167)])
    def test_density_from_PDB(self, sigma, ref_shape, ref_gridsum):
        # simple test that the code runs: use B-factors from PDB or fixed sigma
        D = density.density_from_PDB(PDB_full, delta=2.0, sigma=sigma)

        assert isinstance(D, density.Density)
        assert_equal(D.grid.shape, ref_shape)
        assert_almost_equal(D.grid.sum(), ref_gridsum, decimal=5)

    def test_has_DeprecationWarning(self):
        with pytest.warns(DeprecationWarning,
                          match="will be removed in release 2.0.0"):
            density.density_from_PDB(PDB_full, delta=5.0, sigma=2)

class TestGridImport(object):

    @block_import('gridData')
    def test_absence_griddata(self):
        sys.modules.pop('MDAnalysis.analysis.density', None)
        # if gridData package is missing an ImportError should be raised
        # at the module level of MDAnalysis.analysis.density
        with pytest.raises(ImportError):
            import MDAnalysis.analysis.density

    def test_presence_griddata(self):
        sys.modules.pop('MDAnalysis.analysis.density', None)
        # no ImportError exception is raised when gridData is properly
        # imported by MDAnalysis.analysis.density

        # mock gridData in case there are testing scenarios where
        # it is not available
        mock = Mock()
        with patch.dict('sys.modules', {'gridData':mock}):
            try:
                import MDAnalysis.analysis.density
            except ImportError:
                pytest.fail(msg='''MDAnalysis.analysis.density should not raise
                             an ImportError if gridData is available.''')
