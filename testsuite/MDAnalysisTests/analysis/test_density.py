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

from MDAnalysisTests.datafiles import TPR, XTC, GRO
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


class Test_density_from_Universe(object):
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

    def check_density_from_Universe(self, atomselection, ref_meandensity,
                                    universe, tmpdir, **kwargs):
        import MDAnalysis.analysis.density

        with tmpdir.as_cwd():
            D = MDAnalysis.analysis.density.density_from_Universe(
                universe, select=atomselection,
                delta=self.delta, **kwargs)
            assert_almost_equal(D.grid.mean(), ref_meandensity,
                                err_msg="mean density does not match")

            D.export(self.outfile)

            D2 = MDAnalysis.analysis.density.Density(self.outfile)
            assert_almost_equal(D.grid, D2.grid, decimal=self.precision,
                                err_msg="DX export failed: different grid sizes")

    def test_density_from_Universe(self, universe, tmpdir):
        self.check_density_from_Universe(
            self.selections['static'],
            self.references['static']['meandensity'],
            universe=universe,
            tmpdir=tmpdir
        )

    def test_density_from_Universe_sliced(self, universe, tmpdir):
        self.check_density_from_Universe(
            self.selections['static'],
            self.references['static_sliced']['meandensity'],
            start=1, stop=-1, step=2,
            universe=universe,
            tmpdir=tmpdir
            )

    def test_density_from_Universe_update_selection(self, universe, tmpdir):
        self.check_density_from_Universe(
            self.selections['dynamic'],
            self.references['dynamic']['meandensity'],
            update_selection=True,
            universe=universe,
            tmpdir=tmpdir
        )

    def test_density_from_Universe_notwithin(self, universe, tmpdir):
        self.check_density_from_Universe(
            self.selections['static'],
            self.references['notwithin']['meandensity'],
            soluteselection=self.selections['solute'],
            cutoff=self.cutoffs['notwithin'],
            universe=universe,
            tmpdir=tmpdir
        )

    def test_density_from_Universe_userdefn_eqbox(self, universe, tmpdir):
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

    def test_density_from_Universe_userdefn_neqbox(self, universe, tmpdir):
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

    def test_density_from_Universe_userdefn_boxshape(self, universe):
        import MDAnalysis.analysis.density
        D = MDAnalysis.analysis.density.density_from_Universe(
            universe, select=self.selections['static'],
            delta=1.0, xdim=8.0, ydim=12.0, zdim=17.0,
            gridcenter=self.gridcenters['static_defined'])
        assert D.grid.shape == (8, 12, 17)

    def test_density_from_Universe_userdefn_padding(self, universe):
        import MDAnalysis.analysis.density
        wmsg = ("Box padding (currently set at 1.0) is not used in user "
                "defined grids.")
        with pytest.warns(UserWarning) as record:
            D = MDAnalysis.analysis.density.density_from_Universe(
                universe, select=self.selections['static'],
                delta=1.0, xdim=1.0, ydim=2.0, zdim=2.0, padding=1.0,
                gridcenter=self.gridcenters['static_defined'])
        
        assert len(record) == 3
        assert str(record[1].message.args[0]) == wmsg

    def test_density_from_Universe_userdefn_selwarning(self, universe):
        import MDAnalysis.analysis.density
        wmsg = ("Atom selection does not fit grid --- "
                "you may want to define a larger box")
        with pytest.warns(UserWarning) as record:
            D = MDAnalysis.analysis.density.density_from_Universe(
                universe, select=self.selections['static'],
                delta=1.0, xdim=1.0, ydim=2.0, zdim=2.0, padding=0.0,
                gridcenter=self.gridcenters['static_defined'])

        assert len(record) == 2
        assert str(record[1].message.args[0]) == wmsg

    def test_density_from_Universe_userdefn_ValueErrors(self, universe):
        import MDAnalysis.analysis.density
        # Test len(gridcenter) != 3
        with pytest.raises(ValueError):
            D = MDAnalysis.analysis.density.density_from_Universe(
                universe, select=self.selections['static'],
                delta=self.delta, xdim=10.0, ydim=10.0, zdim=10.0,
                gridcenter=self.gridcenters['error1'])
        # Test gridcenter includes non-numeric strings
        with pytest.raises(ValueError):
            D = MDAnalysis.analysis.density.density_from_Universe(
                universe, select=self.selections['static'],
                delta=self.delta, xdim=10.0, ydim=10.0, zdim=10.0,
                gridcenter=self.gridcenters['error2'])
        # Test xdim != int or float
        with pytest.raises(ValueError):
            D = MDAnalysis.analysis.density.density_from_Universe(
                universe, select=self.selections['static'],
                delta=self.delta, xdim="MDAnalysis", ydim=10.0, zdim=10.0,
                gridcenter=self.gridcenters['static_defined'])


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


class TestNotWithin(object):
    # tests notwithin_coordinates_factory
    # only checks that KDTree and distance_array give same results

    @staticmethod
    @pytest.fixture()
    def u():
        return mda.Universe(GRO)

    def test_within(self, u):
        from MDAnalysis.analysis.density import notwithin_coordinates_factory as ncf

        vers1 = ncf(u, 'resname SOL', 'protein', 2, not_within=False, use_kdtree=True)()
        vers2 = ncf(u, 'resname SOL', 'protein', 2, not_within=False, use_kdtree=False)()

        assert_equal(vers1, vers2)

    def test_not_within(self, u):
        from MDAnalysis.analysis.density import notwithin_coordinates_factory as ncf

        vers1 = ncf(u, 'resname SOL', 'protein', 2, not_within=True, use_kdtree=True)()
        vers2 = ncf(u, 'resname SOL', 'protein', 2, not_within=True, use_kdtree=False)()

        assert_equal(vers1, vers2)
