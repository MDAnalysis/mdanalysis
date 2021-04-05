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
import numpy as np
import pytest
import sys
import warnings
from unittest.mock import Mock, patch

from numpy.testing import assert_equal, assert_almost_equal

import gridData.OpenDX

import MDAnalysis as mda
from MDAnalysis.analysis import density

from MDAnalysisTests.datafiles import TPR, XTC, GRO, PDB_full
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
    selections = {'none': "resname None",
                  'static': "name OW",
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
        return mda.Universe(self.topology, self.trajectory, tpr_resid_from_one=False)

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
        regex = (r"Box padding \(currently set at 1\.0\) is not used "
                 r"in user defined grids\.")
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
        with pytest.raises(ValueError, match="Gridcenter must be a 3D coordinate"):
            D = density.DensityAnalysis(
                universe.select_atoms(self.selections['static']),
                delta=self.delta, xdim=10.0, ydim=10.0, zdim=10.0,
                gridcenter=self.gridcenters['error1']).run(step=5)

    def test_ValueError_userdefn_gridcenter_type(self, universe):
        # Test gridcenter includes non-numeric strings
        with pytest.raises(ValueError, match="Gridcenter must be a 3D coordinate"):
            D = density.DensityAnalysis(
                universe.select_atoms(self.selections['static']),
                delta=self.delta, xdim=10.0, ydim=10.0, zdim=10.0,
                gridcenter=self.gridcenters['error2']).run(step=5)

    def test_ValueError_userdefn_gridcenter_missing(self, universe):
        # Test no gridcenter provided when grid dimensions are given
        regex = ("Gridcenter or grid dimensions are not provided")
        with pytest.raises(ValueError, match=regex):
            D = density.DensityAnalysis(
                universe.select_atoms(self.selections['static']),
                delta=self.delta, xdim=10.0, ydim=10.0, zdim=10.0).run(step=5)

    def test_ValueError_userdefn_xdim_type(self, universe):
        # Test xdim != int or float
        with pytest.raises(ValueError, match="xdim, ydim, and zdim must be numbers"):
            D = density.DensityAnalysis(
                universe.select_atoms(self.selections['static']),
                delta=self.delta, xdim="MDAnalysis", ydim=10.0, zdim=10.0,
                gridcenter=self.gridcenters['static_defined']).run(step=5)

    def test_ValueError_userdefn_xdim_nanvalue(self, universe):
        # Test  xdim set to NaN value
        regex = ("Gridcenter or grid dimensions have NaN element")
        with pytest.raises(ValueError, match=regex):
            D = density.DensityAnalysis(
                universe.select_atoms(self.selections['static']),
                delta=self.delta, xdim=np.NaN, ydim=10.0, zdim=10.0,
                gridcenter=self.gridcenters['static_defined']).run(step=5)

    def test_warn_noatomgroup(self, universe):
        regex = ("No atoms in AtomGroup at input time frame. "
                 "This may be intended; please ensure that "
                 "your grid selection covers the atomic "
                 "positions you wish to capture.")
        with pytest.warns(UserWarning, match=regex):
            D = density.DensityAnalysis(
                universe.select_atoms(self.selections['none']),
                delta=self.delta, xdim=1.0, ydim=2.0, zdim=2.0, padding=0.0,
                gridcenter=self.gridcenters['static_defined']).run(step=5)

    def test_ValueError_noatomgroup(self, universe):
        with pytest.raises(ValueError, match="No atoms in AtomGroup at input"
                                             " time frame. Grid for density"
                                             " could not be automatically"
                                             " generated. If this is"
                                             " expected, a user"
                                             " defined grid will "
                                             "need to be provided instead."):
            D = density.DensityAnalysis(
                universe.select_atoms(self.selections['none'])).run(step=5)

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
