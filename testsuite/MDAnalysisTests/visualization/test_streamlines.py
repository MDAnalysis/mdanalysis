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
from numpy.testing import assert_almost_equal
import MDAnalysis
from MDAnalysis.visualization import (streamlines,
                                      streamlines_3D)
from MDAnalysis.coordinates.XTC import XTCWriter
from MDAnalysisTests.datafiles import Martini_membrane_gro
import pytest
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import os

@pytest.fixture(scope="session")
def univ():
    u = MDAnalysis.Universe(Martini_membrane_gro)
    return u

@pytest.fixture(scope="session")
def membrane_xtc(tmpdir_factory, univ):
    x_delta, y_delta, z_delta  = 0.5, 0.3, 0.2
    tmp_xtc = tmpdir_factory.mktemp('streamlines').join('dummy.xtc')

    with XTCWriter(str(tmp_xtc), n_atoms=univ.atoms.n_atoms) as xtc_writer:
        for i in range(5):
           univ.atoms.translate([x_delta, y_delta, z_delta])
           xtc_writer.write(univ.atoms)
           x_delta += 0.1
           y_delta += 0.08
           z_delta += 0.02
    return str(tmp_xtc)

def test_streamplot_2D(membrane_xtc, univ, tmpdir):
    # simple roundtrip test to ensure that
    # a plot is generated by the 2D streamplot
    # code
    u1, v1, avg, std = streamlines.generate_streamlines(topology_file_path=Martini_membrane_gro,
                                                        trajectory_file_path=membrane_xtc,
                                                        grid_spacing=20,
                                                        MDA_selection='name PO4',
                                                        start_frame=1,
                                                        end_frame=2,
                                                        xmin=univ.atoms.positions[...,0].min(),
                                                        xmax=univ.atoms.positions[...,0].max(),
                                                        ymin=univ.atoms.positions[...,1].min(),
                                                        ymax=univ.atoms.positions[...,1].max(),
                                                        maximum_delta_magnitude=2.0,
                                                        num_cores=1)
    x = np.linspace(univ.atoms.positions[...,0].min(),
                    univ.atoms.positions[...,0].max(),
                    5)
    y = np.linspace(univ.atoms.positions[...,1].min(),
                    univ.atoms.positions[...,1].max(),
                    5)
    speed = np.sqrt(u1*u1 + v1*v1)
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')
    ax.streamplot(x, y, u1, v1, density=(10,10), color=speed, linewidth=3*speed/speed.max(),
                  cmap=plt.cm.viridis)
    plot_outpath = os.path.join(str(tmpdir), 'testing_streamline.png')
    fig.savefig(plot_outpath, dpi=300)

    with open(plot_outpath, 'rb'):
        pass

def test_streamplot_3D(membrane_xtc, univ, tmpdir):
    # because mayavi is too heavy of a dependency
    # for a roundtrip plotting test, simply
    # aim to check for sensible values
    # returned by generate_streamlines_3d
    dx, dy, dz = streamlines_3D.generate_streamlines_3d(topology_file_path=Martini_membrane_gro,
                                                        trajectory_file_path=membrane_xtc,
                                                        grid_spacing=20,
                                                        MDA_selection='name PO4',
                                                        start_frame=1,
                                                        end_frame=2,
                                                        xmin=univ.atoms.positions[...,0].min(),
                                                        xmax=univ.atoms.positions[...,0].max(),
                                                        ymin=univ.atoms.positions[...,1].min(),
                                                        ymax=univ.atoms.positions[...,1].max(),
                                                        zmin=univ.atoms.positions[...,2].min(),
                                                        zmax=univ.atoms.positions[...,2].max(),
                                                        maximum_delta_magnitude=2.0,
                                                        num_cores=1)
    assert dx.shape == (5, 5, 2)
    assert dy.shape == (5, 5, 2)
    assert dz.shape == (5, 5, 2)
    assert_almost_equal(dx[4, 4, 0], 0.700004, decimal=5)
    assert_almost_equal(dy[0, 0, 0], 0.460000, decimal=5)
    assert_almost_equal(dz[2, 2, 0], 0.240005, decimal=5)
