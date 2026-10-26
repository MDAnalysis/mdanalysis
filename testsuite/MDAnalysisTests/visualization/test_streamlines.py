# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
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
from numpy.testing import assert_allclose
import MDAnalysis
from MDAnalysis.visualization import streamlines, streamlines_3D
from MDAnalysis.coordinates.XTC import XTCWriter
from MDAnalysisTests.datafiles import Martini_membrane_gro
import pytest
from pytest import approx
import matplotlib.pyplot as plt
import os


@pytest.fixture(scope="session")
def univ():
    u = MDAnalysis.Universe(Martini_membrane_gro)
    return u


@pytest.fixture(scope="session")
def membrane_xtc(tmpdir_factory, univ):
    x_delta, y_delta, z_delta = 0.5, 0.3, 0.2
    tmp_xtc = tmpdir_factory.mktemp("streamlines").join("dummy.xtc")

    with XTCWriter(str(tmp_xtc), n_atoms=univ.atoms.n_atoms) as xtc_writer:
        for i in range(5):
            univ.atoms.translate([x_delta, y_delta, z_delta])
            xtc_writer.write(univ.atoms)
            x_delta += 0.1
            y_delta += 0.08
            z_delta += 0.02
    return str(tmp_xtc)


def test_produce_list_indices_point_in_polygon_this_frame():
    # Define two squares:
    # square1 covers the area [0,1]x[0,1]
    # square2 covers the area [2,3]x[2,3]
    square1 = [(0, 0), (1, 0), (1, 1), (0, 1)]
    square2 = [(2, 2), (3, 2), (3, 3), (2, 3)]

    # Create a list of vertex coordinates (for two squares)
    vertex_list = [square1, square2]

    # Define points:
    # Point [0.5, 0.5] lies inside square1.
    # Point [1.5, 1.5] lies outside both squares.
    # Point [2.5, 2.5] lies inside square2.
    # Point [3.5, 3.5] lies outside both squares.
    points = np.array([[0.5, 0.5], [1.5, 1.5], [2.5, 2.5], [3.5, 3.5]])

    # Call the function under test.
    result = streamlines._produce_list_indices_point_in_polygon_this_frame(
        vertex_list, points
    )

    # np.where returns a tuple; thus for each square we expect:
    # For square1: point index 0 is inside → (array([0]),)
    # For square2: point index 2 is inside → (array([2]),)
    expected = [(np.array([0]),), (np.array([2]),)]

    # Check that each result matches the expected indices.
    for res_tuple, exp_tuple in zip(result, expected):
        np.testing.assert_array_equal(res_tuple[0], exp_tuple[0])


def test_produce_list_centroids_empty():
    # Simulate an empty index set for one square:
    list_indices = [(np.array([]),)]
    # Dummy particle coordinate array (won't be used since indices is empty)
    pts = np.array([[0, 0], [1, 1]])
    result = streamlines._produce_list_centroids_this_frame(list_indices, pts)
    assert result == [None]


def test_produce_list_centroids_single_square():
    # Create an array of particle coordinates
    pts = np.array([[0, 0], [1, 1], [2, 2], [3, 3]])
    # Choose indices that pick points [1] and [3]
    indices_tuple = (np.array([1, 3]),)
    list_indices = [indices_tuple]
    result = streamlines._produce_list_centroids_this_frame(list_indices, pts)
    expected = np.array([2.0, 2.0])
    np.testing.assert_allclose(result[0], expected)


def test_produce_list_centroids_multiple_squares():
    pts = np.array([[0, 0], [2, 2], [4, 4], [6, 6]])
    # First square will use pts[0] and pts[2] -> average is [2,2]
    indices1 = (np.array([0, 2]),)
    # Second square will use pts[1] and pts[3] -> average is [4,4]
    indices2 = (np.array([1, 3]),)
    list_indices = [indices1, indices2]
    result = streamlines._produce_list_centroids_this_frame(list_indices, pts)
    expected1 = np.array([2, 2])
    expected2 = np.array([4, 4])
    np.testing.assert_array_equal(result[0], expected1)
    np.testing.assert_array_equal(result[1], expected2)


def test_adjacent_squares():
    # Test two adjacent squares that share a boundary.
    # Square1 covers [0,1]x[0,1] and square2 covers [1,2]x[0,1].
    # A point at [0.5, 0.5] should be in square1 and one at [1.5,0.5] should be in square2.
    square1 = [(0, 0), (1, 0), (1, 1), (0, 1)]
    square2 = [(1, 0), (2, 0), (2, 1), (1, 1)]
    vertex_list = [square1, square2]
    points = np.array(
        [
            [0.5, 0.5],
            [1.5, 0.5],
        ]
    )
    result = streamlines._produce_list_indices_point_in_polygon_this_frame(
        vertex_list, points
    )
    expected = [(np.array([0]),), (np.array([1]),)]
    for res_tuple, exp_tuple in zip(result, expected):
        np.testing.assert_array_equal(res_tuple[0], exp_tuple[0])


def test_point_on_boundary():
    # Test that a point exactly on the square's boundary is not considered inside.
    # For a square covering [0,1]x[0,1], a point at [1,0.5] lies exactly on the right edge.
    # By default, matplotlib.path.Path.contains_points includes boundary points.
    square = [(0, 0), (1, 0), (1, 1), (0, 1)]
    vertex_list = [square]
    points = np.array([[1, 0.5]])  # exactly on the boundary
    result = streamlines._produce_list_indices_point_in_polygon_this_frame(
        vertex_list, points
    )
    expected = [(np.array([0], dtype=int),)]
    np.testing.assert_array_equal(result, expected)


def test_points_on_boundary_of_two_adjacent_squares():
    square1 = [(0, 0), (1, 0), (1, 1), (0, 1)]
    square2 = [(1, 0), (2, 0), (2, 1), (1, 1)]
    vertex_list = [square1, square2]
    points = np.array([[1, 0.5], [1, 0.7]])  # exactly on the boundary
    result = streamlines._produce_list_indices_point_in_polygon_this_frame(
        vertex_list, points
    )
    expected = [(np.array([0, 1], dtype=int),), (np.array([0, 1], dtype=int),)]
    np.testing.assert_array_equal(result, expected)


def test_per_core_work_2D(membrane_xtc, univ):
    xmin = univ.atoms.positions[..., 0].min()
    xmax = univ.atoms.positions[..., 0].max()
    ymin = univ.atoms.positions[..., 1].min()
    ymax = univ.atoms.positions[..., 1].max()
    tuple_of_limits = (xmin, xmax, ymin, ymax)
    grid = streamlines.produce_grid(
        tuple_of_limits=tuple_of_limits, grid_spacing=20
    )
    (
        list_square_vertex_arrays_per_core,
        list_parent_index_values,
        _,
        _,
    ) = streamlines.split_grid(grid=grid, num_cores=1)
    values = streamlines.per_core_work(
        topology_file_path=Martini_membrane_gro,
        trajectory_file_path=membrane_xtc,
        list_square_vertex_arrays_this_core=list_square_vertex_arrays_per_core[
            0
        ],
        MDA_selection="name PO4",
        start_frame=1,
        end_frame=2,
        reconstruction_index_list=list_parent_index_values[0],
        maximum_delta_magnitude=2.0,
    )
    for entry in values:
        res = entry[1]
        np.testing.assert_allclose(res[:2], np.array([0.7, 0.46]), atol=1e-1)

    # When start_frame exceeds trajectory length, num_steps == 0
    zero_step_values = streamlines.per_core_work(
        topology_file_path=Martini_membrane_gro,
        trajectory_file_path=membrane_xtc,
        list_square_vertex_arrays_this_core=list_square_vertex_arrays_per_core[
            0
        ],
        MDA_selection="name PO4",
        start_frame=100,
        end_frame=105,
        reconstruction_index_list=list_parent_index_values[0],
        maximum_delta_magnitude=2.0,
    )
    for entry in zero_step_values:
        np.testing.assert_allclose(entry[1], [0.0, 0.0])


def test_streamplot_2D(membrane_xtc, univ):
    # regression test the data structures
    # generated by the 2D streamplot code
    u1, v1, avg, std = streamlines.generate_streamlines(
        topology_file_path=Martini_membrane_gro,
        trajectory_file_path=membrane_xtc,
        grid_spacing=20,
        MDA_selection="name PO4",
        start_frame=1,
        end_frame=2,
        xmin=univ.atoms.positions[..., 0].min(),
        xmax=univ.atoms.positions[..., 0].max(),
        ymin=univ.atoms.positions[..., 1].min(),
        ymax=univ.atoms.positions[..., 1].max(),
        maximum_delta_magnitude=2.0,
        num_cores=1,
    )
    assert_allclose(
        u1,
        np.array(
            [
                [0.70000172, 0.70000267, 0.69999313, 0.70000458, 0.69999695],
                [0.70000076, 0.70000267, 0.70000076, 0.69999695, 0.70000458],
                [0.70000076, 0.70000076, 0.70000458, 0.69999695, 0.70000458],
                [0.69999981, 0.69999886, 0.69999695, 0.69998932, 0.69999695],
                [0.70000076, 0.69999695, 0.69999313, 0.69999695, 0.69999695],
            ]
        ),
        atol=1e-4,
    )
    assert_allclose(
        v1,
        np.array(
            [
                [0.46000004, 0.45999908, 0.46000004, 0.46000099, 0.45999908],
                [0.45999718, 0.45999908, 0.45999908, 0.46000481, 0.46000099],
                [0.45999908, 0.45999908, 0.45999908, 0.4600029, 0.4600029],
                [0.45999908, 0.45999146, 0.46000671, 0.46001434, 0.45999908],
                [0.45999146, 0.45999908, 0.45998383, 0.45999908, 0.45999908],
            ]
        ),
        atol=1e-4,
    )
    assert avg == pytest.approx(0.837614669, abs=1e-4)
    assert std == pytest.approx(4.039906815e-06, abs=1e-5)


def test_streamplot_2D_multiframe_average(membrane_xtc, univ):
    # Regression test for Issue #4549: ensure streamlines accumulates
    # and time-averages across all intermediate frames in [start_frame, end_frame].
    # Step 0->1 delta: [0.6, 0.38]
    # Step 1->2 delta: [0.7, 0.46]
    # Step 2->3 delta: [0.8, 0.54]

    # Two-step window: frames 0 to 2
    u_2steps, v_2steps, _, _ = streamlines.generate_streamlines(
        topology_file_path=Martini_membrane_gro,
        trajectory_file_path=membrane_xtc,
        grid_spacing=20,
        MDA_selection="name PO4",
        start_frame=0,
        end_frame=2,
        xmin=univ.atoms.positions[..., 0].min(),
        xmax=univ.atoms.positions[..., 0].max(),
        ymin=univ.atoms.positions[..., 1].min(),
        ymax=univ.atoms.positions[..., 1].max(),
        maximum_delta_magnitude=2.0,
        num_cores=1,
    )
    assert_allclose(u_2steps, 0.65, atol=1e-3)
    assert_allclose(v_2steps, 0.42, atol=1e-3)

    # Three-step window: frames 0 to 3
    u_3steps, v_3steps, _, _ = streamlines.generate_streamlines(
        topology_file_path=Martini_membrane_gro,
        trajectory_file_path=membrane_xtc,
        grid_spacing=20,
        MDA_selection="name PO4",
        start_frame=0,
        end_frame=3,
        xmin=univ.atoms.positions[..., 0].min(),
        xmax=univ.atoms.positions[..., 0].max(),
        ymin=univ.atoms.positions[..., 1].min(),
        ymax=univ.atoms.positions[..., 1].max(),
        maximum_delta_magnitude=2.0,
        num_cores=1,
    )
    assert_allclose(u_3steps, 0.70, atol=1e-3)
    assert_allclose(v_3steps, 0.46, atol=1e-3)


@pytest.mark.parametrize(
    "kwargs, err_type, err_match",
    [
        (
            {"start_frame": "0"},
            TypeError,
            "start_frame and end_frame must be integers",
        ),
        (
            {"end_frame": 1.5},
            TypeError,
            "start_frame and end_frame must be integers",
        ),
        (
            {"start_frame": -1},
            ValueError,
            "start_frame must be non-negative",
        ),
        (
            {"start_frame": 2, "end_frame": 2},
            ValueError,
            "end_frame must be greater than start_frame",
        ),
        (
            {"start_frame": 3, "end_frame": 1},
            ValueError,
            "end_frame must be greater than start_frame",
        ),
        (
            {"grid_spacing": 0},
            ValueError,
            "grid_spacing must be positive",
        ),
        (
            {"grid_spacing": -5},
            ValueError,
            "grid_spacing must be positive",
        ),
        (
            {"maximum_delta_magnitude": 0},
            ValueError,
            "maximum_delta_magnitude must be positive",
        ),
        (
            {"maximum_delta_magnitude": -1.0},
            ValueError,
            "maximum_delta_magnitude must be positive",
        ),
        (
            {"xmin": 10, "xmax": 5},
            ValueError,
            "xmin must be less than xmax",
        ),
        (
            {"ymin": 20, "ymax": 20},
            ValueError,
            "ymin must be less than ymax",
        ),
        (
            {"num_cores": 0},
            ValueError,
            "num_cores must be 'maximum' or a positive integer",
        ),
        (
            {"num_cores": -1},
            ValueError,
            "num_cores must be 'maximum' or a positive integer",
        ),
        (
            {"num_cores": "all"},
            ValueError,
            "num_cores must be 'maximum' or a positive integer",
        ),
    ],
)
def test_generate_streamlines_invalid_inputs(
    membrane_xtc, univ, kwargs, err_type, err_match
):
    default_kwargs = {
        "topology_file_path": Martini_membrane_gro,
        "trajectory_file_path": membrane_xtc,
        "grid_spacing": 20,
        "MDA_selection": "name PO4",
        "start_frame": 0,
        "end_frame": 1,
        "xmin": univ.atoms.positions[..., 0].min(),
        "xmax": univ.atoms.positions[..., 0].max(),
        "ymin": univ.atoms.positions[..., 1].min(),
        "ymax": univ.atoms.positions[..., 1].max(),
        "maximum_delta_magnitude": 2.0,
        "num_cores": 1,
    }
    default_kwargs.update(kwargs)
    with pytest.raises(err_type, match=err_match):
        streamlines.generate_streamlines(**default_kwargs)


def test_streamplot_2D_zero_return(membrane_xtc, univ, tmpdir):
    # simple roundtrip test to ensure that
    # zeroed arrays are returned by the 2D streamplot
    # code when called with an empty selection
    u1, v1, avg, std = streamlines.generate_streamlines(
        topology_file_path=Martini_membrane_gro,
        trajectory_file_path=membrane_xtc,
        grid_spacing=20,
        MDA_selection="name POX",
        start_frame=1,
        end_frame=2,
        xmin=univ.atoms.positions[..., 0].min(),
        xmax=univ.atoms.positions[..., 0].max(),
        ymin=univ.atoms.positions[..., 1].min(),
        ymax=univ.atoms.positions[..., 1].max(),
        maximum_delta_magnitude=2.0,
        num_cores=1,
    )
    assert_allclose(u1, np.zeros((5, 5)))
    assert_allclose(v1, np.zeros((5, 5)))
    assert avg == approx(0.0)
    assert std == approx(0.0)


def test_streamplot_2D_dual_core(membrane_xtc, univ, tmpdir):
    # simple test to ensure that it runs with multiple cores
    u1, v1, avg, std = streamlines.generate_streamlines(
        topology_file_path=Martini_membrane_gro,
        trajectory_file_path=membrane_xtc,
        grid_spacing=20,
        MDA_selection="name POX",
        start_frame=1,
        end_frame=2,
        xmin=univ.atoms.positions[..., 0].min(),
        xmax=univ.atoms.positions[..., 0].max(),
        ymin=univ.atoms.positions[..., 1].min(),
        ymax=univ.atoms.positions[..., 1].max(),
        maximum_delta_magnitude=2.0,
        num_cores=2,
    )
    assert_allclose(u1, np.zeros((5, 5)))
    assert_allclose(v1, np.zeros((5, 5)))
    assert avg == approx(0.0)
    assert std == approx(0.0)


def test_streamplot_2D_default_cores(membrane_xtc, univ):
    # Test using default num_cores="maximum"
    u1, v1, avg, std = streamlines.generate_streamlines(
        topology_file_path=Martini_membrane_gro,
        trajectory_file_path=membrane_xtc,
        grid_spacing=20,
        MDA_selection="name POX",
        start_frame=1,
        end_frame=2,
        xmin=univ.atoms.positions[..., 0].min(),
        xmax=univ.atoms.positions[..., 0].max(),
        ymin=univ.atoms.positions[..., 1].min(),
        ymax=univ.atoms.positions[..., 1].max(),
        maximum_delta_magnitude=2.0,
    )
    assert_allclose(u1, np.zeros((5, 5)))
    assert_allclose(v1, np.zeros((5, 5)))


def test_streamplot_2D_frames_beyond_trajectory(membrane_xtc, univ):
    # Test when start_frame is beyond trajectory length (0 steps evaluated)
    u_past, v_past, _, _ = streamlines.generate_streamlines(
        topology_file_path=Martini_membrane_gro,
        trajectory_file_path=membrane_xtc,
        grid_spacing=20,
        MDA_selection="name PO4",
        start_frame=100,
        end_frame=105,
        xmin=univ.atoms.positions[..., 0].min(),
        xmax=univ.atoms.positions[..., 0].max(),
        ymin=univ.atoms.positions[..., 1].min(),
        ymax=univ.atoms.positions[..., 1].max(),
        maximum_delta_magnitude=2.0,
        num_cores=1,
    )
    assert_allclose(u_past, np.zeros((5, 5)))
    assert_allclose(v_past, np.zeros((5, 5)))

    # Test when end_frame exceeds trajectory length (loop completes naturally)
    u_end, v_end, _, _ = streamlines.generate_streamlines(
        topology_file_path=Martini_membrane_gro,
        trajectory_file_path=membrane_xtc,
        grid_spacing=20,
        MDA_selection="name PO4",
        start_frame=0,
        end_frame=100,
        xmin=univ.atoms.positions[..., 0].min(),
        xmax=univ.atoms.positions[..., 0].max(),
        ymin=univ.atoms.positions[..., 1].min(),
        ymax=univ.atoms.positions[..., 1].max(),
        maximum_delta_magnitude=2.0,
        num_cores=1,
    )
    assert np.all(u_end > 0)
    assert np.all(v_end > 0)


def test_streamplot_3D(membrane_xtc, univ, tmpdir):
    # because mayavi is too heavy of a dependency
    # for a roundtrip plotting test, simply
    # aim to check for sensible values
    # returned by generate_streamlines_3d
    dx, dy, dz = streamlines_3D.generate_streamlines_3d(
        topology_file_path=Martini_membrane_gro,
        trajectory_file_path=membrane_xtc,
        grid_spacing=20,
        MDA_selection="name PO4",
        start_frame=1,
        end_frame=2,
        xmin=univ.atoms.positions[..., 0].min(),
        xmax=univ.atoms.positions[..., 0].max(),
        ymin=univ.atoms.positions[..., 1].min(),
        ymax=univ.atoms.positions[..., 1].max(),
        zmin=univ.atoms.positions[..., 2].min(),
        zmax=univ.atoms.positions[..., 2].max(),
        maximum_delta_magnitude=2.0,
        num_cores=1,
    )
    assert dx.shape == (5, 5, 2)
    assert dy.shape == (5, 5, 2)
    assert dz.shape == (5, 5, 2)
    assert dx[4, 4, 0] == approx(0.700004, abs=1e-5)
    assert dy[0, 0, 0] == approx(0.460000, abs=1e-5)
    assert dz[2, 2, 0] == approx(0.240005, abs=1e-5)
