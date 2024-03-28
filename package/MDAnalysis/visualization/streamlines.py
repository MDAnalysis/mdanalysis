# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
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

"""Streamplots (2D)  --- :mod:`MDAnalysis.visualization.streamlines`
=================================================================

:Authors: Tyler Reddy and Matthieu Chavent
:Year: 2014
:Copyright: GNU Public License v3


The :func:`generate_streamlines` function can generate a 2D flow field from a
MD trajectory, for instance, lipid molecules in a flat membrane. It can make
use of multiple cores to perform the analyis in parallel (using
:mod:`multiprocessing`).

See Also
--------
MDAnalysis.visualization.streamlines_3D : streamplots in 3D


.. autofunction:: generate_streamlines

"""
import multiprocessing

import numpy as np

try:
    import matplotlib
    import matplotlib.path
except ImportError:
    raise ImportError(
            '2d streamplot module requires: matplotlib.path for its '
            'path.Path.contains_points method. The installation '
            'instructions for the matplotlib module can be found here: '
            'http://matplotlib.org/faq/installing_faq.html?highlight=install'
            ) from None

import MDAnalysis



def produce_grid(tuple_of_limits, grid_spacing):
    """Produce a 2D grid for the simulation system.

    The grid is based on the tuple of Cartesian Coordinate limits calculated in
    an earlier step.

    Parameters
    ----------
    tuple_of_limits : tuple
        ``x_min, x_max, y_min, y_max``
    grid_spacing : float
        grid size in all directions in ångström

    Returns
    -------
    grid : array
       ``numpy.mgrid[x_min:x_max:grid_spacing, y_min:y_max:grid_spacing]``
    """
    x_min, x_max, y_min, y_max = tuple_of_limits
    grid = np.mgrid[x_min:x_max:grid_spacing, y_min:y_max:grid_spacing]
    return grid


def split_grid(grid, num_cores):
    """Split the grid into blocks of vertices.

    Take the overall `grid` for the system and split it into lists of
    square vertices that can be distributed to each core.

    Parameters
    ----------
    grid : numpy.array
        2D array
    num_cores : int
        number of partitions to generate

    Returns
    -------
    list_square_vertex_arrays_per_core : array of arrays
         split the list of square vertices
         ``[[v1,v2,v3,v4],[v1,v2,v3,v4],...,...]`` into roughly equally-sized
         sublists to be distributed over the available cores on the system
    list_parent_index_values : array of arrays
         arrays of `[[row, column], [row, column], ...]`` for each core
    current_row : int
         last row + 1
    current_column : int
         last column + 1

    Note
    ----
    Limited to 2D for now.

    """

    # produce an array containing the cartesian coordinates of all vertices in the grid:
    x_array, y_array = grid
    grid_vertex_cartesian_array = np.dstack((x_array, y_array))
    #the grid_vertex_cartesian_array has N_rows, with each row corresponding to a column of coordinates in the grid (
    # so a given row has shape N_rows, 2); overall shape (N_columns_in_grid, N_rows_in_a_column, 2)
    #although I'll eventually want a pure numpy/scipy/vector-based solution, for now I'll allow loops to simplify the
    #  division of the cartesian coordinates into a list of the squares in the grid
    list_all_squares_in_grid = []  # should eventually be a nested list of all the square vertices in the grid/system
    list_parent_index_values = []  # want an ordered list of assignment indices for reconstructing the grid positions
    # in the parent process
    current_column = 0
    while current_column < grid_vertex_cartesian_array.shape[0] - 1:
        # go through all the columns except the last one and account for the square vertices (the last column
        #  has no 'right neighbour')
        current_row = 0
        while current_row < grid_vertex_cartesian_array.shape[1] - 1:
            # all rows except the top row, which doesn't have a row above it for forming squares
            bottom_left_vertex_current_square = grid_vertex_cartesian_array[current_column, current_row]
            bottom_right_vertex_current_square = grid_vertex_cartesian_array[current_column + 1, current_row]
            top_right_vertex_current_square = grid_vertex_cartesian_array[current_column + 1, current_row + 1]
            top_left_vertex_current_square = grid_vertex_cartesian_array[current_column, current_row + 1]
            #append the vertices of this square to the overall list of square vertices:
            list_all_squares_in_grid.append(
                [bottom_left_vertex_current_square, bottom_right_vertex_current_square, top_right_vertex_current_square,
                    top_left_vertex_current_square])
            list_parent_index_values.append([current_row, current_column])
            current_row += 1
        current_column += 1
    #split the list of square vertices [[v1,v2,v3,v4],[v1,v2,v3,v4],...,...] into roughly equally-sized sublists to
    # be distributed over the available cores on the system:
    list_square_vertex_arrays_per_core = np.array_split(list_all_squares_in_grid, num_cores)
    list_parent_index_values = np.array_split(list_parent_index_values, num_cores)
    return [list_square_vertex_arrays_per_core, list_parent_index_values, current_row, current_column]


def per_core_work(topology_file_path, trajectory_file_path, list_square_vertex_arrays_this_core, MDA_selection,
                  start_frame, end_frame, reconstruction_index_list, maximum_delta_magnitude):
    """Run the analysis on one core.

    The code to perform on a given core given the list of square vertices assigned to it.
    """
    # obtain the relevant coordinates for particles of interest
    universe_object = MDAnalysis.Universe(topology_file_path, trajectory_file_path)
    list_previous_frame_centroids = []
    list_previous_frame_indices = []
    #define some utility functions for trajectory iteration:

    def produce_list_indices_point_in_polygon_this_frame(vertex_coord_list):
        list_indices_point_in_polygon = []
        for square_vertices in vertex_coord_list:
            path_object = matplotlib.path.Path(square_vertices)
            index_list_in_polygon = np.where(path_object.contains_points(relevant_particle_coordinate_array_xy))
            list_indices_point_in_polygon.append(index_list_in_polygon)
        return list_indices_point_in_polygon

    def produce_list_centroids_this_frame(list_indices_in_polygon):
        list_centroids_this_frame = []
        for indices in list_indices_in_polygon:
            if not indices[0].size > 0:  # if there are no particles of interest in this particular square
                list_centroids_this_frame.append(None)
            else:
                current_coordinate_array_in_square = relevant_particle_coordinate_array_xy[indices]
                current_square_indices_centroid = np.average(current_coordinate_array_in_square, axis=0)
                list_centroids_this_frame.append(current_square_indices_centroid)
        return list_centroids_this_frame  # a list of numpy xy centroid arrays for this frame

    for ts in universe_object.trajectory:
        if ts.frame < start_frame:  # don't start until first specified frame
            continue
        relevant_particle_coordinate_array_xy = universe_object.select_atoms(MDA_selection).positions[..., :-1]
        # only 2D / xy coords for now
        #I will need a list of indices for relevant particles falling within each square in THIS frame:
        list_indices_in_squares_this_frame = produce_list_indices_point_in_polygon_this_frame(
            list_square_vertex_arrays_this_core)
        #likewise, I will need a list of centroids of particles in each square (same order as above list):
        list_centroids_in_squares_this_frame = produce_list_centroids_this_frame(list_indices_in_squares_this_frame)
        if list_previous_frame_indices:  # if the previous frame had indices in at least one square I will need to use
            #  those indices to generate the updates to the corresponding centroids in this frame:
            list_centroids_this_frame_using_indices_from_last_frame = produce_list_centroids_this_frame(
                list_previous_frame_indices)
            #I need to write a velocity of zero if there are any 'empty' squares in either frame:
            xy_deltas_to_write = []
            for square_1_centroid, square_2_centroid in zip(list_centroids_this_frame_using_indices_from_last_frame,
                                                            list_previous_frame_centroids):
                if square_1_centroid is None or square_2_centroid is None:
                    xy_deltas_to_write.append([0, 0])
                else:
                    xy_deltas_to_write.append(np.subtract(square_1_centroid, square_2_centroid).tolist())

            #xy_deltas_to_write = np.subtract(np.array(
            # list_centroids_this_frame_using_indices_from_last_frame),np.array(list_previous_frame_centroids))
            xy_deltas_to_write = np.array(xy_deltas_to_write)
            #now filter the array to only contain distances in the range [-8,8] as a placeholder for dealing with PBC
            #  issues (Matthieu seemed to use a limit of 8 as well);
            xy_deltas_to_write = np.clip(xy_deltas_to_write, -maximum_delta_magnitude, maximum_delta_magnitude)

            #with the xy and dx,dy values calculated I need to set the values from this frame to previous frame
            # values in anticipation of the next frame:
            list_previous_frame_centroids = list_centroids_in_squares_this_frame[:]
            list_previous_frame_indices = list_indices_in_squares_this_frame[:]
        else:  # either no points in squares or after the first frame I'll just reset the 'previous' values so they
            # can be used when consecutive frames have proper values
            list_previous_frame_centroids = list_centroids_in_squares_this_frame[:]
            list_previous_frame_indices = list_indices_in_squares_this_frame[:]
        if ts.frame > end_frame:
            break  # stop here
    return list(zip(reconstruction_index_list, xy_deltas_to_write.tolist()))


def generate_streamlines(topology_file_path, trajectory_file_path, grid_spacing, MDA_selection, start_frame,
                         end_frame, xmin, xmax, ymin, ymax, maximum_delta_magnitude, num_cores='maximum'):
    r"""Produce the x and y components of a 2D streamplot data set.

    Parameters
    ----------
    topology_file_path : str
            Absolute path to the topology file
    trajectory_file_path : str
            Absolute path to the trajectory file. It will normally be desirable
            to filter the trajectory with a tool such as GROMACS
            :program:`g_filter` (see :footcite:p:`Chavent2014`)
    grid_spacing : float
            The spacing between grid lines (angstroms)
    MDA_selection : str
            MDAnalysis selection string
    start_frame : int
            First frame number to parse
    end_frame : int
            Last frame number to parse
    xmin : float
            Minimum coordinate boundary for x-axis (angstroms)
    xmax : float
            Maximum coordinate boundary for x-axis (angstroms)
    ymin : float
            Minimum coordinate boundary for y-axis (angstroms)
    ymax : float
            Maximum coordinate boundary for y-axis (angstroms)
    maximum_delta_magnitude : float
            Absolute value of the largest displacement tolerated for the
            centroid of a group of particles ( angstroms). Values above this
            displacement will not count in the streamplot (treated as
            excessively large displacements crossing the periodic boundary)
    num_cores : int or 'maximum' (optional)
            The number of cores to use. (Default 'maximum' uses all available
            cores)

    Returns
    -------
    dx_array : array of floats
            An array object containing the displacements in the x direction
    dy_array : array of floats
            An array object containing the displacements in the y direction
    average_displacement : float
            :math:`\frac{\sum\sqrt[]{dx^2 + dy^2}}{N}`
    standard_deviation_of_displacement : float
            standard deviation of :math:`\sqrt[]{dx^2 + dy^2}`

    Examples
    --------
    Generate 2D streamlines and plot::

        import matplotlib, matplotlib.pyplot, np
        import MDAnalysis, MDAnalysis.visualization.streamlines

        u1, v1, average_displacement, standard_deviation_of_displacement =
            MDAnalysis.visualization.streamlines.generate_streamlines('testing.gro', 'testing_filtered.xtc',
                    grid_spacing=20, MDA_selection='name PO4', start_frame=2, end_frame=3,
                    xmin=-8.73000049591, xmax= 1225.96008301,
                    ymin= -12.5799999237, ymax=1224.34008789,
                    maximum_delta_magnitude=1.0, num_cores=16)
        x = np.linspace(0, 1200, 61)
        y = np.linspace(0, 1200, 61)
        speed = np.sqrt(u1*u1 + v1*v1)
        fig = matplotlib.pyplot.figure()
        ax = fig.add_subplot(111, aspect='equal')
        ax.set_xlabel('x ($\AA$)')
        ax.set_ylabel('y ($\AA$)')
        ax.streamplot(x, y, u1, v1, density=(10,10), color=speed, linewidth=3*speed/speed.max())
        fig.savefig('testing_streamline.png',dpi=300)

    .. image:: testing_streamline.png


    References
    ----------
    .. footbibliography::


    See Also
    --------
    MDAnalysis.visualization.streamlines_3D.generate_streamlines_3d

    """
    # work out the number of cores to use:
    if num_cores == 'maximum':
        num_cores = multiprocessing.cpu_count()  # use all available cores
    else:
        num_cores = num_cores  # use the value specified by the user
        #assert isinstance(num_cores,(int,long)), "The number of specified cores must (of course) be an integer."
    np.seterr(all='warn', over='raise')
    parent_list_deltas = []  # collect all data from child processes here

    def log_result_to_parent(delta_array):
        parent_list_deltas.extend(delta_array)

    tuple_of_limits = (xmin, xmax, ymin, ymax)
    grid = produce_grid(tuple_of_limits=tuple_of_limits, grid_spacing=grid_spacing)
    list_square_vertex_arrays_per_core, list_parent_index_values, total_rows, total_columns = \
        split_grid(grid=grid,
                   num_cores=num_cores)
    pool = multiprocessing.Pool(num_cores)
    for vertex_sublist, index_sublist in zip(list_square_vertex_arrays_per_core, list_parent_index_values):
        pool.apply_async(per_core_work, args=(
            topology_file_path, trajectory_file_path, vertex_sublist, MDA_selection, start_frame, end_frame,
            index_sublist, maximum_delta_magnitude), callback=log_result_to_parent)
    pool.close()
    pool.join()
    dx_array = np.zeros((total_rows, total_columns))
    dy_array = np.zeros((total_rows, total_columns))
    #the parent_list_deltas is shaped like this: [ ([row_index,column_index],[dx,dy]), ... (...),...,]
    for index_array, delta_array in parent_list_deltas:  # go through the list in the parent process and assign to the
        #  appropriate positions in the dx and dy matrices:
        #build in a filter to replace all values at the cap (currently between -8,8) with 0 to match Matthieu's code
        # (I think eventually we'll reduce the cap to a narrower boundary though)
        index_1 = index_array.tolist()[0]
        index_2 = index_array.tolist()[1]
        if abs(delta_array[0]) == maximum_delta_magnitude:
            dx_array[index_1, index_2] = 0
        else:
            dx_array[index_1, index_2] = delta_array[0]
        if abs(delta_array[1]) == maximum_delta_magnitude:
            dy_array[index_1, index_2] = 0
        else:
            dy_array[index_1, index_2] = delta_array[1]

    #at Matthieu's request, we now want to calculate the average and standard deviation of the displacement values:
    displacement_array = np.sqrt(dx_array ** 2 + dy_array ** 2)
    average_displacement = np.average(displacement_array)
    standard_deviation_of_displacement = np.std(displacement_array)

    return (dx_array, dy_array, average_displacement, standard_deviation_of_displacement)
