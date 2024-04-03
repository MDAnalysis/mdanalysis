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

"""Streamplots (3D) --- :mod:`MDAnalysis.visualization.streamlines_3D`
===================================================================

:Authors: Tyler Reddy and Matthieu Chavent
:Year: 2014
:Copyright: GNU Public License v3


The :func:`generate_streamlines_3d` function can generate a 3D flow field from
a MD trajectory, for instance, lipid molecules in a virus capsid. It can make
use of multiple cores to perform the analyis in parallel (using
:mod:`multiprocessing`).

.. rubric: References

.. footbibliography::

See Also
--------
MDAnalysis.visualization.streamlines : streamplots in 2D


.. autofunction:: generate_streamlines_3d

"""
import multiprocessing

import numpy as np
import numpy.testing
import scipy
import scipy.spatial.distance

import MDAnalysis


def determine_container_limits(topology_file_path, trajectory_file_path,
                               buffer_value):
    """Calculate the extent of the atom coordinates + buffer.

    A function for the parent process which should take the input trajectory
    and calculate the limits of the container for the system and return these
    limits.

    Parameters
    ----------
    topology_file_path : str
        topology file name
    trajectory_file_path : str
        trajectory file name
    buffer_value : float
        buffer value (padding) in +/- {x, y, z}
    """
    universe_object = MDAnalysis.Universe(topology_file_path, trajectory_file_path)
    all_atom_selection = universe_object.select_atoms('all')  # select all particles
    all_atom_coordinate_array = all_atom_selection.positions
    x_min, x_max, y_min, y_max, z_min, z_max = [
        all_atom_coordinate_array[..., 0].min(),
        all_atom_coordinate_array[..., 0].max(), all_atom_coordinate_array[..., 1].min(),
        all_atom_coordinate_array[..., 1].max(), all_atom_coordinate_array[..., 2].min(),
        all_atom_coordinate_array[..., 2].max()]
    tuple_of_limits = \
        (
            x_min - buffer_value,
            x_max + buffer_value, y_min - buffer_value, y_max + buffer_value, z_min - buffer_value,
            z_max + buffer_value)  # using buffer_value to catch particles near edges
    return tuple_of_limits


def produce_grid(tuple_of_limits, grid_spacing):
    """Produce a 3D grid for the simulation system.

    The partitioning is based on the tuple of Cartesian Coordinate limits
    calculated in an earlier step.

    Parameters
    ----------
    tuple_of_limits : tuple
        ``x_min, x_max, y_min, y_max, z_min, z_max``
    grid_spacing : float
        grid size in all directions in ångström

    Returns
    -------
    grid : array
        ``numpy.mgrid[x_min:x_max:grid_spacing, y_min:y_max:grid_spacing, z_min:z_max:grid_spacing]``

    """
    x_min, x_max, y_min, y_max, z_min, z_max = tuple_of_limits
    grid = np.mgrid[x_min:x_max:grid_spacing, y_min:y_max:grid_spacing, z_min:z_max:grid_spacing]
    return grid


def split_grid(grid, num_cores):
    """Split the grid into blocks of vertices.

    Take the overall `grid` for the system and split it into lists of cube
    vertices that can be distributed to each core.

    Parameters
    ----------
    grid : numpy.array
        3D array
    num_cores : int
        number of partitions to generate

    Returns
    -------
    list_dictionaries_for_cores : list of dict
    total_cubes : int
    num_sheets : int
    delta_array_shape : tuple

    """
    # unpack the x,y,z mgrid arrays
    x, y, z = grid
    num_z_values = z.shape[-1]
    num_sheets = z.shape[0]
    delta_array_shape = tuple(
        [n - 1 for n in x.shape])  # the final target shape for return delta arrays is n-1 in each dimension

    ordered_list_per_sheet_x_values = []
    for x_sheet in x:  # each x_sheet should have shape (25,23) and the same x value in each element
        array_all_x_values_current_sheet = x_sheet.flatten()
        ordered_list_per_sheet_x_values.append(array_all_x_values_current_sheet)
    ordered_list_per_sheet_y_values = []
    for y_columns in y:
        array_all_y_values_current_sheet = y_columns.flatten()
        ordered_list_per_sheet_y_values.append(array_all_y_values_current_sheet)
    ordered_list_per_sheet_z_values = []
    for z_slices in z:
        array_all_z_values_current_sheet = z_slices.flatten()
        ordered_list_per_sheet_z_values.append(array_all_z_values_current_sheet)

    ordered_list_cartesian_coordinates_per_sheet = []
    for x_sheet_coords, y_sheet_coords, z_sheet_coords in zip(ordered_list_per_sheet_x_values,
                                                              ordered_list_per_sheet_y_values,
                                                              ordered_list_per_sheet_z_values):
        ordered_list_cartesian_coordinates_per_sheet.append(list(zip(x_sheet_coords, y_sheet_coords, z_sheet_coords)))
    array_ordered_cartesian_coords_per_sheet = np.array(ordered_list_cartesian_coordinates_per_sheet)
    #now I'm going to want to build cubes in an ordered fashion, and in such a way that I can track the index /
    # centroid of each cube for domain decomposition / reconstruction and mayavi mlab.flow() input
    #cubes will be formed from N - 1 base sheets combined with subsequent sheets
    current_base_sheet = 0
    dictionary_cubes_centroids_indices = {}
    cube_counter = 0
    while current_base_sheet < num_sheets - 1:
        current_base_sheet_array = array_ordered_cartesian_coords_per_sheet[current_base_sheet]
        current_top_sheet_array = array_ordered_cartesian_coords_per_sheet[
            current_base_sheet + 1]  # the points of the sheet 'to the right' in the grid
        current_index = 0
        while current_index < current_base_sheet_array.shape[0] - num_z_values:
            # iterate through all the indices in each of the sheet arrays (careful to avoid extra
            # points not needed for cubes)
            column_z_level = 0  # start at the bottom of a given 4-point column and work up
            while column_z_level < num_z_values - 1:
                current_list_cube_vertices = []
                first_two_vertices_base_sheet = current_base_sheet_array[current_index:current_index + 2, ...].tolist()
                first_two_vertices_top_sheet = current_top_sheet_array[current_index:current_index + 2, ...].tolist()
                next_two_vertices_base_sheet = current_base_sheet_array[current_index +
                                                                        num_z_values: 2 +
                                                                        num_z_values + current_index, ...].tolist()
                next_two_vertices_top_sheet = current_top_sheet_array[current_index +
                                                                      num_z_values: 2 +
                                                                      num_z_values + current_index, ...].tolist()
                for vertex_set in [
                    first_two_vertices_base_sheet, first_two_vertices_top_sheet,
                    next_two_vertices_base_sheet, next_two_vertices_top_sheet
                ]:
                    current_list_cube_vertices.extend(vertex_set)
                vertex_array = np.array(current_list_cube_vertices)
                assert vertex_array.shape == (8, 3), "vertex_array has incorrect shape"
                cube_centroid = np.average(np.array(current_list_cube_vertices), axis=0)
                dictionary_cubes_centroids_indices[cube_counter] = {
                    'centroid': cube_centroid,
                    'vertex_list': current_list_cube_vertices}
                cube_counter += 1
                current_index += 1
                column_z_level += 1
                if column_z_level == num_z_values - 1:  # the loop will break but I should also increment the
                    # current_index
                    current_index += 1
        current_base_sheet += 1
    total_cubes = len(dictionary_cubes_centroids_indices)

    #produce an array of pseudo cube indices (actually the dictionary keys which are cube numbers in string format):
    pseudo_cube_indices = np.arange(0, total_cubes)
    sublist_of_cube_indices_per_core = np.array_split(pseudo_cube_indices, num_cores)
    #now, the split of pseudoindices seems to work well, and the above sublist_of_cube_indices_per_core is a list of
    # arrays of cube numbers / keys in the original dictionary
    #now I think I'll try to produce a list of dictionaries that each contain their assigned cubes based on the above
    #  per core split
    list_dictionaries_for_cores = []
    subdictionary_counter = 0
    for array_cube_indices in sublist_of_cube_indices_per_core:
        current_core_dictionary = {}
        items_to_pop = len(array_cube_indices)
        items_popped = 0
        while items_popped < items_to_pop:
            key, value = dictionary_cubes_centroids_indices.popitem()
            current_core_dictionary.update({key: value})
            items_popped += 1
        list_dictionaries_for_cores.append(current_core_dictionary)
        subdictionary_counter += 1
    return list_dictionaries_for_cores, total_cubes, num_sheets, delta_array_shape


def per_core_work(start_frame_coord_array, end_frame_coord_array, dictionary_cube_data_this_core, MDA_selection,
                  start_frame, end_frame):
    """Run the analysis on one core.

    The code to perform on a given core given the dictionary of cube data.
    """
    list_previous_frame_centroids = []
    list_previous_frame_indices = []
    # define some utility functions for trajectory iteration:

    def point_in_cube(array_point_coordinates, list_cube_vertices, cube_centroid):
        """Determine if an array of coordinates are within a cube."""
        #the simulation particle point can't be more than half the cube side length away from the cube centroid in
        # any given dimension:
        array_cube_vertices = np.array(list_cube_vertices)
        cube_half_side_length = scipy.spatial.distance.pdist(array_cube_vertices, 'euclidean').min() / 2.0
        array_cube_vertex_distances_from_centroid = scipy.spatial.distance.cdist(array_cube_vertices,
                                                                                 cube_centroid[np.newaxis, :])
        np.testing.assert_allclose(array_cube_vertex_distances_from_centroid.min(),
                                          array_cube_vertex_distances_from_centroid.max(), rtol=0, atol=1.5e-4,
                                          err_msg="not all cube vertex to centroid distances are the same, "
                                                  "so not a true cube")
        absolute_delta_coords = np.absolute(np.subtract(array_point_coordinates, cube_centroid))
        absolute_delta_x_coords = absolute_delta_coords[..., 0]
        indices_delta_x_acceptable = np.where(absolute_delta_x_coords <= cube_half_side_length)
        absolute_delta_y_coords = absolute_delta_coords[..., 1]
        indices_delta_y_acceptable = np.where(absolute_delta_y_coords <= cube_half_side_length)
        absolute_delta_z_coords = absolute_delta_coords[..., 2]
        indices_delta_z_acceptable = np.where(absolute_delta_z_coords <= cube_half_side_length)
        intersection_xy_acceptable_arrays = np.intersect1d(indices_delta_x_acceptable[0],
                                                              indices_delta_y_acceptable[0])
        overall_indices_points_in_current_cube = np.intersect1d(intersection_xy_acceptable_arrays,
                                                                   indices_delta_z_acceptable[0])
        return overall_indices_points_in_current_cube

    def update_dictionary_point_in_cube_start_frame(array_simulation_particle_coordinates,
                                                    dictionary_cube_data_this_core):
        """Basically update the cube dictionary objects assigned to this core to contain a new key/value pair
        corresponding to the indices of the relevant particles that fall within a given cube. Also, for a given cube,
        store a key/value pair for the centroid of the particles that fall within the cube."""
        cube_counter = 0
        for key, cube in dictionary_cube_data_this_core.items():
            index_list_in_cube = point_in_cube(array_simulation_particle_coordinates, cube['vertex_list'],
                                               cube['centroid'])
            cube['start_frame_index_list_in_cube'] = index_list_in_cube
            if len(index_list_in_cube) > 0:  # if there's at least one particle in this cube
                centroid_particles_in_cube = np.average(array_simulation_particle_coordinates[index_list_in_cube],
                                                           axis=0)
                cube['centroid_of_particles_first_frame'] = centroid_particles_in_cube
            else:  # empty cube
                cube['centroid_of_particles_first_frame'] = None
            cube_counter += 1

    def update_dictionary_end_frame(array_simulation_particle_coordinates, dictionary_cube_data_this_core):
        """Update the cube dictionary objects again as appropriate for the second and final frame."""
        cube_counter = 0
        for key, cube in dictionary_cube_data_this_core.items():
            # if there were no particles in the cube in the first frame, then set dx,dy,dz each to 0
            if cube['centroid_of_particles_first_frame'] is None:
                cube['dx'] = 0
                cube['dy'] = 0
                cube['dz'] = 0
            else:  # there was at least one particle in the starting cube so we can get dx,dy,dz centroid values
                new_coordinate_array_for_particles_starting_in_this_cube = array_simulation_particle_coordinates[
                    cube['start_frame_index_list_in_cube']]
                new_centroid_for_particles_starting_in_this_cube = np.average(
                    new_coordinate_array_for_particles_starting_in_this_cube, axis=0)
                cube['centroid_of_paticles_final_frame'] = new_centroid_for_particles_starting_in_this_cube
                delta_centroid_array_this_cube = new_centroid_for_particles_starting_in_this_cube - cube[
                    'centroid_of_particles_first_frame']
                cube['dx'] = delta_centroid_array_this_cube[0]
                cube['dy'] = delta_centroid_array_this_cube[1]
                cube['dz'] = delta_centroid_array_this_cube[2]
            cube_counter += 1

    #now that the parent process is dealing with the universe object & grabbing required coordinates, each child
    # process only needs to take the coordinate arrays & perform the operations with its assigned cubes (no more file
    #  opening and trajectory iteration on each core--which I'm hoping will substantially reduce the physical memory
    # footprint of my 3D streamplot code)
    update_dictionary_point_in_cube_start_frame(start_frame_coord_array, dictionary_cube_data_this_core)
    update_dictionary_end_frame(end_frame_coord_array, dictionary_cube_data_this_core)
    return dictionary_cube_data_this_core


def produce_coordinate_arrays_single_process(topology_file_path, trajectory_file_path, MDA_selection, start_frame,
                                             end_frame):
    """Generate coordinate arrays.

    To reduce memory footprint produce only a single MDA selection and get
    desired coordinate arrays; can later send these coordinate arrays to all
    child processes rather than having each child process open a trajectory and
    waste memory.

    """
    universe_object = MDAnalysis.Universe(topology_file_path, trajectory_file_path)
    relevant_particles = universe_object.select_atoms(MDA_selection)
    # pull out coordinate arrays from desired frames:
    for ts in universe_object.trajectory:
        if ts.frame > end_frame:
            break  # stop here
        if ts.frame == start_frame:
            start_frame_relevant_particle_coordinate_array_xyz = relevant_particles.positions
        elif ts.frame == end_frame:
            end_frame_relevant_particle_coordinate_array_xyz = relevant_particles.positions
        else:
            continue
    return (start_frame_relevant_particle_coordinate_array_xyz, end_frame_relevant_particle_coordinate_array_xyz)


def generate_streamlines_3d(topology_file_path, trajectory_file_path, grid_spacing, MDA_selection, start_frame,
                            end_frame, xmin, xmax, ymin, ymax, zmin, zmax, maximum_delta_magnitude=2.0,
                            num_cores='maximum'):
    r"""Produce the x, y and z components of a 3D streamplot data set.

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
    dz_array : array of floats
            An array object containing the displacements in the z direction

    Examples
    --------
    Generate 3D streamlines and visualize in `mayavi`_::

        import numpy as np

        import MDAnalysis
        import MDAnalysis.visualization.streamlines_3D

        import mayavi, mayavi.mlab

        # assign coordinate system limits and grid spacing:
        x_lower,x_upper = -8.73, 1225.96
        y_lower,y_upper = -12.58, 1224.34
        z_lower,z_upper = -300, 300
        grid_spacing_value = 20

        x1, y1, z1 = MDAnalysis.visualization.streamlines_3D.generate_streamlines_3d(
                        'testing.gro', 'testing_filtered.xtc',
                         xmin=x_lower, xmax=x_upper,
                         ymin=y_lower, ymax=y_upper,
                         zmin=z_lower, zmax=z_upper,
                         grid_spacing=grid_spacing_value, MDA_selection = 'name PO4',
                         start_frame=2, end_frame=3, num_cores='maximum')

        x, y, z = np.mgrid[x_lower:x_upper:x1.shape[0]*1j,
                          y_lower:y_upper:y1.shape[1]*1j,
                          z_lower:z_upper:z1.shape[2]*1j]

        # plot with mayavi:
        fig = mayavi.mlab.figure(bgcolor=(1.0, 1.0, 1.0), size=(800, 800), fgcolor=(0, 0, 0))
        for z_value in np.arange(z_lower, z_upper, grid_spacing_value):
            st = mayavi.mlab.flow(x, y, z, x1, y1, z1, line_width=1,
                                  seedtype='plane', integration_direction='both')
            st.streamline_type = 'tube'
            st.tube_filter.radius = 2
            st.seed.widget.origin = np.array([ x_lower,  y_upper,   z_value])
            st.seed.widget.point1 = np.array([ x_upper, y_upper,  z_value])
            st.seed.widget.point2 = np.array([ x_lower, y_lower,  z_value])
            st.seed.widget.resolution = int(x1.shape[0])
            st.seed.widget.enabled = False
        mayavi.mlab.axes(extent = [0, 1200, 0, 1200, -300, 300])
        fig.scene.z_plus_view()
        mayavi.mlab.savefig('test_streamplot_3D.png')
        # more compelling examples can be produced for vesicles and other spherical systems

    .. image:: test_streamplot_3D.png

    See Also
    --------
    MDAnalysis.visualization.streamlines.generate_streamlines


    .. _mayavi: http://docs.enthought.com/mayavi/mayavi/
    """
    # work out the number of cores to use:
    if num_cores == 'maximum':
        num_cores = multiprocessing.cpu_count()  # use all available cores
    else:
        num_cores = num_cores  # use the value specified by the user
        # assert isinstance(num_cores,(int,long)), "The number of specified cores must (of course) be an integer."
    np.seterr(all='warn', over='raise')
    parent_cube_dictionary = {}  # collect all data from child processes here

    def log_result_to_parent(process_dict):
        parent_cube_dictionary.update(process_dict)

    #step 1: produce tuple of cartesian coordinate limits for the first frame
    #tuple_of_limits = determine_container_limits(topology_file_path = topology_file_path,trajectory_file_path =
    # trajectory_file_path,buffer_value=buffer_value)
    tuple_of_limits = (xmin, xmax, ymin, ymax, zmin, zmax)
    #step 2: produce a suitable grid (will assume that grid size / container size does not vary during simulation--or
    #  at least not beyond the buffer limit, such that this grid can be used for all subsequent frames)
    grid = produce_grid(tuple_of_limits=tuple_of_limits, grid_spacing=grid_spacing)
    #step 3: split the grid into a dictionary of cube information that can be sent to each core for processing:
    list_dictionaries_for_cores, total_cubes, num_sheets, delta_array_shape = split_grid(grid=grid, num_cores=num_cores)
    #step 3b: produce required coordinate arrays on a single core to avoid making a universe object on each core:
    start_frame_coord_array, end_frame_coord_array = produce_coordinate_arrays_single_process(topology_file_path,
                                                                                              trajectory_file_path,
                                                                                              MDA_selection,
                                                                                              start_frame, end_frame)
    #step 4: per process work using the above grid data split
    pool = multiprocessing.Pool(num_cores)
    for sub_dictionary_of_cube_data in list_dictionaries_for_cores:
        pool.apply_async(per_core_work, args=(
            start_frame_coord_array, end_frame_coord_array, sub_dictionary_of_cube_data, MDA_selection, start_frame,
            end_frame), callback=log_result_to_parent)
    pool.close()
    pool.join()
    #so, at this stage the parent process now has a single dictionary with all the cube objects updated from all
    # available cores
    #the 3D streamplot (i.e, mayavi flow() function) will require separate 3D np arrays for dx,dy,dz
    #the shape of each 3D array will unfortunately have to match the mgrid data structure (bit of a pain): (
    # num_sheets - 1, num_sheets - 1, cubes_per_column)
    cubes_per_sheet = int(float(total_cubes) / float(num_sheets - 1))
    #produce dummy zero arrays for dx,dy,dz of the appropriate shape:
    dx_array = np.zeros(delta_array_shape)
    dy_array = np.zeros(delta_array_shape)
    dz_array = np.zeros(delta_array_shape)
    #now use the parent cube dictionary to correctly substitute in dx,dy,dz values
    current_sheet = 0  # which is also the current row
    y_index_current_sheet = 0  # sub row
    z_index_current_column = 0  # column
    total_cubes_current_sheet = 0
    for cube_number in range(0, total_cubes):
        dx_array[current_sheet, y_index_current_sheet, z_index_current_column] = parent_cube_dictionary[cube_number][
            'dx']
        dy_array[current_sheet, y_index_current_sheet, z_index_current_column] = parent_cube_dictionary[cube_number][
            'dy']
        dz_array[current_sheet, y_index_current_sheet, z_index_current_column] = parent_cube_dictionary[cube_number][
            'dz']
        z_index_current_column += 1
        total_cubes_current_sheet += 1
        if z_index_current_column == delta_array_shape[2]:
            # done building current y-column so iterate y value and reset z
            z_index_current_column = 0
            y_index_current_sheet += 1
            if y_index_current_sheet == delta_array_shape[1]:  # current sheet is complete
                current_sheet += 1
                y_index_current_sheet = 0  # restart for new sheet
                z_index_current_column = 0
                total_cubes_current_sheet = 0
    # now set velocity component values greater than a certain cutoff to 0,
    # because they tend to reflect spurious values (i.e., PBC jumping)
    dx_array[abs(dx_array) >= maximum_delta_magnitude] = 1.0
    dy_array[abs(dy_array) >= maximum_delta_magnitude] = 1.0
    dz_array[abs(dz_array) >= maximum_delta_magnitude] = 1.0
    return (dx_array, dy_array, dz_array)
