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
# Python implementation of FORTRAN HELANAL
# [Bansal et al, J Biomol Struct Dyn. 17 (2000), 811]
# Copyright (c) 2009 Benjamin Hall <benjamin.hall@bioch.ox.ac.uk>
# Published under the GNU General Public License v2 or higher
#
# Copyright (c) 2011 Oliver Beckstein <orbeckst@gmail.com>
# Integrated into MDAnalysis and NumPy-fied


"""
HELANAL (Deprecated) --- analysis of protein helices
====================================================

:Author:  Benjamin Hall <benjamin.a.hall@ucl.ac.uk>, Oliver Beckstein, Xavier Deupi
:Year:    2009, 2011, 2013
:License: GNU General Public License v2 (or higher)

.. note::

    This module was deprecated in 1.0 and will be removed in 2.0.
    Please use MDAnalysis.analysis.helix_analysis instead
    (available in 2.0.0).


The :mod:`MDAnalysis.analysis.helanal` module is a Python implementation of the
HELANAL_ algorithm [Bansal2000]_ in `helanal.f`_, which is also available
through the `HELANAL webserver`_.

Please cite the paper [Bansal2000]_ (and possibly [Kumar1996]_ and
[Kumar1998]_) in published work when using
:mod:`~MDAnalysis.analysis.helanal.helanal_trajectory` or
:mod:`~MDAnalysis.analysis.helanal.helanal_main`.

HELANAL_ quantifies the geometry of helices in proteins on the basis of their Cα
atoms alone. It can extract the helices from the structure files and then
characterises the overall geometry of each helix as being linear, curved or
kinked, in terms of its local structural features, viz. local helical twist and
rise, virtual torsion angle, local helix origins and bending angles between
successive local helix axes. Even helices with large radius of curvature are
unambiguously identified as being linear or curved. The program can also be
used to differentiate a kinked helix and other motifs, such as helix-loop-helix
or a helix-turn-helix (with a single residue linker) with the help of local
bending angles. In addition to these, the program can also be used to
characterise the helix start and end as well as other types of secondary
structures.

.. _HELANAL: http://www.ccrnp.ncifcrf.gov/users/kumarsan/HELANAL/helanal.html
.. _Helanal webserver: http://nucleix.mbu.iisc.ernet.in/helanal/helanal.shtml
.. `helanal.f`: http://www.webcitation.org/5y1RpVJtF
.. _`helanal.f`: http://www.ccrnp.ncifcrf.gov/users/kumarsan/HELANAL/helanal.f

Background
----------

From the HELANAL_ home page:

HELANAL_ can be used to characterize the geometry of helices with a minimum 9
residues. The geometry of an alpha helix is characterized by computing local
helix axes and local helix origins for four contiguous C-Alpha atoms, using the
procedure of Sugeta and Miyazawa [Sugeta1967]_ and sliding this window over
the length of the helix in steps of one C-alpha atom.

The angles between successive local helix axes can identify *local bends* or
*kinks* as well as occurrence of *smooth curvature* in the helix. A matrix, whose
elements *M(I, J)* are the *bending angles* between local helix axes *I* and *J*,
is obtained to get an idea about the overall geometry of the helix.

*Unit twist* and *unit height* of the alpha helix are also computed to analyze the
uniformity of the helix. The *local helix origins* trace out the path described
by the helix in three dimensional space. The local helix origins are reoriented
in *X-Y* plane and the reoriented points are used to fit a circle as well as a
line, by least squares method. Based on the relative goodness of line and
circle fit to local helix origins, the helix is *classified as being linear or
curved*. A helix is classified as being *kinked*, if at least one local bending
angle in the middle of the helix is greater than 20 degrees.


References
----------

.. [Sugeta1967] Sugeta, H. and Miyazawa, T. 1967. General method for
   calculating helical parameters of polymer chains from bond lengths, bond
   angles and internal rotation angles. *Biopolymers* 5 673 - 679

.. [Kumar1996] Kumar, S. and Bansal, M. 1996. Structural and sequence
   characteristics of long alpha-helices in globular proteins. *Biophysical
   Journal* 71(3):1574-1586.

.. [Kumar1998] Kumar, S. and Bansal, M. 1998. Geometrical and sequence
   characteristics of alpha helices in globular proteins. *Biophysical Journal*
   75(4):1935-1944.

.. [Bansal2000] Bansal M, Kumar S, Velavan R. 2000.
   HELANAL - A program to characterise helix geometry in proteins.
   *J Biomol Struct Dyn.*  17(5):811-819.

Functions
---------

.. autofunction:: helanal_trajectory
.. autofunction:: helanal_main

"""
from __future__ import print_function, division, absolute_import
from six.moves import range, zip

import os

import numpy as np

import MDAnalysis
from ..lib.log import ProgressBar
from ..lib import mdamath, util

import warnings
import logging
logger = logging.getLogger("MDAnalysis.analysis.helanal")

warnings.warn("This module is deprecated as of MDAnalysis version 1.0. "
              "It will be removed in MDAnalysis version 2.0."
              "Please use (available in 2.0.0) "
              "MDAnalysis.analysis.helix_analysis instead.",
              category=DeprecationWarning)

def center(coordinates):
    """Return the geometric center (centroid) of the coordinates.

    Coordinates must be "list of cartesians", i.e. a Nx3 array.
    """
    return np.mean(coordinates, axis=0)

def vecnorm(a):
    """Return a/|a|"""
    return a / mdamath.norm(a)

def wrapangle(angle):
    """Wrap angle (in radians) to be within -pi < angle =< pi"""
    if angle > np.pi:
        angle -= 2 * np.pi
    elif angle <= -np.pi:
        angle += 2 * np.pi
    return angle

def sample_sd(a, dummy):
    return np.std(a, ddof=1)

def mean_abs_dev(a, mean_a=None):
    if mean_a is None:
        mean_a = np.mean(a)
    return np.mean(np.fabs(a - mean_a))


@util.deprecate(release="1.0.1", remove="2.0",
                message="In 2.0 use the upcoming helix_analysis module")
def helanal_trajectory(universe, select="name CA",
                       begin=None, finish=None,
                       matrix_filename="bending_matrix.dat",
                       origin_pdbfile="origin.pdb",
                       summary_filename="summary.txt",
                       screw_filename="screw.xvg",
                       tilt_filename="local_tilt.xvg",
                       fitted_tilt_filename="fit_tilt.xvg",
                       bend_filename="local_bend.xvg",
                       twist_filename="unit_twist.xvg",
                       prefix="helanal_", ref_axis=None,
                       verbose=False):
    """Perform HELANAL helix analysis on all frames in `universe`.

    Parameters
    ----------
    universe : Universe
    select : str (optional)
        selection string that selects Calpha atoms [``"name CA"``]
    begin : float (optional)
        start analysing for time (ps) >= *begin*; ``None`` starts from the
        beginning [``None``]
    finish : float (optional)
        stop analysis for time (ps) =< *finish*; ``None`` goes to the
        end of the trajectory [``None``]
    matrix_filename : str (optional)
        Output file- bending matrix [``"bending_matrix.dat"``]
    origin_pdbfile : str (optional)
        Output file- origin pdb file [``"origin.pdb"``]
    summary_filename : str (optional)
        Output file- all of the basic data [``"summary.txt"``]
    screw_filename : str (optional)
        Output file- local tilts of individual residues from 2 to n-1
        [``"screw.xvg"``]
    tilt_filename : str (optional)
        Output file- tilt of line of best fit applied to origin axes
        [``"local_tilt.xvg"``]
    bend_filename : str (optional)
        Output file- local bend angles between successive local helix axes
        [``"local_bend.xvg"``]
    twist_filename : str (optional)
        Output file- local unit twist between successive helix turns
        [``"unit_twist.xvg"``]
    prefix : str (optional)
        Prefix to add to all output file names; set to ``None`` to disable
        [``"helanal__"``]
    ref_axis : array_like (optional)
        Calculate tilt angle relative to the axis; if ``None`` then ``[0,0,1]``
        is chosen [``None``]
    verbose : bool (optional)
        Toggle diagnostic outputs. [``True``]

    Raises
    ------
    ValueError
          If the specified start (begin) time occurs after the end of the
          trajectory object.
          If the specified finish time precedes the specified start time or
          current time stamp of trajectory object.

    Notes
    -----
    Only a single helix is analyzed. Use the selection to specify the helix,
    e.g. with "name CA and resid 1:20" or use start=1, stop=20.


    .. versionchanged:: 0.13.0
       New `quiet` keyword to silence frame progress output and most of the
       output that used to be printed to stdout is now logged to the logger
       *MDAnalysis.analysis.helanal* (at logelevel *INFO*).

    .. versionchanged:: 0.16.0
       Removed the `start` and `end` keywords for selecting residues because this can
       be accomplished more transparently with `select`. The first and last resid
       are directly obtained from the selection.

    .. deprecated:: 0.16.0
       The `quiet` keyword argument is deprecated in favor of the new
       `verbose` one.

    .. versionchanged:: 0.20.0
       ProgressMeter now iterates over the number of frames analysed.

    .. versionchanged:: 1.0.0
       Changed `selection` keyword to `select`

    .. deprecated:: 1.0.1
       helanal_trajectory is deprecated and will be removed in 2.0.0.
       Please use helix_analysis.HELANAL instead (available in 2.0.0).
    """

    if ref_axis is None:
        ref_axis = np.array([0., 0., 1.])
    else:
        # enable MDA API so that one can use a tuple of atoms or AtomGroup with
        # two atoms
        ref_axis = np.asarray(ref_axis)

    ca = universe.select_atoms(select)
    start, end = ca.resids[[0, -1]]
    trajectory = universe.trajectory

    # Validate user supplied begin / end times
    traj_end_time = trajectory.ts.time + trajectory.totaltime

    if begin is not None:
        if traj_end_time < begin:
            # Begin occurs after the end of the trajectory, throw error
            msg = ("The input begin time ({0} ps) occurs after the end "
                   "of the trajectory ({1} ps)".format(begin, traj_end_time))
            raise ValueError(msg)
        elif trajectory.ts.time > begin:
            # Begin occurs before trajectory start, warn and reset
            msg = ("The input begin time ({0} ps) precedes the starting "
                   "trajectory time --- Setting starting frame to 0".format(
                    begin))
            warnings.warn(msg)
            logger.warning(msg)
            start_frame = None
        else:
            start_frame = int(np.ceil((begin - trajectory.ts.time)
                                      / trajectory.ts.dt))
    else:
        start_frame = None

    if finish is not None:
        if (begin is not None) and (begin > finish):
            # finish occurs before begin time
            msg = ("The input finish time ({0} ps) precedes the input begin "
                   "time ({1} ps)".format(finish, begin))
            raise ValueError(msg)
        elif trajectory.ts.time > finish:
            # you'd be starting with a finish time(in ps) that has already
            # passed or is not available
            msg = ("The input finish time ({0} ps) precedes the current "
                   "trajectory time ({1} ps)".format(finish, trajectory.time))
            raise ValueError(msg)
        elif traj_end_time < finish:
            # finish time occurs after the end of trajectory, warn
            msg = ("The input finish time ({0} ps) occurs after the end of "
                   "the trajectory ({1} ps). Finish time will be set to the "
                   "end of the trajectory".format(finish, traj_end_time))
            warnings.warn(msg)
            logger.warning(msg)
            end_frame = None
        else:
            # To replicate the original behaviour of break when
            # trajectory.time > finish, we add 1 here.
            end_frame = int(np.floor((finish - trajectory.ts.time)
                            // trajectory.ts.dt) + 1)
    else:
        end_frame = None

    start_frame, end_frame, frame_step = trajectory.check_slice_indices(
                                          start_frame, end_frame, 1)
    n_frames = len(range(start_frame, end_frame, frame_step))

    if start is not None and end is not None:
        logger.info("Analysing from residue %d to %d", start, end)
    elif start is not None and end is None:
        logger.info("Analysing from residue %d to the C termini", start)
    elif start is None and end is not None:
        logger.info("Analysing from the N termini to %d", end)
    logger.info("Analysing %d/%d residues", ca.n_atoms, universe.atoms.n_residues)

    if prefix is not None:
        prefix = str(prefix)
        matrix_filename = prefix + matrix_filename
        origin_pdbfile = prefix + origin_pdbfile
        summary_filename = prefix + summary_filename
        screw_filename = prefix + screw_filename
        tilt_filename = prefix + tilt_filename
        fitted_tilt_filename = prefix + fitted_tilt_filename
        bend_filename = prefix + bend_filename
        twist_filename = prefix + twist_filename
    backup_file(matrix_filename)
    backup_file(origin_pdbfile)
    backup_file(summary_filename)
    backup_file(screw_filename)
    backup_file(tilt_filename)
    backup_file(fitted_tilt_filename)
    backup_file(bend_filename)
    backup_file(twist_filename)

    global_height = []
    global_twist = []
    global_rnou = []
    global_bending = []
    global_bending_matrix = []
    global_tilt = []
    global_fitted_tilts = []
    global_screw = []

    for ts in ProgressBar(trajectory[start_frame:end_frame:frame_step],
                          verbose=verbose, desc="Helix analysis"):

        frame = ts.frame
        ca_positions = ca.positions
        twist, bending_angles, height, rnou, origins, local_helix_axes, local_screw_angles = \
            main_loop(ca_positions, ref_axis=ref_axis)

        origin_pdb(origins, origin_pdbfile)

        #calculate local bending matrix( it is looking at all i, j combinations)
        if len(global_bending_matrix) == 0:
            global_bending_matrix = [[[] for item in local_helix_axes] for item in local_helix_axes]

        for i in range(len(local_helix_axes)):
            for j in range(i + 1, len(local_helix_axes)):
                angle = np.rad2deg(np.arccos(np.dot(local_helix_axes[i], local_helix_axes[j])))
                global_bending_matrix[i][j].append(angle)
                #global_bending_matrix[j][i].append(angle)
                #global_bending_matrix[i][i].append(0.)

        fit_vector, fit_tilt = vector_of_best_fit(origins)
        global_height += height
        global_twist += twist
        global_rnou += rnou
        #global_screw.append(local_screw_angles)
        global_fitted_tilts.append(np.rad2deg(fit_tilt))

        #print out rotations across the helix to a file
        with open(twist_filename, "a") as twist_output:
            print(frame, end=' ', file=twist_output)
            for loc_twist in twist:
                print(loc_twist, end=' ', file=twist_output)
            print("", file=twist_output)

        with open(bend_filename, "a") as bend_output:
            print(frame, end=' ', file=bend_output)
            for loc_bend in bending_angles:
                print(loc_bend, end=' ', file=bend_output)
            print("", file=bend_output)

        with open(screw_filename, "a") as rot_output:
            print(frame, end=' ', file=rot_output)
            for rotation in local_screw_angles:
                print(rotation, end=' ', file=rot_output)
            print("", file=rot_output)

        with open(tilt_filename, "a") as tilt_output:
            print(frame, end=' ', file=tilt_output)
            for tilt in local_helix_axes:
                print(np.rad2deg(mdamath.angle(tilt, ref_axis)),
                      end=' ', file=tilt_output)
            print("", file=tilt_output)

        with open(fitted_tilt_filename, "a") as tilt_output:
            print(frame, np.rad2deg(fit_tilt), file=tilt_output)

        if len(global_bending) == 0:
            global_bending = [[] for item in bending_angles]
            #global_tilt = [ [] for item in local_helix_axes ]
        for store, tmp in zip(global_bending, bending_angles):
            store.append(tmp)
        #for store,tmp in zip(global_tilt,local_helix_axes): store.append(mdamath.angle(tmp,ref_axis))

    twist_mean, twist_sd, twist_abdev = stats(global_twist)
    height_mean, height_sd, height_abdev = stats(global_height)
    rnou_mean, rnou_sd, rnou_abdev = stats(global_rnou)
    ftilt_mean, ftilt_sd, ftilt_abdev = stats(global_fitted_tilts)

    bending_statistics = [stats(item) for item in global_bending]
    #tilt_statistics =    [ stats(item) for item in global_tilt]

    bending_statistics_matrix = [[stats(col) for col in row] for row in global_bending_matrix]
    with open(matrix_filename, 'w') as mat_output:
        print("Mean", file=mat_output)
        for row in bending_statistics_matrix:
            for col in row:
                formatted_angle = "{0:6.1f}".format(col[0])
                print(formatted_angle, end=' ', file=mat_output)
            print('', file=mat_output)

        print('\nSD', file=mat_output)
        for row in bending_statistics_matrix:
            for col in row:
                formatted_angle = "{0:6.1f}".format(col[1])
                print(formatted_angle, end=' ', file=mat_output)
            print('', file=mat_output)

        print("\nABDEV", file=mat_output)
        for row in bending_statistics_matrix:
            for col in row:
                formatted_angle = "{0:6.1f}".format(col[2])
                print(formatted_angle, end=' ', file=mat_output)
            print('', file=mat_output)

    logger.info("Height: %g  SD: %g  ABDEV: %g  (Angstroem)", height_mean, height_sd, height_abdev)
    logger.info("Twist: %g  SD: %g  ABDEV: %g", twist_mean, twist_sd, twist_abdev)
    logger.info("Residues/turn: %g  SD: %g  ABDEV: %g", rnou_mean, rnou_sd, rnou_abdev)
    logger.info("Fitted tilt: %g  SD: %g  ABDEV: %g", ftilt_mean, ftilt_sd, ftilt_abdev)
    logger.info("Local bending angles:")
    residue_statistics = list(zip(*bending_statistics))
    measure_names = ["Mean ", "SD   ", "ABDEV"]
    if start is None:
        output = " ".join(["{0:8d}".format(item)
                           for item in range(4, len(residue_statistics[0]) + 4)])
    else:
        output = " ".join(["{0:8d}".format(item)
                           for item in range(start + 3, len(residue_statistics[0]) + start + 3)])
    logger.info("ResID %s", output)
    for measure, name in zip(residue_statistics, measure_names):
        output = str(name) + " "
        output += " ".join(["{0:8.1f}".format(residue) for residue in measure])
        logger.info(output)

    with open(summary_filename, 'w') as summary_output:
        print("Height:", height_mean, "SD", height_sd, "ABDEV", height_abdev, '(nm)', file=summary_output)
        print("Twist:", twist_mean, "SD", twist_sd, "ABDEV", twist_abdev,
              file=summary_output)
        print("Residues/turn:", rnou_mean, "SD", rnou_sd, "ABDEV", rnou_abdev,
              file=summary_output)
        print("Local bending angles:", file=summary_output)
        residue_statistics = list(zip(*bending_statistics))
        measure_names = ["Mean ", "SD   ", "ABDEV"]
        print("ResID", end='', file=summary_output)
        if start is None:
            for item in range(4, len(residue_statistics[0]) + 4):
                output = "{0:8d}".format(item)
                print(output, end=' ', file=summary_output)
        else:
            for item in range(start + 3, len(residue_statistics[0]) + start + 3):
                output = "{0:8d}".format(item)
                print(output, end=' ', file=summary_output)
        print('', file=summary_output)

        for measure, name in zip(residue_statistics, measure_names):
            print(name, end=' ', file=summary_output)
            for residue in measure:
                output = "{0:8.1f}".format(residue)
                print(output, end=' ', file=summary_output)
            print('', file=summary_output)


def tilt_correct(number):
    """Changes an angle (in degrees) so that it is between 0º and 90º"""
    if number < 90.:
        return number
    else:
        return 180. - number


def backup_file(filename):
    if os.path.exists(filename):
        target_name = "#" + filename
        failure = True
        if not os.path.exists(target_name):
            os.rename(filename, target_name)
            failure = False
        else:
            for i in range(20):
                alt_target_name = target_name + "." + str(i)
                if os.path.exists(alt_target_name):
                    continue
                else:
                    os.rename(filename, alt_target_name)
                    failure = False
                    break
        if failure:
            raise IOError("Too many backups. Clean up and try again")


def stats(some_list):
    if len(some_list) == 0:
        return [0, 0, 0]
    list_mean = np.mean(some_list)
    list_sd = sample_sd(some_list, list_mean)
    list_abdev = mean_abs_dev(some_list, list_mean)
    return [list_mean, list_sd, list_abdev]


@util.deprecate(release="1.0.1", remove="2.0",
                message="In 2.0 use the upcoming helix_analysis module")
def helanal_main(pdbfile, select="name CA", ref_axis=None):
    """Simple HELANAL_ run on a single frame PDB/GRO.

    Computed data are returned as a dict and also logged at level INFO to the
    logger *MDAnalysis.analysis.helanal*. A simple way to enable a logger is to
    use :func:`~MDAnalysis.lib.log.start_logging`.

    Parameters
    ----------
    pdbfile : str
        filename of the single-frame input file
    select : str (optional)
        selection string, default is "name CA" to select all C-alpha atoms.
    ref_axis : array_like (optional)
        Calculate tilt angle relative to the axis; if ``None`` then ``[0,0,1]``
        is chosen [``None``]

    Returns
    -------
    result : dict
       The `result` contains keys
        * Height: mean, stdev, abs dev
        * Twist: mean, stdev, abs dev
        * Residues/turn: mean, stdev, abs dev
        * Local bending angles: array for computed angles (per residue)
        * Unit twist angles: array for computed angles (per residue)
        * Best fit tilt
        * Rotation angles: local screw angles (per residue)


    Notes
    -----
    Only a single helix is analyzed. Use the selection to specify the
    helix, e.g. with "name CA and resid 1:20".


    Example
    -------
    Analyze helix 8 in AdK (PDB 4AKE); the standard logger is started and
    writes output to the file ``MDAnalysis.log``::

       MDAnalysis.start_logging()
       data = MDAnalysis.analysis.helanal_main("4ake_A.pdb", select="name CA and resnum 161-187")


    .. versionchanged:: 0.13.0
       All output is returned as a dict and logged to the logger
       *MDAnalysis.analysis.helanal* instead of being printed to stdout.

    .. versionchanged:: 0.16.0
       Removed the `start` and `end` keywords for selecting residues because this can
       be accomplished more transparently with `select`.

    .. versionchanged:: 1.0.0
       Changed `selection` keyword to `select`

    .. deprecated:: 1.0.1
       helanal_trajectory is deprecated and will be removed in 2.0.0.
       Please use helix_analysis.HELANAL instead (available in 2.0.0).

    """

    universe = MDAnalysis.Universe(pdbfile)
    ca = universe.select_atoms(select)

    logger.info("Analysing %d/%d residues", ca.n_atoms, universe.atoms.n_residues)

    twist, bending_angles, height, rnou, origins, local_helix_axes, local_screw_angles = \
        main_loop(ca.positions, ref_axis=ref_axis)

    #TESTED- origins are correct
    #print current_origin
    #print origins

    max_angle = np.max(bending_angles)
    mean_angle = np.mean(bending_angles)
    #sd calculated using n-1 to replicate original fortran- assumes a limited sample so uses the sample standard
    # deviation
    sd_angle = sample_sd(bending_angles, mean_angle)
    mean_absolute_deviation_angle = mean_abs_dev(bending_angles, mean_angle)
    #TESTED- stats correct
    #print max_angle, mean_angle, sd_angle, mean_absolute_deviation_angle

    #calculate local bending matrix(now it is looking at all i, j combinations)
    # (not used for helanal_main())
#     for i in local_helix_axes:
#         for j in local_helix_axes:
#             if (i == j).all():
#                 angle = 0.
#             else:
#                 angle = np.rad2deg(np.arccos(np.dot(i, j)))
#             #string_angle = "%6.0f\t" % angle
#             #print string_angle,
#             #print ''
#             #TESTED- local bending matrix!

    #Average helical parameters
    mean_twist = np.mean(twist)
    sd_twist = sample_sd(twist, mean_twist)
    abdev_twist = mean_abs_dev(twist, mean_twist)
    #TESTED-average twists
    #print mean_twist, sd_twist, abdev_twist
    mean_rnou = np.mean(rnou)
    sd_rnou = sample_sd(rnou, mean_rnou)
    abdev_rnou = mean_abs_dev(rnou, mean_rnou)
    #TESTED-average residues per turn
    #print mean_rnou, sd_rnou, abdev_rnou
    mean_height = np.mean(height)
    sd_height = sample_sd(height, mean_height)
    abdev_height = mean_abs_dev(height, mean_height)
    #TESTED- average rises

    #calculate best fit vector and tilt of said vector
    fit_vector, fit_tilt = vector_of_best_fit(origins)

    data = {
        'Height': np.array([mean_height, sd_height, abdev_height]),
        'Twist': np.array([mean_twist, sd_twist, abdev_twist]),
        'Residues/turn': np.array([mean_rnou, sd_rnou, abdev_rnou]),
        'Local bending angles': np.asarray(bending_angles),
        'Unit twist angles': np.asarray(twist),
        'Best fit tilt': fit_tilt,
        'Rotation Angles': np.asarray(local_screw_angles),
        }

    logger.info("Height: %g  SD: %g  ABDEV: %g  (Angstroem)", mean_height, sd_height, abdev_height)
    logger.info("Twist: %g  SD: %g  ABDEV: %g", mean_twist, sd_twist, abdev_twist)
    logger.info("Residues/turn: %g  SD: %g  ABDEV: %g", mean_rnou, sd_rnou, abdev_rnou)

    output = " ".join(["{0:8.1f}\t".format(angle) for angle in bending_angles])
    logger.info("Local bending angles: %s", output)

    output = " ".join(["{0:8.1f}\t".format(twist_ang) for twist_ang in twist])
    logger.info("Unit twist angles: %s", output)

    logger.info("Best fit tilt: %g", fit_tilt)

    output = " ".join(["{0:.1f}".format(item) for item in local_screw_angles])
    logger.info("Rotation Angles from 1 to n-1 (local screw angles): %s", output)

    return data

def origin_pdb(origins, pdbfile):
    """Write origins to PDB (multi-frame).

    This PDB can be loaded together with the trajectory into, e.g. VMD_, to
    view the helix axis together with all the atoms.
    """
    with open(pdbfile, 'a') as output:
        i = 1
        for xyz in origins:
            tmp = "ATOM    {0:3d}  CA  ALA   {1:3d}    {2:8.3f}{3:8.3f}{4:8.3f}  1.00  0.00".format(i, i, xyz[0], xyz[1], xyz[2])
            print(tmp, file=output)
            i += 1
        print("TER\nENDMDL", file=output)


def main_loop(positions, ref_axis=None):
    # rewrite in cython?

    if ref_axis is None:
        ref_axis = np.array([0., 0., 1.])
    else:
        ref_axis = np.asarray(ref_axis)
    twist = []
    rnou = []
    height = []
    origins = [[0., 0., 0.] for item in positions[:-2]]
    local_helix_axes = []
    location_rotation_vectors = []
    for i in range(len(positions) - 3):
        vec12 = positions[i + 1] - positions[i]
        vec23 = positions[i + 2] - positions[i + 1]
        vec34 = positions[i + 3] - positions[i + 2]

        dv13 = vec12 - vec23
        dv24 = vec23 - vec34

        #direction of the local helix axis
        current_uloc = vecnorm(np.cross(dv13, dv24))
        local_helix_axes.append(current_uloc)

        #TESTED- Axes correct
        #print current_uloc

        dmag = mdamath.norm(dv13)
        emag = mdamath.norm(dv24)

        costheta = np.dot(dv13, dv24) / (dmag * emag)
        #rnou is the number of residues per turn
        current_twist = np.arccos(costheta)
        twist.append(np.rad2deg(current_twist))
        rnou.append(2 * np.pi / current_twist)
        #radius of local helix cylinder radmag

        costheta1 = 1.0 - costheta
        radmag = (dmag * emag) ** 0.5 / (2 * costheta1)

        #Height of local helix cylinder
        current_height = np.dot(vec23, current_uloc)
        height.append(current_height)
        #TESTED- Twists etc correct
        #print current_twist*180/np.pi, 2*np.pi/current_twist, height

        dv13 = vecnorm(dv13)
        dv24 = vecnorm(dv24)

        #record local rotation
        location_rotation_vectors.append(dv13)

        rad = [radmag * item for item in dv13]
        current_origin = [(item[0] - item[1]) for item in zip(positions[i + 1], rad)]
        origins[i] = current_origin

        #TESTED- origins are correct
        #print current_origin

        rad = [radmag * item for item in dv24]
        current_origin = [(item[0] - item[1]) for item in zip(positions[i + 2], rad)]
        origins[i + 1] = current_origin
    #Record final rotation vector
    location_rotation_vectors.append(dv24)

    #local bending angles (eg i > i+3, i+3 > i+6)

    bending_angles = [0 for item in range(len(local_helix_axes) - 3)]
    for axis in range(len(local_helix_axes) - 3):
        angle = np.arccos(np.dot(local_helix_axes[axis], local_helix_axes[axis + 3]))
        bending_angles[axis] = np.rad2deg(angle)
        #TESTED- angles are correct
        #print np.rad2deg(angle)

    local_screw_angles = []
    #Calculate rotation angles for (+1) to (n-1)
    fit_vector, fit_tilt = vector_of_best_fit(origins)
    for item in location_rotation_vectors:
        local_screw_tmp = np.rad2deg(rotation_angle(fit_vector, ref_axis, item))
        #print local_screw_tmp
        local_screw_angles.append(local_screw_tmp)

    return twist, bending_angles, height, rnou, origins, local_helix_axes, local_screw_angles


def rotation_angle(helix_vector, axis_vector, rotation_vector):
    reference_vector = np.cross(np.cross(helix_vector, axis_vector), helix_vector)
    second_reference_vector = np.cross(axis_vector, helix_vector)
    screw_angle = mdamath.angle(reference_vector, rotation_vector)
    alt_screw_angle = mdamath.angle(second_reference_vector, rotation_vector)
    updown = np.cross(reference_vector, rotation_vector)

    if not (np.pi < screw_angle < 3 * np.pi / 4):
        if screw_angle < np.pi / 4 and alt_screw_angle < np.pi / 2:
            screw_angle = np.pi / 2 - alt_screw_angle
        elif screw_angle < np.pi / 4 and alt_screw_angle > np.pi / 2:
            screw_angle = alt_screw_angle - np.pi / 2
        elif screw_angle > 3 * np.pi / 4 and alt_screw_angle < np.pi / 2:
            screw_angle = np.pi / 2 + alt_screw_angle
        elif screw_angle > 3 * np.pi / 4 and alt_screw_angle > np.pi / 2:
            screw_angle = 3 * np.pi / 2 - alt_screw_angle
        else:
            logger.debug("Big Screw Up: screw_angle=%g degrees", np.rad2deg(screw_angle))

    if mdamath.norm(updown) == 0:
        logger.warning("PROBLEM (vector is at 0 or 180)")

    helix_dot_rehelix = mdamath.angle(updown, helix_vector)

    #if ( helix_dot_rehelix < np.pi/2 and helix_dot_rehelix >= 0 )or helix_dot_rehelix <-np.pi/2:
    if (-np.pi / 2 < helix_dot_rehelix < np.pi / 2) or (helix_dot_rehelix > 3 * np.pi / 2):
        screw_angle = -screw_angle

    return screw_angle


def vector_of_best_fit(origins):
    origins = np.asarray(origins)
    centroids = center(origins)
    M = origins - centroids
    A = np.dot(M.transpose(), M)
    u, s, vh = np.linalg.linalg.svd(A)
    vector = vh[0]
    #Correct vector to face towards first residues
    rough_helix = origins[0] - centroids
    agreement = mdamath.angle(rough_helix, vector)
    if not (-np.pi / 2 < agreement < np.pi / 2):
        vector *= -1
    best_fit_tilt = mdamath.angle(vector, [0, 0, 1])
    return vector, best_fit_tilt
