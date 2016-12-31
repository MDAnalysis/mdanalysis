# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.mdanalysis.org
# Copyright (c) 2006-2016 The MDAnalysis Development Team and contributors
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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
"""
Distance Matrix calculation --- :mod:`MDAnalysis.analysis.ensemble.confdistmatrix`
==================================================================================


The module contains a base class to easily compute, using
parallelization and shared memory, matrices of conformational
distance between the structures stored as frames in a Universe. A
class to compute an RMSD matrix in such a way is also available.

:Author: Matteo Tiberti, Wouter Boomsma, Tone Bengtsen

.. versionadded:: 0.16.0

"""

import numpy as np
from getpass import getuser
from socket import gethostname
from datetime import datetime
from time import sleep
import logging

from sklearn.externals.joblib import Parallel, delayed

from ...core.universe import Universe

from ..align import rotation_matrix

from .cutils import PureRMSD
from .utils import TriangularMatrix, trm_indeces



def conformational_distance_matrix(ensemble,
                                   conf_dist_function, selection="",
                                   superimposition_selection="", ncores=1, pairwise_align=True,
                                   mass_weighted=True, metadata=True, verbose=False):
    """
    Run the conformational distance matrix calculation.
    args and kwargs are passed to conf_dist_function.

    Parameters
    ----------

    ensemble : Universe object
        Universe object for which the conformational distance matrix will
        be computed.

    conf_dist_function : function object
        Function that fills the matrix with conformational distance
        values. See set_rmsd_matrix_elements for an example.

    pairwise_align : bool
        Whether to perform pairwise alignment between conformations.
        Default is True (do the superimposition)

    mass_weighted : bool
        Whether to perform mass-weighted superimposition and metric
        calculation. Default is True.

    metadata : bool
        Whether to build a metadata dataset for the calculated matrix.
        Default is True.

    ncores : int
        Number of cores to be used for parallel calculation
        Default is 1.

    Returns
    -------

    conf_dist_matrix : encore.utils.TriangularMatrix object
        Conformational distance matrix in triangular representation.

    """

    # Decide how many cores have to be used. Since the main process is
    # stopped while the workers do their job, ncores workers will be
    # spawned.

    if ncores < 1:
        ncores = 1

    # framesn: number of frames
    framesn = len(ensemble.trajectory.timeseries(
        ensemble.select_atoms(selection), format='fac'))

    # Prepare metadata recarray
    if metadata:
        metadata = np.array([(gethostname(),
                           getuser(),
                           str(datetime.now()),
                           ensemble.filename,
                           framesn,
                           pairwise_align,
                           selection,
                           mass_weighted)],
                         dtype=[('host', object),
                                ('user', object),
                                ('date', object),
                                ('topology file', object),
                                ('number of frames', int),
                                ('pairwise superimposition', bool),
                                ('superimposition subset', object),
                                ('mass-weighted', bool)])

    # Prepare alignment subset coordinates as necessary

    rmsd_coordinates = ensemble.trajectory.timeseries(
            ensemble.select_atoms(selection),
            format='fac')

    if pairwise_align:
        if superimposition_selection:
            subset_selection = superimposition_selection
        else:
            subset_selection = selection

        fitting_coordinates = ensemble.trajectory.timeseries(
            ensemble.select_atoms(subset_selection),
            format='fac')
    else:
        fitting_coordinates = None

    # Prepare masses as necessary
    if mass_weighted:
        masses = ensemble.select_atoms(selection).masses.astype(np.float64)
        if pairwise_align:
            subset_masses = ensemble.select_atoms(subset_selection).masses.astype(np.float64)
        else:
            subset_masses = None
    else:
        masses = np.ones((ensemble.trajectory.timeseries(
            ensemble.select_atoms(selection))[0].shape[0])).astype(np.float64)
        if pairwise_align:
            subset_masses = np.ones((fit_coords[0].shape[0])).astype(np.float64)
        else:
            subset_masses = None

    # Allocate for output matrix
    matsize = framesn * (framesn + 1) / 2
    distmat = np.empty(matsize, np.float64)


    # Initialize workers. Simple worker doesn't perform fitting,
    # fitter worker does.
    indices = trm_indeces((0, 0), (framesn - 1, framesn - 1))
    Parallel(n_jobs=ncores, verbose=verbose)(delayed(conf_dist_function)(
        element,
        rmsd_coordinates,
        distmat,
        masses,
        fitting_coordinates,
        subset_masses,
        masses) for element in indices)


    # When the workers have finished, return a TriangularMatrix object
    return TriangularMatrix(distmat, metadata=metadata)


def set_rmsd_matrix_elements(tasks, coords, rmsdmat, masses, fit_coords=None,
                             fit_masses=None, pbar_counter=None, *args, **kwargs):

    '''
    RMSD Matrix calculator

    Parameters
    ----------

    tasks : iterator of int of length 2
        Given a triangular matrix, this function will calculate RMSD
        values from element tasks[0] to tasks[1]. Since the matrix
        is triangular, the trm_indeces matrix automatically
        calculates the corrisponding i,j matrix indices.
        The matrix is written as an array in a row-major
        order (see the TriangularMatrix class for details).

        If fit_coords and fit_masses are specified, the structures
        will be superimposed before calculating RMSD, and fit_coords and fit_masses
        will be used to place both structures at their center of mass and
        compute the rotation matrix. In this case, both fit_coords and fit_masses
        must be specified.

    coords : numpy.array
        Array of the ensemble coordinates

    masses : numpy.array
        Array of atomic masses, having the same order as the
        coordinates array

    rmsdmat : encore.utils.TriangularMatrix
        Memory-shared triangular matrix object

    fit_coords : numpy.array or None
        Array of the coordinates used for fitting

    fit_masses : numpy.array
        Array of atomic masses, having the same order as the
        fit_coords array
        '''
    i, j = tasks

    if fit_coords is None and fit_masses is None:
        summasses = np.sum(masses)
        rmsdmat[(i + 1) * i / 2 + j] = PureRMSD(coords[i],
                                                coords[j],
                                                coords[j].shape[0],
                                                masses,
                                                summasses)

    elif fit_coords is not None and fit_coords is not None:
        summasses = np.sum(masses)
        subset_weights = np.asarray(fit_masses) / np.mean(fit_masses)
        com_i = np.average(fit_coords[i], axis=0,
                           weights=fit_masses)
        translated_i = coords[i] - com_i
        subset1_coords = fit_coords[i] - com_i
        com_j = np.average(fit_coords[j], axis=0,
                           weights=fit_masses)
        translated_j = coords[j] - com_j
        subset2_coords = fit_coords[j] - com_j
        rotamat = rotation_matrix(subset1_coords, subset2_coords,
                                  subset_weights)[0]
        rotated_i = np.transpose(np.dot(rotamat, np.transpose(translated_i)))
        rmsdmat[(i + 1) * i / 2 + j] = PureRMSD(
            rotated_i.astype(np.float64), translated_j.astype(np.float64),
            coords[j].shape[0], masses, summasses)
    else:
        raise TypeError("Both fit_coords and fit_masses must be specified \
                        if one of them is given")


def get_distance_matrix(ensemble,
                        selection="name CA",
                        load_matrix=None,
                        save_matrix=None,
                        superimpose=True,
                        superimposition_subset="name CA",
                        mass_weighted=True,
                        ncores=1,
                        verbose=False,
                        *conf_dist_args,
                        **conf_dist_kwargs):
    """
    Retrieves or calculates the conformational distance (RMSD)
    matrix. The distance matrix is calculated between all the frames of all
    the :class:`~MDAnalysis.core.universe.Universe` objects given as input.
    The order of the matrix elements depends on the order of the coordinates
    of the ensembles and on the order of the input ensembles themselves,
    therefore the order of the input list is significant.

    The distance matrix can either be calculated from input ensembles or
    loaded from an input numpy binary file.

    Please notice that the .npz file does not contain a bidimensional array,
    but a flattened representation that is meant to represent the elements of
    an encore.utils.TriangularMatrix object.


    Parameters
    ----------
    ensemble : Universe
    selection : str
        Atom selection string in the MDAnalysis format. Default is "name CA"
    load_matrix : str, optional
        Load similarity/dissimilarity matrix from numpy binary file instead
        of calculating it (default is None). A filename is required.
    save_matrix : bool, optional
        Save calculated matrix as numpy binary file (default is None). A
        filename is required.
    superimpose : bool, optional
        Whether to superimpose structures before calculating distance
        (default is True).
    superimposition_subset : str, optional
        Group for superimposition using MDAnalysis selection syntax
        (default is CA atoms: "name CA")
    mass_weighted : bool, optional
        calculate a mass-weighted RMSD (default is True). If set to False
        the superimposition will also not be mass-weighted.
    ncores : int, optional
        Maximum number of cores to be used (default is 1)
    verbose : bool, optional
        print progress

    Returns
    -------

    confdistmatrix : encore.utils.TriangularMatrix
        Conformational distance matrix. .
    """

    # Load the matrix if required
    if load_matrix:
        logging.info(
            "        Loading similarity matrix from: {0}".format(load_matrix))
        confdistmatrix = \
            TriangularMatrix(
                size=ensemble.trajectory.timeseries(
                    ensemble.select_atoms(selection),
                    format='fac').shape[0],
                loadfile=load_matrix)
        logging.info("        Done!")
        for key in confdistmatrix.metadata.dtype.names:
            logging.info("        {0} : {1}".format(
                key, str(confdistmatrix.metadata[key][0])))

        # Check matrix size for consistency
        if not confdistmatrix.size == \
                ensemble.trajectory.timeseries(
                    ensemble.select_atoms(selection),
                    format='fac').shape[0]:
            logging.error(
                "ERROR: The size of the loaded matrix and of the ensemble"
                " do not match")
            return None


    # Calculate the matrix
    else:
        logging.info(
            "        Perform pairwise alignment: {0}".format(str(superimpose)))
        logging.info("        Mass-weighted alignment and RMSD: {0}"
            .format(str(mass_weighted)))
        if superimpose:
            logging.info(
                "        Atoms subset for alignment: {0}"
                    .format(superimposition_subset))
        logging.info("    Calculating similarity matrix . . .")

        # Use superimposition subset, if necessary. If the pairwise alignment
        # is not required, it will not be performed anyway.
        confdistmatrix = conformational_distance_matrix(ensemble,
                                                        conf_dist_function=set_rmsd_matrix_elements,
                                                        selection=selection,
                                                        pairwise_align=superimpose,
                                                        mass_weighted=mass_weighted,
                                                        ncores=ncores,
                                                        verbose=verbose)

        logging.info("    Done!")

        if save_matrix:
            confdistmatrix.savez(save_matrix)

    return confdistmatrix
