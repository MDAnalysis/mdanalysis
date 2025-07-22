# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
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
"""
Distance Matrix calculation --- :mod:`MDAnalysis.analysis.ensemble.confdistmatrix`
==================================================================================


The module contains a base class to easily compute, using
parallelization and shared memory, matrices of conformational
distance between the structures stored as frames in a Universe. A
class to compute an RMSD matrix in such a way is also available.

:Author: Matteo Tiberti, Wouter Boomsma, Tone Bengtsen

.. versionadded:: 0.16.0

.. deprecated:: 2.8.0
   This module is deprecated in favour of the 
   MDAKit `mdaencore <https://mdanalysis.org/mdaencore/>`_ and will be removed
   in MDAnalysis 3.0.0.

"""
from joblib import Parallel, delayed
import numpy as np
from getpass import getuser
from socket import gethostname
from datetime import datetime
from time import sleep
import logging
import warnings

from ...core.universe import Universe

from ..align import rotation_matrix

from .cutils import PureRMSD
from .utils import TriangularMatrix, trm_indices


def conformational_distance_matrix(
    ensemble,
    conf_dist_function,
    select="",
    superimposition_select="",
    n_jobs=1,
    pairwise_align=True,
    weights="mass",
    metadata=True,
    verbose=False,
    max_nbytes=None,
):
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
    select : str, optional
        use this selection for the calculation of conformational distance
    superimposition_select : str, optional
        use atoms from this selection for fitting instead of those of
        `select`
    pairwise_align : bool, optional
        Whether to perform pairwise alignment between conformations.
        Default is True (do the superimposition)
    weights : str/array_like, optional
       weights to be used for fit. Can be either 'mass' or an array_like
    metadata : bool, optional
        Whether to build a metadata dataset for the calculated matrix.
        Default is True.
    n_jobs : int, optional
        Number of cores to be used for parallel calculation
        Default is 1. -1 uses all available cores
    max_nbytes : str, optional
        Threshold on the size of arrays passed to the workers that triggers automated memory mapping in temp_folder (default is None).
        See https://joblib.readthedocs.io/en/latest/generated/joblib.Parallel.html for detailed documentation.
    verbose : bool, optional
        enable verbose output

    Returns
    -------
    conf_dist_matrix : encore.utils.TriangularMatrix object
        Conformational distance matrix in triangular representation.

    """

    # framesn: number of frames
    framesn = len(
        ensemble.trajectory.timeseries(
            ensemble.select_atoms(select), order="fac"
        )
    )

    # Prepare metadata recarray
    if metadata:
        metadata = np.array(
            [
                (
                    gethostname(),
                    getuser(),
                    str(datetime.now()),
                    ensemble.filename,
                    framesn,
                    pairwise_align,
                    select,
                    weights == "mass",
                )
            ],
            dtype=[
                ("host", object),
                ("user", object),
                ("date", object),
                ("topology file", object),
                ("number of frames", int),
                ("pairwise superimposition", bool),
                ("superimposition subset", object),
                ("mass-weighted", bool),
            ],
        )

    # Prepare alignment subset coordinates as necessary

    rmsd_coordinates = ensemble.trajectory.timeseries(
        ensemble.select_atoms(select), order="fac"
    )

    if pairwise_align:
        if superimposition_select:
            subset_select = superimposition_select
        else:
            subset_select = select

        fitting_coordinates = ensemble.trajectory.timeseries(
            ensemble.select_atoms(subset_select), order="fac"
        )
    else:
        fitting_coordinates = None

    if (
        not isinstance(weights, (list, tuple, np.ndarray))
        and weights == "mass"
    ):
        weights = ensemble.select_atoms(select).masses.astype(np.float64)
        if pairwise_align:
            subset_weights = ensemble.select_atoms(
                subset_select
            ).masses.astype(np.float64)
        else:
            subset_weights = None
    elif weights is None:
        weights = np.ones(
            (
                ensemble.trajectory.timeseries(ensemble.select_atoms(select))[
                    0
                ].shape[0]
            )
        ).astype(np.float64)
        if pairwise_align:
            subset_weights = np.ones((fit_coords[0].shape[0])).astype(
                np.float64
            )
        else:
            subset_weights = None
    else:
        if pairwise_align:
            if len(weights) != 2:
                raise RuntimeError(
                    "used pairwise alignment with custom "
                    "weights. Please provide 2 tuple with "
                    "weights for 'select' and "
                    "'superimposition_select'"
                )
            subset_weights = weights[1]
            weights = weights[0]
        else:
            subset_weights = None

    # Allocate for output matrix
    matsize = framesn * (framesn + 1) // 2
    distmat = np.empty(matsize, np.float64)

    # Initialize workers. Simple worker doesn't perform fitting,
    # fitter worker does.
    indices = trm_indices((0, 0), (framesn - 1, framesn - 1))
    Parallel(
        n_jobs=n_jobs,
        verbose=verbose,
        require="sharedmem",
        max_nbytes=max_nbytes,
    )(
        delayed(conf_dist_function)(
            np.int64(element),
            rmsd_coordinates,
            distmat,
            weights,
            fitting_coordinates,
            subset_weights,
        )
        for element in indices
    )

    # When the workers have finished, return a TriangularMatrix object
    return TriangularMatrix(distmat, metadata=metadata)


def set_rmsd_matrix_elements(
    tasks,
    coords,
    rmsdmat,
    weights,
    fit_coords=None,
    fit_weights=None,
    *args,
    **kwargs,
):
    """
    RMSD Matrix calculator

    Parameters
    ----------
    tasks : iterator of int of length 2
        Given a triangular matrix, this function will calculate RMSD
        values from element tasks[0] to tasks[1]. Since the matrix
        is triangular, the trm_indices matrix automatically
        calculates the corrisponding i,j matrix indices.
        The matrix is written as an array in a row-major
        order (see the TriangularMatrix class for details).

        If fit_coords and fit_weights are specified, the structures
        will be superimposed before calculating RMSD, and fit_coords and fit_weights
        will be used to place both structures at their center of mass and
        compute the rotation matrix. In this case, both fit_coords and fit_weights
        must be specified.
    coords : numpy.array
        Array of the ensemble coordinates
    weights : numpy.array
        Array of atomic weights, having the same order as the
        coordinates array
    rmsdmat : encore.utils.TriangularMatrix
        Memory-shared triangular matrix object
    fit_coords : numpy.array or None, optional
        Array of the coordinates used for fitting
    fit_weights : numpy.array. optional
        Array of atomic weights, having the same order as the
        fit_coords array
    """
    i, j = tasks

    if fit_coords is None and fit_weights is None:
        sumweights = np.sum(weights)
        rmsdmat[(i + 1) * i // 2 + j] = PureRMSD(
            coords[i].astype(np.float64),
            coords[j].astype(np.float64),
            coords[j].shape[0],
            weights,
            sumweights,
        )

    elif fit_coords is not None and fit_weights is not None:
        sumweights = np.sum(weights)
        subset_weights = np.asarray(fit_weights) / np.mean(fit_weights)
        com_i = np.average(fit_coords[i], axis=0, weights=fit_weights)
        translated_i = coords[i] - com_i
        subset1_coords = fit_coords[i] - com_i
        com_j = np.average(fit_coords[j], axis=0, weights=fit_weights)
        translated_j = coords[j] - com_j
        subset2_coords = fit_coords[j] - com_j
        rotamat = rotation_matrix(
            subset1_coords, subset2_coords, subset_weights
        )[0]
        rotated_i = np.transpose(np.dot(rotamat, np.transpose(translated_i)))
        rmsdmat[(i + 1) * i // 2 + j] = PureRMSD(
            rotated_i.astype(np.float64),
            translated_j.astype(np.float64),
            coords[j].shape[0],
            weights,
            sumweights,
        )
    else:
        raise TypeError(
            "Both fit_coords and fit_weights must be specified "
            "if one of them is given"
        )


def get_distance_matrix(
    ensemble,
    select="name CA",
    load_matrix=None,
    save_matrix=None,
    superimpose=True,
    superimposition_subset="name CA",
    weights="mass",
    n_jobs=1,
    max_nbytes=None,
    verbose=False,
    *conf_dist_args,
    **conf_dist_kwargs,
):
    """
    Retrieves or calculates the conformational distance (RMSD)
    matrix. The distance matrix is calculated between all the frames of all
    the :class:`~MDAnalysis.core.universe.Universe` objects given as input.
    The order of the matrix elements depends on the order of the coordinates
    of the ensembles and on the order of the input ensembles themselves,
    therefore the order of the input list is significant.

    The distance matrix can either be calculated from input ensembles or
    loaded from an input numpy binary file.

    Please notice that the .npz file does not contain a bi-dimensional array,
    but a flattened representation that is meant to represent the elements of
    an encore.utils.TriangularMatrix object.


    Parameters
    ----------
    ensemble : Universe
    select : str
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
    weights : str/array_like, optional
        weights to be used for fit. Can be either 'mass' or an array_like
    n_jobs : int, optional
        Maximum number of cores to be used (default is 1). If -1 use all cores.
    max_nbytes : str, optional
        Threshold on the size of arrays passed to the workers that triggers automated memory mapping in temp_folder (default is None).
        See https://joblib.readthedocs.io/en/latest/generated/joblib.Parallel.html for detailed documentation.
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
            "        Loading similarity matrix from: {0}".format(load_matrix)
        )
        confdistmatrix = TriangularMatrix(
            size=ensemble.trajectory.timeseries(
                ensemble.select_atoms(select), order="fac"
            ).shape[0],
            loadfile=load_matrix,
        )
        logging.info("        Done!")
        for key in confdistmatrix.metadata.dtype.names:
            logging.info(
                "        {0} : {1}".format(
                    key, str(confdistmatrix.metadata[key][0])
                )
            )

        # Check matrix size for consistency
        if (
            not confdistmatrix.size
            == ensemble.trajectory.timeseries(
                ensemble.select_atoms(select), order="fac"
            ).shape[0]
        ):
            logging.error(
                "ERROR: The size of the loaded matrix and of the ensemble"
                " do not match"
            )
            return None

    # Calculate the matrix
    else:

        # Transfer universe to memory to ensure timeseries() support
        ensemble.transfer_to_memory()

        if (
            not isinstance(weights, (list, tuple, np.ndarray))
            and weights == "mass"
        ):
            weight_type = "Mass"
        elif weights is None:
            weight_type = "None"
        else:
            weight_type = "Custom"
        logging.info(
            "        Perform pairwise alignment: {0}".format(str(superimpose))
        )
        logging.info(
            "        weighted alignment and RMSD: {0}".format(weight_type)
        )
        if superimpose:
            logging.info(
                "        Atoms subset for alignment: {0}".format(
                    superimposition_subset
                )
            )
        logging.info("    Calculating similarity matrix . . .")

        # Use superimposition subset, if necessary. If the pairwise alignment
        # is not required, it will not be performed anyway.
        confdistmatrix = conformational_distance_matrix(
            ensemble,
            conf_dist_function=set_rmsd_matrix_elements,
            select=select,
            pairwise_align=superimpose,
            weights=weights,
            n_jobs=n_jobs,
            max_nbytes=max_nbytes,
            verbose=verbose,
        )

        logging.info("    Done!")

        if save_matrix:
            confdistmatrix.savez(save_matrix)

    return confdistmatrix
