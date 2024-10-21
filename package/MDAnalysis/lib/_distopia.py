# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
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

"""Stub module for distopia --- :mod:`MDAnalysis.analysis.distopia`
===================================================================

This module is a stub to provide distopia distance functions to `distances.py`
as a selectable backend.
"""
import warnings
from packaging.version import Version

MIN_DISTOPIA_VERSION = Version("0.3.1")

# check for distopia
try:
    import distopia
except ImportError:
    HAS_DISTOPIA = False
else:
    HAS_DISTOPIA = True

    # check for compatibility: currently needs to be >=0.3.1,

    # some versions of `distopia` don't have a version attribute
    try:
        distopia_version = Version(distopia.__version__)
    except AttributeError:
        warnings.warn("distopia version cannot be determined, assuming 0.0.0")
        distopia_version = Version("0.0.0")
    else:
        if distopia_version < MIN_DISTOPIA_VERSION:
            warnings.warn(
                f"distopia version {distopia_version} is too old; "
                f"need at least {MIN_DISTOPIA_VERSION}, Your installed version of "
                "distopia will NOT be used."
            )
            HAS_DISTOPIA = False


from .c_distances import (
    calc_bond_distance_triclinic as _calc_bond_distance_triclinic_serial,
)
import numpy as np


def calc_bond_distance_ortho(
    coords1, coords2: np.ndarray, box: np.ndarray, results: np.ndarray
) -> None:
    distopia.calc_bonds_ortho(
        coords1, coords2, box[:3], results=results
    )
    # upcast is currently required, change for 3.0, see #3927


def calc_bond_distance(
    coords1: np.ndarray, coords2: np.ndarray, results: np.ndarray
) -> None:
    distopia.calc_bonds_no_box(
        coords1, coords2, results=results
    )
    # upcast is currently required, change for 3.0, see #3927


def calc_bond_distance_triclinic(
    coords1: np.ndarray, coords2: np.ndarray, box: np.ndarray, results: np.ndarray
) -> None:
    distopia.calc_bonds_triclinic(coords1, coords2, box, results=results)


def calc_angle(
    coords1: np.ndarray, coords2: np.ndarray, coords3: np.ndarray, results: np.ndarray
) -> None:
    distopia.calc_angles_no_box(coords1, coords2, coords3, results=results)

def calc_angle_ortho(
    coords1: np.ndarray, coords2: np.ndarray, coords3: np.ndarray, box: np.ndarray, results: np.ndarray
) -> None:
    distopia.calc_angles_ortho(coords1, coords2, coords3, box[:3], results=results)


def calc_angle_triclinic(
    coords1: np.ndarray, coords2: np.ndarray, coords3: np.ndarray, box: np.ndarray, results: np.ndarray
) -> None:

    distopia.calc_angles_triclinic(coords1, coords2, coords3, box, results=results)


def calc_dihedral(
    coords1: np.ndarray, coords2: np.ndarray, coords3: np.ndarray, coords4: np.ndarray, results: np.ndarray
) -> None:
    distopia.calc_dihedrals_no_box(coords1, coords2, coords3, coords4, results=results)

def calc_dihedral_ortho(
    coords1: np.ndarray, coords2: np.ndarray, coords3: np.ndarray, coords4: np.ndarray, box: np.ndarray, results: np.ndarray
) -> None:
    distopia.calc_dihedrals_ortho(coords1, coords2, coords3, coords4, box[:3], results=results)

def calc_dihedral_triclinic(
    coords1: np.ndarray, coords2: np.ndarray, coords3: np.ndarray, coords4: np.ndarray, box: np.ndarray, results: np.ndarray
) -> None:
    distopia.calc_dihedrals_triclinic(coords1, coords2, coords3, coords4, box, results=results)

def calc_distance_array(
    coords1: np.ndarray, coords2: np.ndarray, results: np.ndarray
) -> None:
    distopia.calc_distance_array_no_box(coords1, coords2, results=results)

def calc_distance_array_ortho(
    coords1: np.ndarray, coords2: np.ndarray, box: np.ndarray, results: np.ndarray
) -> None:
    distopia.calc_distance_array_ortho(coords1, coords2, box[:3], results=results)

def calc_distance_array_triclinic(
    coords1: np.ndarray, coords2: np.ndarray, box: np.ndarray, results: np.ndarray
) -> None:
    distopia.calc_distance_array_triclinic(coords1, coords2, box, results=results)

def calc_self_distance_array(
    coords: np.ndarray, results: np.ndarray
) -> None:
    distopia.calc_self_distance_array_no_box(coords, results=results)

def calc_self_distance_array_ortho(
    coords: np.ndarray, box: np.ndarray, results: np.ndarray
) -> None:
    distopia.calc_self_distance_array_ortho(coords, box[:3], results=results)

def calc_self_distance_array_triclinic(
    coords: np.ndarray, box: np.ndarray, results: np.ndarray
) -> None:
    distopia.calc_self_distance_array_triclinic(coords, box, results=results)

