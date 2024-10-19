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

# check for distopia
try:
    import distopia
except ImportError:
    HAS_DISTOPIA = False
else:
    HAS_DISTOPIA = True

    # check for compatibility: currently needs to be >=0.2.0,<0.3.0 (issue
    # #4740) No distopia.__version__ available so we have to do some probing.
    needed_funcs = ['calc_bonds_no_box_float', 'calc_bonds_ortho_float']
    has_distopia_020 = all([hasattr(distopia, func) for func in needed_funcs])
    if not has_distopia_020:
        warnings.warn("Install 'distopia>=0.2.0,<0.3.0' to be used with this "
                      "release of MDAnalysis. Your installed version of "
                      "distopia >=0.3.0 will NOT be used.",
                      category=RuntimeWarning)
        del distopia
        HAS_DISTOPIA = False


from .c_distances import (
    calc_bond_distance_triclinic as _calc_bond_distance_triclinic_serial,
)
import numpy as np


def calc_bond_distance_ortho(
    coords1, coords2: np.ndarray, box: np.ndarray, results: np.ndarray
) -> None:
    distopia.calc_bonds_ortho_float(
        coords1, coords2, box[:3], results=results
    )
    # upcast is currently required, change for 3.0, see #3927


def calc_bond_distance(
    coords1: np.ndarray, coords2: np.ndarray, results: np.ndarray
) -> None:
    distopia.calc_bonds_no_box_float(
        coords1, coords2, results=results
    )
    # upcast is currently required, change for 3.0, see #3927


def calc_bond_distance_triclinic(
    coords1: np.ndarray, coords2: np.ndarray, box: np.ndarray, results: np.ndarray
) -> None:
    # redirect to serial backend
    warnings.warn(
        "distopia does not support triclinic boxes, using serial backend instead."
    )
    _calc_bond_distance_triclinic_serial(coords1, coords2, box, results)
