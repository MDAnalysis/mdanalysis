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
from .similarity import hes, ces, dres, ces_convergence, dres_convergence

from .clustering.ClusterCollection import ClusterCollection, Cluster
from .clustering.ClusteringMethod import *
from .clustering.cluster import cluster
from .dimensionality_reduction.DimensionalityReductionMethod import *
from .dimensionality_reduction.reduce_dimensionality import (
    reduce_dimensionality,
)
from .confdistmatrix import get_distance_matrix
from .utils import merge_universes

__all__ = ["covariance", "similarity", "confdistmatrix", "clustering"]

import warnings

wmsg = (
    "Deprecation in version 2.8.0\n"
    "MDAnalysis.analysis.encore is deprecated in favour of the "
    "MDAKit mdaencore (https://www.mdanalysis.org/mdaencore/) "
    "and will be removed in MDAnalysis version 3.0.0."
)
warnings.warn(wmsg, category=DeprecationWarning)

from ...due import due, Doi

due.cite(
    Doi("10.1371/journal.pcbi.1004415"),
    description="ENCORE Ensemble Comparison",
    path="MDAnalysis.analysis.encore",
    cite_module=True,
)
