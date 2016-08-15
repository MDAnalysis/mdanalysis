# __init__.py
# Copyright (C) 2014 Wouter Boomsma, Matteo Tiberti
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

__all__ = [
    'covariance',
    'similarity',
    'confdistmatrix',
    'clustering'
]

from .similarity import hes, ces, dres, \
    ces_convergence, dres_convergence

from .clustering.ClusterCollection import ClusterCollection, Cluster
from .clustering.ClusteringMethod import *
from .clustering.cluster import cluster
from .dimensionality_reduction.DimensionalityReductionMethod import *
from .dimensionality_reduction.reduce_dimensionality import (
    reduce_dimensionality)
from .confdistmatrix import get_distance_matrix
from utils import merge_universes