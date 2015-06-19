# Placeholder because KDTree moved
__all__ = ['KDTree', 'NeighborSearch']

import warnings

with warnings.catch_warnings():
    warnings.simplefilter('always', DeprecationWarning)
    warnings.warn(("KDTree has moved to MDAnalysis.lib.KDTree "
                   "and will be removed from here in release 1.0"),
                  DeprecationWarning)

from .lib.KDTree import KDTree, NeighborSearch
