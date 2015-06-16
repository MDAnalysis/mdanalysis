# Placeholder because KDTree moved
__all__ = ['KDTree', 'NeighborSearch']

import warnings
warnings.simplefilter('always', DeprecationWarning)
warnings.warn("KDTree has moved to MDAnalysis.lib.KDTree",
              DeprecationWarning)

from .lib.KDTree import KDTree, NeighborSearch
