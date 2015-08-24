# Placeholder because KDTree moved
# Remove this in version 1.0
__all__ = ['NeighborSearch']

import warnings

with warnings.catch_warnings():
    warnings.simplefilter('always', DeprecationWarning)
    warnings.warn(('KDTree has been removed in 0.11. Instead you can use the '
                   'BioPython or scikit-learn implementation directly. The '
                   '"AtomNeighborSearch" class is still available in the '
                   'NeighborSearch module which is moved to MDAnalysis.lib.NeighborSearch.'
                   'This KDTree module will be removed in the 1.0 release.'),
                  DeprecationWarning)

from .lib import NeighborSearch
