# Placeholder because distances moved
# Remove this in version 1.0
import warnings

with warnings.catch_warnings():
    warnings.simplefilter('always', DeprecationWarning)
    warnings.warn(("distances has moved to MDAnalysis.lib.distances "
                   "and will be removed from here in release 1.0"),
                  DeprecationWarning)

from ..lib.distances import *
