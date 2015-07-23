# Placeholder because parallel moved
# Remove this in version 1.0
import warnings

with warnings.catch_warnings():
    warnings.simplefilter('always', DeprecationWarning)
    warnings.warn(("parallel has moved to MDAnalysis.lib.parallel "
                   "and will be removed from here in release 1.0"),
                  DeprecationWarning)

from ..lib.parallel import *
