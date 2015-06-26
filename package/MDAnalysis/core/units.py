# Placeholder because units moved
# Remove this in version 1.0
import warnings

with warnings.catch_warnings():
    warnings.simplefilter('always', DeprecationWarning)
    warnings.warn(("units has moved to MDAnalysis.units "
                   "and will be removed from here in release 1.0"),
                  DeprecationWarning)

from ..units import *
