# Placeholder because util moved
# Remove this in version 1.0
import warnings

with warnings.catch_warnings():
    warnings.simplefilter('always', DeprecationWarning)
    warnings.warn(("util has moved to MDAnalysis.lib.util "
                   "and will be removed from here in release 1.0"),
                  DeprecationWarning)

from ..lib.util import *
