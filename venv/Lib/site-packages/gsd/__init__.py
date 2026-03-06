# Copyright (c) 2016-2025 The Regents of the University of Michigan
# Part of GSD, released under the BSD 2-Clause License.

"""The GSD main module.

The main package :py:mod:`gsd` is the root package. It holds the submodules
:py:mod:`gsd.fl` and :py:mod:`gsd.hoomd`, but does not import them by default.
You must explicitly import these modules before use::

    import gsd.fl

    import gsd.hoomd
"""

import signal
import sys

from . import version

# Install a SIGTERM handler that gracefully exits, allowing open files to flush
# buffered writes and close. Catch ValueError and pass as there is no way to
# determine if this is the main interpreter running the main thread prior to
# the call.
try:
    signal.signal(signal.SIGTERM, lambda n, f: sys.exit(1))
except ValueError:
    pass
