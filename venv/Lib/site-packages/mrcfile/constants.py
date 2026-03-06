# Copyright (c) 2016, Science and Technology Facilities Council
# This software is distributed under a BSD licence. See LICENSE.txt.
"""
constants
---------

Constants used by the ``mrcfile.py`` library.

"""

# Import Python 3 features for future-proofing
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

MAP_ID = b'MAP '
MAP_ID_OFFSET_BYTES = 208  # location of 'MAP ' string in an MRC file

IMAGE_STACK_SPACEGROUP = 0
VOLUME_SPACEGROUP = 1
VOLUME_STACK_SPACEGROUP = 401
