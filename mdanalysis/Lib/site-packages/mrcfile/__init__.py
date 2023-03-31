# Copyright (c) 2016, Science and Technology Facilities Council
# This software is distributed under a BSD licence. See LICENSE.txt.

"""
mrcfile
=======

A pure Python implementation of the MRC2014 file format.

For a full introduction and documentation, see http://mrcfile.readthedocs.io/

Functions
---------

* :func:`new`: Create a new MRC file.
* :func:`open`: Open an MRC file.
* :func:`open_async`: Open an MRC file asynchronously.
* :func:`mmap`: Open a memory-mapped MRC file (fast for large files).
* :func:`new_mmap`: Create a new empty memory-mapped MRC file (fast for large files).
* :func:`validate`: Validate an MRC file

Basic usage
-----------

Examples assume that this package has been imported as ``mrcfile`` and numpy
has been imported as ``np``.

To open an MRC file and read a slice of data:

>>> with mrcfile.open('tests/test_data/EMD-3197.map') as mrc:
...     mrc.data[10,10]
...
array([ 2.58179283,  3.1406002 ,  3.64495397,  3.63812137,  3.61837363,
        4.0115056 ,  3.66981959,  2.07317996,  0.1251585 , -0.87975615,
        0.12517013,  2.07319379,  3.66982722,  4.0115037 ,  3.61837196,
        3.6381247 ,  3.64495087,  3.14059472,  2.58178973,  1.92690361], dtype=float32)

To create a new file with a 2D data array, and change some values:

>>> with mrcfile.new('tmp.mrc') as mrc:
...     mrc.set_data(np.zeros((5, 5), dtype=np.int8))
...     mrc.data[1:4,1:4] = 10
...     mrc.data
...
array([[ 0,  0,  0,  0,  0],
       [ 0, 10, 10, 10,  0],
       [ 0, 10, 10, 10,  0],
       [ 0, 10, 10, 10,  0],
       [ 0,  0,  0,  0,  0]], dtype=int8)

Background
----------

The MRC2014 format was described in the Journal of Structural Biology:
http://dx.doi.org/10.1016/j.jsb.2015.04.002

The format specification is available on the CCP-EM website:
http://www.ccpem.ac.uk/mrc_format/mrc2014.php

"""

# Import Python 3 features for future-proofing
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from .load_functions import new, open, read, write, open_async, mmap, new_mmap
from .validator import validate
from .version import __version__
