"""Homogeneous Transformation Matrices and Quaternions --- :mod:`MDAnalysis.lib.transformations`
==============================================================================================

.. deprecated:: 2.10.0
   This module will be removed in 3.0.0

The ``lib.transformations`` module makes the functions from the
`transformations package`_ available.

.. SeeAlso::

   `transformations package`_ by Christoph Gohlke.

   The function :func:`MDAnalysis.lib.transformations.rotaxis` was
   moved to :func:`MDAnalysis.lib.mdamath.rotaxis`.

.. _`transformations package`:
   https://github.com/cgohlke/transformations?tab=readme-ov-file#homogeneous-transformation-matrices-and-quaternions

"""

from transformations import *

from .mdamath import rotaxis

from .util import deprecate

rotaxis = deprecate(
    rotaxis,
    old_name="MDAnalysis.lib.transformations.rotaxis",
    new_name="MDAnalysis.lib.mdamath.rotaxis",
    release="2.10.0",
    remove="3.0.0",
)
