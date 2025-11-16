"""
MDAnalysis topology tables
==========================

.. deprecated:: 2.8.0
   The :mod:`MDAnalysis.topology.tables` module has been moved to
   :mod:`MDAnalysis.guesser.tables`. This import point will
   be removed in release 3.0.0.

The module contains static lookup tables for atom typing etc. The
tables are dictionaries that are indexed by the element.

.. autodata:: atomelements
.. autodata:: masses
.. autodata:: vdwradii

The original raw data are stored as multi-line strings that are
translated into dictionaries with :func:`kv2dict`. In the future,
these tables might be moved into external data files; see
:func:`kv2dict` for explanation of the file format.

.. autofunction:: kv2dict

The raw tables are stored in the strings

.. autodata:: TABLE_ATOMELEMENTS
.. autodata:: TABLE_MASSES
.. autodata:: TABLE_VDWRADII
"""

import warnings
from MDAnalysis.guesser.tables import (
    kv2dict,
    TABLE_ATOMELEMENTS,
    atomelements,
    elements,
    TABLE_MASSES,
    masses,
    TABLE_VDWRADII,
    vdwradii,
    Z2SYMB,
    SYMB2Z,
    SYBYL2SYMB,
)

wmsg = (
    "Deprecated in version 2.8.0\n"
    "MDAnalysis.topology.tables has been moved to "
    "MDAnalysis.guesser.tables. This import point "
    "will be removed in MDAnalysis version 3.0.0"
)
warnings.warn(wmsg, category=DeprecationWarning)
