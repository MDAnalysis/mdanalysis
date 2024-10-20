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
