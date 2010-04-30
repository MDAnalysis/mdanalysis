"""
NumPy access to Gromacs trajectories
====================================

The module provides classes to iterate over frames of Gromacs_ xtc or
trr trajectories. Relevant data (positions, velocities, forces, box
dimensions) are made available as numpy_ arrays.

"""

__all__ = ['XTC', 'TRR']

import XTC, TRR


