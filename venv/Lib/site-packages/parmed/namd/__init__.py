"""
This package contains code for reading, writing, and modifying NAMD files 
needed to start and/or restart simulations with NAMD.  The most common formats
are the so-called "namdbin" files heavily utilized by VMD; these are often
given the extensions .coor and .vel (coordinates and velocities, respectively).
"""
__all__ = ['NamdBinCoor', 'NamdBinVel']
__authors__ = 'Brian Radak'
__license__ = 'GPL v.3'


from .namdbinfiles import NamdBinCoor, NamdBinVel
