# $Id$
"""The formats sub module contains code to read coordinates, either
single frames (e.g. the PDB module) or trajectories (such as the DCD
reader). All readers are supposed to dlivere a Reader object that
represents a common API to the other code.
"""
# TODO: define the API here (for the moment, look at DCD.DCDReader)

__all__ = ['DCD', 'PDB', 'XTC']
