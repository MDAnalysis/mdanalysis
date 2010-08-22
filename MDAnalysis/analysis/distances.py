# fast distance analysis
# Part of MDAnalysis http://mdanalysis.googlecode.com
# Copyright (c) 2006-2010 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# Released under the GNU Public Licence, v2

"""
Distance analysis
=================

This module provides functions to rapidly compute distances between
atoms or groups of atoms.

.. autofunction:: distance_array
.. autofunction:: self_distance_array

"""

__all__ = ['distance_array', 'self_distance_array']
from MDAnalysis.core.distances import distance_array, self_distance_array

