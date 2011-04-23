# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2011 Naveen Michaud-Agrawal,
#               Elizabeth J. Denning, Oliver Beckstein,
#               and contributors (see website for details)
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
#     N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
#     O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
#     Molecular Dynamics Simulations. J. Comput. Chem. (2011),
#     doi:10.1002/jcc.21787
#

"""
Distance analysis --- :mod:`MDAnalysis.distances`
==================================================

This module provides functions to rapidly compute distances between
atoms or groups of atoms.

"""

__all__ = ['distance_array', 'self_distance_array', 'contact_matrix']
from MDAnalysis.core.distances import distance_array, self_distance_array
from MDAnalysis.analysis.util import progress_meter
from scipy import sparse


def contact_matrix(coord, cutoff=15.0, returntype="numpy", box=None, progress_meter_freq=100, suppress_progmet=False):
    '''
    Calculates a matrix of contacts between a list of coordinates.
    There is a fast, high-memory-usage version for small systems
    (returntype='numpy'), and a slower, low-memory-usage version
    for larger systems (returntype='sparse').

    If box dimensions are passed, then periodic boundary conditions
    are applied.

    Change progress_meter_freq to alter frequency of progress meter
    updates. Or switch suppress_progmet to True to suppress it completely.
    '''
    if returntype=="numpy":
        adj = (distance_array(coord,coord,box=box) < cutoff)
        return adj

    elif returntype=="sparse":
        # Initialize square List of Lists matrix of dimensions equal to number of coordinates passed
        sparse_contacts = sparse.lil_matrix((len(coord),len(coord))  , dtype='bool')
        # if PBC
        if box != None:
            print box
            box_half = box[0:3] / 2.
            print box_half

            for i in range(len(coord)):
                # Print progress every hundred atoms
                # TODO progress_meter will be changed to a class
                if i % progress_meter_freq == 0 and suppress_progmet == False:
                    progress_meter(i , len(coord))

                for j in range(len(coord)):
                    diff = coord[j] - coord[i]
                    if abs(diff[0]) > box_half[0]:
                        if diff[0] < 0.:
                            diff[0] += box[0]
                        else:
                            diff[0] -= box[0]
                    if abs(diff[1]) > box_half[1]:
                        if diff[1] < 0.:
                            diff[1] += box[1]
                        else:
                            diff[1] -= box[1]
                    if abs(diff[2]) > box_half[2]:
                        if diff[2] < 0.:
                            diff[2] += box[2]
                        else:
                            diff[2] -= box[2]

                    dist = ( ( diff[0] ** 2. ) + ( diff[1] ** 2. ) + ( diff[2] ** 2. ) ) ** 0.5
                    if dist != 0.0 and dist < cutoff:
                        sparse_contacts[i,j] = True

            # if no PBC
        else:
            for i in range(len(coord)):
                # Print progress every hundred atoms
                # TODO progress_meter will be changed to a class
                if i % progress_meter_freq == 0 and suppress_progmet == False:
                    progress_meter(i , len(coord))

                for j in range(len(coord)):
                    diff = coord[j] - coord[i]
                    dist = ( ( diff[0] ** 2. ) + ( diff[1] ** 2. ) + ( diff[2] ** 2. ) ) ** 0.5
                    if dist != 0.0 and dist < cutoff:
                        sparse_contacts[i,j] = True

        return sparse_contacts

