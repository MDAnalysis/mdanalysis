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
#     Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
#     doi:10.1002/jcc.21787
#

"""
Distance analysis --- :mod:`MDAnalysis.distances`
==================================================

This module provides functions to rapidly compute distances between
atoms or groups of atoms.

"""

__all__ = ['distance_array', 'self_distance_array', 'contact_matrix', 'dist']

import numpy
from scipy import sparse
from scipy import weave
from scipy.weave import converters

from MDAnalysis.core.distances import distance_array, self_distance_array


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

        # TODO Jan: this distance matrix will be symmetric, hence some of the iterations could be skipped.
        if box != None:
            contact_matrix_pbc(coord, sparse_contacts, box, cutoff, progress_meter_freq, suppress_progmet )

        # if no PBC
        else:
            contact_matrix_no_pbc(coord, sparse_contacts, cutoff, progress_meter_freq, suppress_progmet )

        return sparse_contacts

def contact_matrix_pbc(coord, sparse_contacts, box, cutoff, progress_meter_freq, suppress_progmet ):
    print box
    box_half = numpy.array([x / 2. for x in box] )
    print box_half

    c_code = """
    #include <math.h>

    int rows = Ncoord[0];

    float cutoff2 = powf(cutoff, 2);

    py::tuple args(2);

    bool b = 1;

    args[1] = b;

    for (int i=0; i < rows; i++) {

        // Print progress meter
        if (i % progress_meter_freq == 0 ) {
            printf("%.1f percent done \\n", (100.0 * i / rows));
        }

        for (int j=0; j < rows; j++) {
            float x = coord(i,0) - coord(j,0);
            float y = coord(i,1) - coord(j,1);
            float z = coord(i,2) - coord(j,2);


            // Handle the periodicity
            if (fabs(x) > box_half(0) ) {
                if (x < 0.0) {x += box(0); }
                else { x -= box(0); }
            }

            if (fabs(y) > box_half(1) ) {
                if (y < 0.0) {y += box(1); }
                else { y -= box(1); }
            }

            if (fabs(z) > box_half(2) ) {
                if (z < 0.0) {z += box(2); }
                else { z -= box(2); }
            }


            float dist = powf(x, 2) + powf(y, 2) + powf(z, 2);

            if (dist != 0.0 && dist < cutoff2) {

                py::tuple idx(2);

                idx[0] = i;

                idx[1] = j;

                args[0] = idx;

                sparse_contacts.mcall("__setitem__",args);
            }
        }
    }
    """

    weave.inline(c_code, ['coord', 'sparse_contacts', 'box', 'box_half', 'cutoff', 'progress_meter_freq', 'suppress_progmet'], type_converters=converters.blitz)

def contact_matrix_no_pbc(coord, sparse_contacts, cutoff, progress_meter_freq, suppress_progmet ):
    """

    Examples of python.weave usage http://github.com/scipy/scipy/tree/master/scipy/weave/examples
    """

    c_code = """
    #include <math.h>

    int rows = Ncoord[0];

    float cutoff2 = powf(cutoff, 2.);

    py::tuple args(2);

    bool b = 1;

    args[1] = b;

    for (int i=0; i < rows; i++) {

        // Print progress meter
        if (i % progress_meter_freq == 0 ) {
            printf("%.1f percent done \\n", (100.0 * i / rows));
        }

        for (int j=0; j < rows; j++) {
            float x = coord(i,0) - coord(j,0);
            float y = coord(i,1) - coord(j,1);
            float z = coord(i,2) - coord(j,2);

            float dist = powf(x, 2.) + powf(y, 2.) + powf(z, 2.);

            if (dist != 0.0 && dist < cutoff2) {

                py::tuple idx(2);

                idx[0] = i;

                idx[1] = j;

                args[0] = idx;

                sparse_contacts.mcall("__setitem__",args);
            }
        }
    }
    """

    weave.inline(c_code, ['coord', 'sparse_contacts', 'cutoff', 'progress_meter_freq', 'suppress_progmet'], type_converters=converters.blitz)


def dist(A, B, offset=0):
    """Return distance between atoms in two atom groups.

    The distance is calculated atom-wise. The residue ids are also
    returned because a typical use case is to look at CA distances
    before and after an alignment. Using the *offset* keyword one can
    also add a constant offset to the resids which facilitates
    comparison with PDB numbering.

    :Arguments:
       *A*, *B*
          :class:`~MDAnalysis.core.AtomGroup.AtomGroup` with the
          same number of atoms

    :Keywords:
       *offset* : integer
          The *offset* is added to *resids_A* and *resids_B* (see
          below) in order to produce PDB numbers. The default is 0.

       *offset* : tuple
          *offset[0]* is added to *resids_A* and *offset[1]* to
          *resids_B*.

    :Returns: NumPy `array([resids_A, resids_B, distances])`

    """
    if A.atoms.numberOfAtoms() != B.atoms.numberOfAtoms():
        raise ValueError("AtomGroups A and B do not have the same number of atoms")
    try:
        off_A, off_B = offset
    except (TypeError, ValueError):
        off_A = off_B = int(offset)
    residues_A = numpy.array(A.resids()) + off_A
    residues_B = numpy.array(B.resids()) + off_B
    r = A.coordinates() - B.coordinates()
    d = numpy.sqrt(numpy.sum(r*r, axis=1))
    return numpy.array([residues_A, residues_B, d])


def between(group, A, B, distance):
  """Return sub group of *group* that is within *distance* of both *A* and *B*.

  *group*, *A*, and *B* must be
  :class:`~MDAnalysis.core.AtomGroup.AtomGroup` instances.  Works best
  if *group* is bigger than either *A* or *B*. This function is not
  aware of periodic boundary conditions.

  Can be used to find bridging waters or molecules in an interface.

  Similar to "*group* and (AROUND *A* *distance* and AROUND *B* *distance*)".

  .. SeeAlso:: Makes use of :mod:`MDAnalysis.KDTree.NeighborSearch`.

  .. versionadded: 0.7.5
  """
  from MDAnalysis.KDTree.NeighborSearch import AtomNeighborSearch
  from MDAnalysis.core.AtomGroup import AtomGroup
  ns_group = AtomNeighborSearch(group)
  resA = set(ns_group.search_list(A, distance))
  resB = set(ns_group.search_list(B, distance))
  return AtomGroup(resB.intersection(resA))
