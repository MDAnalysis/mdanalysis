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
Primitive PDB topology parser
=============================

Use a PDB file to build a minimum internal structure representation (list of atoms).

Reads a PDB file line by line and is not fuzzy about numbering.

.. Warning:: Only cares for atoms and their names; neither
             connectivity nor (partial) charges are deduced. Masses
             are guessed and set to 0 if unknown.
"""


from MDAnalysis.topology.core import guess_atom_type, guess_atom_mass, guess_atom_charge
from MDAnalysis.core.distances import distance_array
import MDAnalysis.coordinates.PDB
import numpy as np

class PDBParseError(Exception):
    """Signifies an error during parsing a PDB file."""
    pass

def parse(filename):
    """Parse atom information from PDB file *filename*.

    :Returns: MDAnalysis internal *structure* dict

    .. SeeAlso:: The *structure* dict is defined in
                 :func:`MDAnalysis.topology.PSFParser.parse` and the file is read with
                 :class:`MDAnalysis.coordinates.PDB.PrimitivePDBReader`.
    """
    structure = {}
    primitive_pdb_reader =  MDAnalysis.coordinates.PDB.PrimitivePDBReader(filename)
    __parseatoms_(primitive_pdb_reader, structure)
    # TODO: reconstruct bonds from CONECT or guess from distance search
    #       (e.g. like VMD)
    __parsebonds_(filename, primitive_pdb_reader, structure)
    return structure

def __parseatoms_(pdb, structure):
    from MDAnalysis.core.AtomGroup import Atom
    attr = "_atoms"  # name of the atoms section
    atoms = []       # list of Atom objects

    # translate list of atoms to MDAnalysis Atom.
    for iatom,atom in enumerate(pdb._atoms):
        atomname = atom.name
        atomtype = atom.element or guess_atom_type(atomname)
        resname = atom.resName
        resid = atom.resSeq
        chain = atom.chainID.strip()
        segid = atom.segID.strip() or chain or "SYSTEM"  # no empty segids (or Universe throws IndexError)
        mass = guess_atom_mass(atomname)
        charge = guess_atom_charge(atomname)
        bfactor = atom.tempFactor
        occupancy = atom.occupancy

        atoms.append(Atom(iatom,atomname,atomtype,resname,int(resid),segid,float(mass),float(charge),
                          bfactor=bfactor))

    structure[attr] = atoms


def _guess_bonds(atoms, coords, fudge_factor=0.7, vdwradii=None):
    """
    Bond between two atoms is created, if the two atoms are within R1 * R2 * 0.6
    of each other, where R1 and R2 are the VdW radii of the atoms and 0.6 is an 
    ad-hoc factor. This is false (and the reference provided below is wrong).
    
    Here the bond is created, when sum of the radii multiplied by some fudge_factor
    (0.7 by default) is greater than the distance between the two atoms.
    
    The VMD radii table is taken from GROMACS (/usr/share/gromacs/top/vdwradii.dat)
    
    No check is done after the bonds are guesses to see if Lewis structre is 
    correct. This is wrong and will burn somebody.
    
    The code is also in pure python now, so it's slow. 
    
    Reference: http://www.ks.uiuc.edu/Research/vmd/vmd-1.7/ug/node23.html
    Author: Jan Domanski
    """
    # Taken from GROMACS gromacs/top/vdwradii.dat; in nm
    # FIXME by JD: these should be stored in an asset file, rather than in 
    # source code.
    # FIXME this is not the whole periodic table... (eg halogens are missing)
    if not vdwradii:
      vdwradii = {  "C":     0.15,
                    "F":     0.12,
                    "H":     0.04,
                    "N":     0.110,
                    "O":     0.105,
                    "S":     0.16,}
    
    assert len(atoms) == coords.shape[0]
    
    # compute all-2-all distance array
    # FIXME by JD only the upper rigth triangle of the distance matrix, without 
    # the diagonal is needed
    dist = distance_array(coords, coords)
    bonds = set()
    # FIXME by JD optimize the code below in cython/scipy.weave
    for i in range(dist.shape[0]):
      for j in range(i+1, dist.shape[0]):
        a1, a2 = atoms[i], atoms[j]
        r1, r2 = vdwradii[a1.type], vdwradii[a2.type] 

        # 10 comes from scaling nm to A        
        if not ((r1 + r2) * 10 * fudge_factor) > dist[i,j] : continue
        #print "BOND", ((r1 + r2) * 10 * fudge_factor), dist[i,j]
        bonds.add(frozenset([i+1,j+1]))

    return bonds

def __parsebonds_(filename, primitive_pdb_reader, structure):
  attr = "_bonds"
  bonds = set()
  
  with open(filename , "r") as filename:
    for num,line in enumerate(filename):
      if line[:6] != "CONECT": continue
      bond = line[6:].split()
      atom, atoms = int(bond[0]) , map(int,bond[1:])
      for a in atoms:
          bond = frozenset([atom, a ])
          bonds.add(bond) 
  # FIXME by JD: we could use a BondsGroup class perhaps 
  structure[attr] = bonds