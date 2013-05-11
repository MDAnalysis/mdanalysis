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
Common functions for topology building --- :mod:`MDAnalysis.topology.core`
==========================================================================

The various topology parsers make use of functions and classes in this
module. They are mostly of use to developers.

.. SeeAlso:: :mod:`MDAnalysis.topology.tables` for some hard-coded atom
   information that is used by functions such as :func:`guess_atom_type` and
   :func:`guess_atom_mass`.

"""
import os.path
import MDAnalysis.topology
import tables
from MDAnalysis.core.distances import distance_array

def build_segments(atoms):
    """Create all :class:`~MDAnalysis.core.AtomGroup.Segment` instancess from a list of :class:`~MDAnalysis.core.AtomGroup.Atom` instances.

    The function also builds the :class:`~MDAnalysis.core.AtomGroup.Residue`
    instances by tracking residue numbers.
    """
    from MDAnalysis.core.AtomGroup import Residue, Segment
    struc = {}
    residues = []
    resatomlist = []
    curr_segname = atoms[0].segid
    curr_resnum = atoms[0].resid
    curr_resname = atoms[0].resname
    for a in atoms:
        if (a.segid == curr_segname):
            if (a.resid == curr_resnum):
                resatomlist.append(a)
            else:
                # New residue
                residues.append(Residue(curr_resname, curr_resnum, resatomlist))
                resatomlist = [a]
                curr_resnum = a.resid
                curr_resname = a.resname
        else:
            # We've come to a new segment
            residues.append(Residue(curr_resname, curr_resnum, resatomlist))
            struc[curr_segname] = Segment(curr_segname, residues)
            residues = []
            resatomlist = [a]
            curr_resnum = a.resid
            curr_resname = a.resname
            curr_segname = a.segid
    # Add the last segment
    residues.append(Residue(curr_resname, curr_resnum, resatomlist))
    struc[curr_segname] = Segment(curr_segname, residues)
    return struc

class Bond(object):
    """A bond between two :class:`~MDAnalysis.core.AtomGroup.Atom` instances."""
    def __init__(self, a1, a2, order=None):
        self.atom1 = a1
        self.atom2 = a2
        a1.bonds.append(self)
        a2.bonds.append(self)
        self.order = order 
        self.is_guessed = False
    def partner(self, atom):
        if atom is self.atom1:
            return self.atom2
        else: return self.atom1
    def set_is_guessed(self, b):
        self.is_guessed = b
    def get_is_guessed(self):
        return self.is_guessed
    def length(self):
        """Length of the bond."""
        bond = self.atom1.pos - self.atom2.pos
        import math
        return math.sqrt((bond[0]**2)+(bond[1]**2)+(bond[2]**2))
    def __repr__(self):
        a1 = self.atom1
        a2 = self.atom2
        return "< Bond between: Atom %d (%s of %s-%d) and Atom %s (%s of %s-%d), length %.2f A >" % \
          (a1.number+1, a1.name, a1.resname, a1.resid, a2.number+1, a2.name, a2.resname, a2.resid, self.length())


def build_bondlists(atoms, bonds):
    """Construct the bond list of each :class:`~MDAnalysis.core.AtomGroup.Atom`.

    The bond list is stored in the attribute
    :attr:`MDAnalysis.core.AtomGroup.Atom.bonds` and consists of a list of
    :class:`Bond` instances.
    """
    for a in atoms:
        a.bonds = []
    for a1, a2 in bonds:
        atom1 = atoms[a1]
        atom2 = atoms[a2]
        b = Bond(atom1, atom2)
        atom1.bonds.append(b)
        atom2.bonds.append(b)

def get_parser_for(filename, permissive=False, bonds=False, format=None):
    """Return the appropriate topology parser for *filename*.

    Automatic detection is disabled when an explicit *format* is
    provided.
    """
    format = guess_format(filename, format=format)
    if permissive and not bonds:
        return MDAnalysis.topology._topology_parsers_permissive[format]
    if permissive and bonds:
        return MDAnalysis.topology._topology_parsers_bonds[format]
    return MDAnalysis.topology._topology_parsers[format]

def guess_format(filename, format=None):
    """Returns the type of topology file *filename*.

    The current heuristic simply looks at the filename extension but
    more complicated probes could be implemented here or in the
    individual packages (e.g. as static methods).

    If *format* is supplied then it overrides the auto detection.
    """

    if format is None:
        # simple extension checking... something more complicated is left
        # for the ambitious
        root, ext = os.path.splitext(filename)
        try:
            if ext.startswith('.'):
                ext = ext[1:]
            format = ext.upper()
        except:
            raise TypeError("Cannot determine topology type for %r" % filename)
    else:
        # internally, formats are all uppercase
        format = str(format).upper()

    # sanity check
    if not format in MDAnalysis.topology._topology_parsers:
        raise TypeError("Unknown topology format %r for %r; only %r are implemented in MDAnalysis." %
                        (format, filename, MDAnalysis.topology._topology_parsers.keys()))
    return format

# following guess_* used by PDB parser

def guess_atom_type(atomname):
    """Guess atom type from the name.

    At the moment, this function simply returns the element, as
    guessed by :func:`guess_atom_element`.

    .. SeeAlso:: :func:`guess_atom_element` and :mod:`MDAnalysis.topology.tables`
    """
    return guess_atom_element(atomname)

def guess_atom_element(atomname):
    """Guess the element of the atom from the name.

    Looks in dict to see if element is found, otherwise it uses the first character in the atomname.
    The table comes from CHARMM and AMBER atom types, where the first character is not sufficient to
    determine the atom type. Some GROMOS ions have also been added.

    .. Warning: The translation table is incomplete. This will probably result
                in some mistakes, but it still better than nothing!

    .. SeeAlso:: :func:`guess_atom_type` and
                 :mod:`MDAnalysis.topology.tables` (where the data are stored)
    """
    try:
        return tables.atomelements[atomname]
    except KeyError:
        if atomname.startswith(('1', '2', '3', '4', '5', '6', '7', '8', '9')):
            # catch 1HH etc
            return atomname[1]
        return atomname[0]

def guess_bonds(atoms, coords, fudge_factor=0.72, vdwradii=None):
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
    
    assert len([a for a in atoms if a]) == coords.shape[0]
    
    bonds = set()
    
    for a in atoms:
        if not a.type in vdwradii.keys(): 
            print a.type + "has no defined vdw radius"
            return bonds()
    
    # compute all-2-all distance array
    # FIXME by JD only the upper rigth triangle of the distance matrix, without 
    # the diagonal is needed
    dist = distance_array(coords, coords)
    
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

def get_atom_mass(element):
    """Return the atomic mass in u for *element*.

    Masses are looked up in :data:`MDAnalysis.topology.tables.masses`.

    .. Warning:: Unknown masses are set to 0.
    """
    try:
        return tables.masses[element]
    except KeyError:
        return 0.000

def guess_atom_mass(atomname):
    """Guess a mass based on the atom name.

    :func:`guess_atom_element` is used to determine the kind of atom.

    .. warning:: Anything not recognized is simply set to 0; if you rely on the
                 masses you might want to double check.
    """
    return get_atom_mass(guess_atom_element(atomname))

def guess_atom_charge(atomname):
    """Guess atom charge from the name.

    .. Warning:: Not implemented; simply returns 0.
    """
    # TODO: do something slightly smarter, at least use name/element
    return 0.0
