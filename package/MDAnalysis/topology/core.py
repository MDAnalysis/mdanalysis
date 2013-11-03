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
import MDAnalysis.core.distances as distances
from MDAnalysis.core.util import norm
import numpy
import MDAnalysis.core.AtomGroup as AtomGroup

def build_segments(atoms):
    """Create all :class:`~MDAnalysis.core.AtomGroup.Segment` instancess from a list of :class:`~MDAnalysis.core.AtomGroup.Atom` instances.

    The function also builds the :class:`~MDAnalysis.core.AtomGroup.Residue`
    instances by tracking residue numbers.

    Updating segments also changes the underlying
    :class:`~MDAnalysis.core.AtomGroup.Atom` instances, which record
    to which residue and segment an atom belongs.

    :Returns: structure dict, which associates a segname with a
              :class:`~MDAnalysis.core.AtomGroup.Segment`

    """
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
                residues.append(AtomGroup.Residue(curr_resname, curr_resnum, resatomlist))
                resatomlist = [a]
                curr_resnum = a.resid
                curr_resname = a.resname
        else:
            # We've come to a new segment
            residues.append(AtomGroup.Residue(curr_resname, curr_resnum, resatomlist))
            struc[curr_segname] = AtomGroup.Segment(curr_segname, residues)
            residues = []
            resatomlist = [a]
            curr_resnum = a.resid
            curr_resname = a.resname
            curr_segname = a.segid
    # Add the last segment
    residues.append(AtomGroup.Residue(curr_resname, curr_resnum, resatomlist))
    struc[curr_segname] = AtomGroup.Segment(curr_segname, residues)
    return struc

def build_residues(atoms):
    """Create a list :class:`~MDAnalysis.core.AtomGroup.Residue` instances from a list of :class:`~MDAnalysis.core.AtomGroup.Atom` instances.

    Updating residues also changes the underlying
    :class:`~MDAnalysis.core.AtomGroup.Atom` instances, which record
    to which residue an atom belongs.

    :Returns: List of :class:`~MDAnalysis.core.AtomGroup.Residue` instances

    .. versionadded:: 0.8
    """
    struc = {}
    residues = []
    resatomlist = []
    curr_resnum = atoms[0].resid
    curr_resname = atoms[0].resname
    for a in atoms:
        if (a.resid == curr_resnum):
            resatomlist.append(a)
        else:
            # New residue
            residues.append(AtomGroup.Residue(curr_resname, curr_resnum, resatomlist))
            resatomlist = [a]
            curr_resnum = a.resid
            curr_resname = a.resname
    # Add the last residue
    residues.append(AtomGroup.Residue(curr_resname, curr_resnum, resatomlist))
    return residues

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
    
    .. warning:: 

       No check is done after the bonds are guesses to see if Lewis
       structure is correct. This is wrong and will burn somebody.
    
    The code is also in pure python now, so it's slow. 
    
    * Reference: http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.1/ug/node26.html
    * Author: Jan Domanski

    .. versionadded:: 0.7.7
    """
    # Taken from GROMACS gromacs/top/vdwradii.dat; in nm
    # FIXME by JD: these should be stored in an asset file, rather than in 
    # source code.
    # FIXME this is not the whole periodic table... (eg halogens are missing)
    if not vdwradii:
        vdwradii = tables.vdwradii
    
    assert len([a for a in atoms if a]) == coords.shape[0]
    
    bonds = set()
    
    for a in atoms:
        if not a.type in vdwradii.keys(): 
            print a.type + " has no defined vdw radius"
            return bonds

    # 1-D vector of the upper-triangle of all-to-all distance matrix    
    dist = distances.self_distance_array(coords)
    N = len(coords)
    
    pairs = list()
    [[pairs.append((i,j)) for j in range(i + 1, N)] for i in range(N)]
       
    for x, (d, (i, j)) in enumerate(zip(dist, pairs)):
        a1, a2 = atoms[i], atoms[j]
        r1, r2 = vdwradii[a1.type], vdwradii[a2.type] 
        # 10 comes from scaling nm to A        
        if not ((r1 + r2) * fudge_factor) > dist[x] : continue
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

class topologyDict(object):
    """
    TopDict is a wrapper around a dictionary, which contains lists of different types of bonds.

    One is needed for each type of topological group (bonds, angles, torsions etc).

    The keys in the dictionary represent a type of bond ie, a C - C - H angle, with the values
    being a list of the indices of each angle of this type.

    Getting and setting types of bonds is done smartly, so a C - C - H angle is considered identical 
    to a H - C - C angle

    Duplicate entries are prevented on initialisation, (but not during setitem/getitem) 
    but might need a method to purge duplicates!
    """
    def __init__(self, members):
        self.dict = {}
        for bonds in members: #loop through all atoms
            for b in bonds: #loop through all bonds this atom has
                pair = (b.atom1.number, b.atom2.number)
                btype = (b.atom1.type, b.atom2.type)
                if not pair in self[btype] or tuple(reversed(pair)) in self[btype]:
                    self[btype] += [pair]

    def __len__(self):
        """Returns the number of types in the topology dictionary"""
        return len(self.dict.keys())
    def types(self):
        """Returns a list of the different types of available"""
        return self.dict.keys()
    def __repr__(self):
        return '<'+self.__class__.__name__+' with '+repr(len(self))+' types>'
    def __getitem__(self, key):
        """Returns list of atoms that match a given type (or reverse of type)"""
        if key in self.dict:
            return self.dict[key]
        elif tuple(reversed(key)) in self.dict:
            return self.dict[tuple(reversed(key))]
        else:
            #raise KeyError(key)
            return [] #allows using += from start
    def __setitem__(self, key, value):
        #Adding items to the dictionary
        if not key in self: #If it doesn't exist at all, add it
            k = key
        else: #If it does exist, is it "forwards" or "backwards"?
            if key in self.dict: #This check doesn't do the reversal trick, so is "exact" match
                k = key
            else:
                k = tuple(reversed(key))
                
            #Check against adding duplicates, ie a bond with atoms (1,2) is the same as with atoms (2,1)
        self.dict[k] = value

  #              self.dict[tuple(reversed(key))] = value
    def __contains__(self, other):
        """Returns boolean on whether a topology group exists within this dictionary"""
        # For topology groups, 1-2-3 is considered the same as 3-2-1
        return other in self.dict.keys() or tuple(reversed(other)) in self.dict.keys()


class topologyGroup(object):
    def __init__(self, atomgroups):
        if len(atomgroups) == 2:
            self.type = "bond"
            self.atom1 = atomgroups[0] #will fail if Atoms not atomGroups are passed currently
            self.atom2 = atomgroups[1]
            if not len(self.atom1) == len(self.atom2):
                raise TypeError("atomGroup lengths mismatched")
            self.len = len(self.atom1)
        else:
            raise TypeError("Only bonds implemented currently")

    def __len__(self):
        """Number of bonds in the topology group"""
        return self.len #checked on initialisation that all atomgroups are same length

    def __getitem__(self,item):
        """Returns a particular bond as an atomGroup"""
        if numpy.dtype(type(item)) == numpy.dtype(int):
            return self.atom1[item]  + self.atom2[item]
        elif type(item) == slice:
            return self.atom1[item] + self.atom2[item]

    def __iter__(self):

    def bondslow(self, pbc=False):
        """Slow version of bond"""
        if not self.type == 'bond':
            return "This topologyGroup is not a bond group!"
        else:
            if not pbc:
                return numpy.array([norm(a) for a in self.atom1.coordinates() - self.atom2.coordinates()])
            else:
                box = self.atom1.dimensions #assuming both atom1 and atom2 are same universe, maybe check this on initialisation?
                if (box[6:9] == 90.).all() and not (box[0:3] == 0).any(): #orthogonal and divide by zero check
                    distances = self.atom1.coordinates() - self.atom2.coordinates()
                    return numpy.array([distances - numpy.rint(distances/box[0:3])*box[0:3]])

    def bond(self, pbc=False, result=None):
        """
        Calculates the distance between all bonds in this topologyGroup
        """
        if not result:
            result = numpy.zeros((self.len,), numpy.float64)
        if not pbc:
            return distances.bond_distance(self.atom1.coordinates(), self.atom2.coordinates(), result=result)
        else:
            return distances.bond_distance(self.atom1.coordinates(), self.atom2.coordinates(), box=self.atom1.dimensions[0:3], result=result)
        
