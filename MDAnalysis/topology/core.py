"""
:mod:`MDAnalysis.topology.core` --- Common functions for topology building
==========================================================================

.. function:: build_segments
.. function:: build_bondlist
.. class:: Bond

.. function:: get_parser_for
.. function:: guess_format

"""
import os.path
import MDAnalysis.topology

def build_segments(atoms):
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
    def __init__(self, a1, a2):
        self.atom1 = a1
        self.atom2 = a2
    def partner(self, atom):
        if atom is self.atom1:
            return self.atom2
        else: return self.atom1
    def length(self):
        bond = self.atom1.pos - self.atom2.pos
        import math
        return math.sqrt((bond[0]**2)+(bond[1]**2)+(bond[2]**2))

def build_bondlists(atoms, bonds):
    for a in atoms:
        a.bonds = []
    for a1, a2 in bonds:
        atom1 = atoms[a1]
        atom2 = atoms[a2]
        b = Bond(atom1, atom2)
        atom1.bonds.append(b)
        atom2.bonds.append(b)

def get_parser_for(filename):
    """Return the appropriate topology parser for *filename*."""
    return MDAnalysis.topology._topology_parsers[guess_format(filename)]

def guess_format(filename):
    """Returns the type of topology file *filename*.

    The current heuristic simply looks at the filename extension but
    more complicated probes could be implemented here or in the
    individual packages (e.g. as static methods).
    """

    # simple extension checking... something more complicated is left
    # for the ambitious
    root, ext = os.path.splitext(filename)
    try:
        if ext.startswith('.'):
            ext = ext[1:]
        ext = ext.lower()
    except:
        raise TypeError("Cannot determine topology type for %r" % filename)
    
    if not ext in MDAnalysis.topology._topology_parsers:
        raise TypeError("Unknown topology extension %r from %r; only %r are implemented in MDAnalysis." % 
                        (ext, filename, MDAnalysis.topology._topology_parsers.keys()))
    return ext
