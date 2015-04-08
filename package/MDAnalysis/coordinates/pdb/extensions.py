# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#


# pdb.extensions
# original file: edPDB.xpdb but only kept content needed for MDAnalysis
"""
Extensions to :mod:`Bio.PDB` --- :mod:`pdb.extensions`
======================================================

:Author:  Oliver Beckstein
:Year:    2009
:License: Biopython

Extension to :mod:`Bio.PDB` to handle large pdb files.

Partly published on http://biopython.org/wiki/Reading_large_PDB_files
and more code at
http://github.com/orbeckst/GromacsWrapper/tree/master/edPDB/

Classes
-------

.. autoclass:: SloppyStructureBuilder
.. autoclass:: SloppyPDBIO

Functions
---------

.. autofunction:: get_structure
.. autofunction:: write_pdb
"""

import Bio.PDB
import Bio.PDB.StructureBuilder

import logging

logger = logging.getLogger('MDAnalysis.pdb.extensions')


class SloppyStructureBuilder(Bio.PDB.StructureBuilder.StructureBuilder):
    """Cope with resSeq < 10,000 limitation by just incrementing internally.

    Solves the follwing problem with :class:`Bio.PDB.StructureBuilder.StructureBuilder`:
    
    Q: What's wrong here??
       Some atoms or residues will be missing in the data structure.
       WARNING: Residue (' ', 8954, ' ') redefined at line 74803.
       PDBConstructionException: Blank altlocs in duplicate residue SOL (' ', 8954, ' ') at line 74803.

    A: resSeq only goes to 9999 --> goes back to 0 (PDB format is not really good here)

    .. warning::  H and W records are probably not handled yet (don't have examples to test)
    """

    def __init__(self, verbose=False):
        Bio.PDB.StructureBuilder.StructureBuilder.__init__(self)
        self.max_resseq = -1
        self.verbose = verbose

    def init_residue(self, resname, field, resseq, icode):
        """
        Initiate a new Residue object.

        Arguments:
        o resname - string, e.g. "ASN"
        o field - hetero flag, "W" for waters, "H" for
            hetero residues, otherwise blanc.
        o resseq - int, sequence identifier
        o icode - string, insertion code
        """
        if field != " ":
            if field == "H":
                # The hetero field consists of H_ + the residue name (e.g. H_FUC)
                field = "H_" + resname
        res_id = (field, resseq, icode)

        if resseq > self.max_resseq:
            self.max_resseq = resseq

        if field == " ":
            fudged_resseq = False
            while self.chain.has_id(res_id) or resseq == 0:
                # There already is a residue with the id (field, resseq, icode).
                # resseq == 0 catches already wrapped residue numbers which do not
                # trigger the has_id() test.
                #
                # Be sloppy and just increment...
                # (This code will not leave gaps in resids... I think)
                #
                # XXX: shouldn't we also do this for hetero atoms and water??
                self.max_resseq += 1
                resseq = self.max_resseq
                res_id = (field, resseq, icode)  # use max_resseq!
                fudged_resseq = True

            if fudged_resseq and self.verbose:
                logger.debug("Residues are wrapping (Residue ('%s', %i, '%s') at line %i)."
                             % (field, resseq, icode, self.line_counter) +
                             ".... assigning new resid %d.\n" % self.max_resseq)
        residue = Bio.PDB.Residue.Residue(res_id, resname, self.segid)
        self.chain.add(residue)
        self.residue = residue


class SloppyPDBIO(Bio.PDB.PDBIO):
    """PDBIO class that can deal with large pdb files as used in MD simulations.

    - resSeq simply wrap and are printed modulo 10,000.
    - atom numbers wrap at 99,999 and are printed modulo 100,000
    """
    # directly copied from PDBIO.py
    # (has to be copied because of the package layout it is not externally accessible)
    _ATOM_FORMAT_STRING = "%s%5i %-4s%c%3s %c%4i%c   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s\n"

    def _get_atom_line(self, atom, hetfield, segid, atom_number, resname,
                       resseq, icode, chain_id, element="  ", charge="  "):
        """
        Returns an ATOM PDB string that is guaranteed to fit into the ATOM format.

        - Resid (resseq) is wrapped (modulo 10,000) to fit into %4i (4I) format
        - Atom number (atom_number) is wrapped (modulo 100,000) to fit into %4i (4I) format
        """
        if hetfield != " ":
            record_type = "HETATM"
        else:
            record_type = "ATOM  "
        name = atom.get_fullname()
        altloc = atom.get_altloc()
        x, y, z = atom.get_coord()
        bfactor = atom.get_bfactor()
        occupancy = atom.get_occupancy()
        args = (
            record_type, atom_number % 100000, name, altloc, resname, chain_id,
            resseq % 10000, icode, x, y, z, occupancy, bfactor, segid, element, charge)
        return self._ATOM_FORMAT_STRING % args


sloppyparser = Bio.PDB.PDBParser(PERMISSIVE=True, structure_builder=SloppyStructureBuilder())


def get_structure(pdbfile, pdbid='system'):
    """Read the *pdbfilename*  and return a Bio.PDB structure.

    This function ignores duplicate atom numbers and resids from the
    file and simply increments them.

    .. Note::

       The call signature is reversed compared to the one of
       :meth:`Bio.PDB.PDBParser.get_structure`.
    """
    return sloppyparser.get_structure(pdbid, pdbfile)


def write_pdb(structure, filename, **kwargs):
    """Write Bio.PDB molecule *structure* to *filename*.

    :Arguments:
       *structure*
         Bio.PDB structure instance
       *filename*
         pdb file
       *selection*
         Bio.PDB.Selection
    """
    selection = kwargs.pop('selection', None)

    io = SloppyPDBIO()  # deals with resSeq > 9999
    io.set_structure(structure)
    io.save(filename, select=selection)
