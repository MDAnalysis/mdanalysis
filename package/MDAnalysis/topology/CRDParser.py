# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""
CRD topology parser
===================

Read a list of atoms from a CHARMM CARD coordinate file (CRD_)
to build a basic topology.  Reads atom ids (ATOMNO), atom names (TYPES),
resids (RESID), residue numbers (RESNO), residue names (RESNames), segment ids
(SEGID) and tempfactor (Weighting).  Atom element and mass are guessed based on
the name of the atom.

Residues are detected through a change is either resid or resname
while segments are detected according to changes in segid.

.. _CRD: https://www.charmmtutorial.org/index.php/CHARMM:The_Basics


Classes
-------

.. autoclass:: CRDParser
   :members:
   :inherited-members:

"""
from __future__ import absolute_import
from six import raise_from

import numpy as np

from ..lib.util import openany, FORTRANReader
from .base import TopologyReaderBase, change_squash
from . import guessers
from ..core.topology import Topology
from ..core.topologyattrs import (
    Atomids,
    Atomnames,
    Atomtypes,
    Masses,
    Resids,
    Resnames,
    Resnums,
    Segids,
    Tempfactors,
)


class CRDParser(TopologyReaderBase):
    """Parse a CHARMM CARD coordinate file for topology information.

    Reads the following Attributes:
     - Atomids
     - Atomnames
     - Tempfactors
     - Resids
     - Resnames
     - Resnums
     - Segids

    Guesses the following attributes:
     - Atomtypes
     - Masses
    """
    format = 'CRD'

    def parse(self, **kwargs):
        """Create the Topology object

        Returns
        -------
        MDAnalysis Topology object

        Todo
        ----
        Could use the resnum and temp factor better
        """
        extformat = FORTRANReader('2I10,2X,A8,2X,A8,3F20.10,2X,A8,2X,A8,F20.10')
        stdformat = FORTRANReader('2I5,1X,A4,1X,A4,3F10.5,1X,A4,1X,A4,F10.5')

        atomids = []
        atomnames = []
        tempfactors = []
        resids = []
        resnames = []
        resnums = []
        segids = []

        with openany(self.filename) as crd:
            for linenum, line in enumerate(crd):
                # reading header
                if line.split()[0] == '*':
                    continue
                elif line.split()[-1] == 'EXT' and int(line.split()[0]):
                    r = extformat
                    continue
                elif line.split()[0] == line.split()[-1] and line.split()[0] != '*':
                    r = stdformat
                    continue
                # anything else should be an atom
                try:
                    (serial, resnum, resName, name,
                     x, y, z, segid, resid, tempFactor) = r.read(line)
                except Exception:
                    raise_from(ValueError("Check CRD format at line {0}: {1}"
                                     "".format(linenum + 1, line.rstrip())),
                               None)

                atomids.append(serial)
                atomnames.append(name)
                tempfactors.append(tempFactor)
                resids.append(resid)
                resnames.append(resName)
                resnums.append(resnum)
                segids.append(segid)

        # Convert to np arrays
        atomids = np.array(atomids, dtype=np.int32)
        atomnames = np.array(atomnames, dtype=object)
        tempfactors = np.array(tempfactors, dtype=np.float32)
        resids = np.array(resids, dtype=np.int32)
        resnames = np.array(resnames, dtype=object)
        resnums = np.array(resnums, dtype=np.int32)
        segids = np.array(segids, dtype=object)

        # Guess some attributes
        atomtypes = guessers.guess_types(atomnames)
        masses = guessers.guess_masses(atomtypes)

        atom_residx, (res_resids, res_resnames, res_resnums, res_segids) = change_squash(
            (resids, resnames), (resids, resnames, resnums, segids))
        res_segidx, (seg_segids,) = change_squash(
            (res_segids,), (res_segids,))

        top = Topology(len(atomids), len(res_resids), len(seg_segids),
                       attrs=[
                           Atomids(atomids),
                           Atomnames(atomnames),
                           Atomtypes(atomtypes, guessed=True),
                           Masses(masses, guessed=True),
                           Tempfactors(tempfactors),
                           Resids(res_resids),
                           Resnames(res_resnames),
                           Resnums(res_resnums),
                           Segids(seg_segids),
                       ],
                       atom_resindex=atom_residx,
                       residue_segindex=res_segidx)

        return top
