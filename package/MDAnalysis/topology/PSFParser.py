# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
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
PSF topology parser
===================

Reads a CHARMM/NAMD/XPLOR PSF_ file to build the system. The topology will
contain atom IDs, segids, residue IDs, residue names, atom names, atom types,
charges and masses. Bonds, angles, dihedrals and impropers are also read from
the file.

It reads both standard and extended ("EXT") PSF formats and can also parse NAMD
space-separated "PSF" file variants.

.. _PSF: http://www.charmm.org/documentation/c35b1/struct.html

Classes
-------

.. autoclass:: PSFParser
   :members:
   :inherited-members:

"""
import logging
import functools
from math import ceil
import numpy as np

from ..lib.util import openany
from .base import TopologyReaderBase, squash_by, change_squash
from ..core.topologyattrs import (
    Atomids,
    Atomnames,
    Atomtypes,
    Masses,
    Charges,
    Resids,
    Resnums,
    Resnames,
    Segids,
    Bonds,
    Angles,
    Dihedrals,
    Impropers
)
from ..core.topology import Topology

logger = logging.getLogger("MDAnalysis.topology.PSF")


class PSFParser(TopologyReaderBase):
    """Read topology information from a CHARMM/NAMD/XPLOR PSF_ file.

    Creates a Topology with the following Attributes:
    - ids
    - names
    - types
    - masses
    - charges
    - resids
    - resnames
    - segids
    - bonds
    - angles
    - dihedrals
    - impropers

    .. _PSF: http://www.charmm.org/documentation/c35b1/struct.html
    """
    format = 'PSF'

    def parse(self, **kwargs):
        """Parse PSF file into Topology

        Returns
        -------
        MDAnalysis *Topology* object
        """
        # Open and check psf validity
        with openany(self.filename) as psffile:
            header = next(psffile)
            if not header.startswith("PSF"):
                err = ("{0} is not valid PSF file (header = {1})"
                       "".format(self.filename, header))
                logger.error(err)
                raise ValueError(err)
            header_flags = header[3:].split()

            if "NAMD" in header_flags:
                self._format = "NAMD"        # NAMD/VMD
            elif "EXT" in header_flags:
                self._format = "EXTENDED"    # CHARMM
            else:
                self._format = "STANDARD"    # CHARMM

            next(psffile)
            title = next(psffile).split()
            if not (title[1] == "!NTITLE"):
                err = "{0} is not a valid PSF file".format(self.filename)
                logger.error(err)
                raise ValueError(err)
            # psfremarks = [psffile.next() for i in range(int(title[0]))]
            for _ in range(int(title[0])):
                next(psffile)
            logger.debug("PSF file {0}: format {1}"
                         "".format(self.filename, self._format))

            # Atoms first and mandatory
            top = self._parse_sec(
                psffile, ('NATOM', 1, 1, self._parseatoms))
            # Then possibly other sections
            sections = (
                #("atoms", ("NATOM", 1, 1, self._parseatoms)),
                (Bonds, ("NBOND", 2, 4, self._parsesection)),
                (Angles, ("NTHETA", 3, 3, self._parsesection)),
                (Dihedrals, ("NPHI", 4, 2, self._parsesection)),
                (Impropers, ("NIMPHI", 4, 2, self._parsesection)),
                #("donors", ("NDON", 2, 4, self._parsesection)),
                #("acceptors", ("NACC", 2, 4, self._parsesection))
            )

            try:
                for attr, info in sections:
                    next(psffile)
                    top.add_TopologyAttr(
                        attr(self._parse_sec(psffile, info)))
            except StopIteration:
                # Reached the end of the file before we expected
                for attr in (Bonds, Angles, Dihedrals, Impropers):
                    if not hasattr(top, attr.attrname):
                        top.add_TopologyAttr(attr([]))

        return top

    def _parse_sec(self, psffile, section_info):
        """Parse a single section of the PSF

        Returns
        -------
        A list of Attributes from this section
        """
        desc, atoms_per, per_line, parsefunc = section_info
        header = next(psffile)
        while header.strip() == "":
            header = next(psffile)
        header = header.split()
        # Get the number
        num = float(header[0])
        sect_type = header[1].strip('!:')
        # Make sure the section type matches the desc
        if not sect_type == desc:
            err = ("Expected section {0} but found {1}"
                   "".format(desc, sect_type))
            logger.error(err)
            raise ValueError(err)
        # Now figure out how many lines to read
        numlines = int(ceil(num/per_line))

        psffile_next = functools.partial(next, psffile)
        return parsefunc(psffile_next, atoms_per, numlines)

    def _parseatoms(self, lines, atoms_per, numlines):
        """Parses atom section in a Charmm PSF file.

        Normal (standard) and extended (EXT) PSF format are
        supported. CHEQ is supported in the sense that CHEQ data is simply
        ignored.


        CHARMM Format from ``source/psffres.src``:

        CHEQ::
          II,LSEGID,LRESID,LRES,TYPE(I),IAC(I),CG(I),AMASS(I),IMOVE(I),ECH(I),EHA(I)

          standard format:
            (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8,2G14.6)
            (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,2G14.6,I8,2G14.6)  XPLOR
          expanded format EXT:
            (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8,2G14.6)
            (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A4,1X,2G14.6,I8,2G14.6) XPLOR

        no CHEQ::
          II,LSEGID,LRESID,LRES,TYPE(I),IAC(I),CG(I),AMASS(I),IMOVE(I)

         standard format:
            (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8)
            (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,2G14.6,I8)  XPLOR
          expanded format EXT:
            (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8)
            (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A4,1X,2G14.6,I8) XPLOR

        NAMD PSF

        space separated, see release notes for VMD 1.9.1, psfplugin at
        http://www.ks.uiuc.edu/Research/vmd/current/devel.html :

        psfplugin: Added more logic to the PSF plugin to determine cases where the
        CHARMM "EXTended" PSF format cannot accomodate long atom types, and we add
        a "NAMD" keyword to the PSF file flags line at the top of the file. Upon
        reading, if we detect the "NAMD" flag there, we know that it is possible
        to parse the file correctly using a simple space-delimited scanf() format
        string, and we use that strategy rather than holding to the inflexible
        column-based fields that are a necessity for compatibility with CHARMM,
        CNS, X-PLOR, and other formats. NAMD and the psfgen plugin already assume
        this sort of space-delimited formatting, but that's because they aren't
        expected to parse the PSF variants associated with the other programs. For
        the VMD PSF plugin, having the "NAMD" tag in the flags line makes it
        absolutely clear that we're dealing with a NAMD-specific file so we can
        take the same approach.

        """
        # how to partition the line into the individual atom components
        atom_parsers = {
            'STANDARD': lambda l:
            (l[:8], l[9:13].strip() or "SYSTEM", l[14:18],
             l[19:23].strip(), l[24:28].strip(),
             l[29:33].strip(), l[34:48], l[48:62]),
            # l[62:70], l[70:84], l[84:98] ignore IMOVE, ECH and EHA,
            'EXTENDED': lambda l:
            (l[:10], l[11:19].strip() or "SYSTEM", l[20:28],
             l[29:37].strip(), l[38:46].strip(),
             l[47:51].strip(), l[52:66], l[66:70]),
            # l[70:78],  l[78:84], l[84:98] ignore IMOVE, ECH and EHA,
            'NAMD': lambda l: l.split()[:8],
        }
        atom_parser = atom_parsers[self._format]
        # once partitioned, assigned each component the correct type
        set_type = lambda x: (int(x[0]) - 1, x[1] or "SYSTEM", int(x[2]), x[3],
                              x[4], x[5], float(x[6]), float(x[7]))

        # Oli: I don't think that this is the correct OUTPUT format:
        #   psf_atom_format = "   %5d %4s %4d %4s %-4s %-4s %10.6f      %7.4f%s\n"
        # It should be rather something like:
        #   psf_ATOM_format = '%(iatom)8d %(segid)4s %(resid)-4d %(resname)4s '+\
        #                     '%(name)-4s %(type)4s %(charge)-14.6f%(mass)-14.4f%(imove)8d\n'

        # source/psfres.src (CHEQ and now can be used for CHEQ EXTended), see comments above
        #   II,LSEGID,LRESID,LRES,TYPE(I),IAC(I),CG(I),AMASS(I),IMOVE(I),ECH(I),EHA(I)
        #  (I8,1X,A4, 1X,A4,  1X,A4,  1X,A4,  1X,I4,  1X,2G14.6,     I8,   2G14.6)
        #   0:8   9:13   14:18   19:23   24:28   29:33   34:48 48:62 62:70 70:84 84:98

        # Allocate arrays
        atomids = np.zeros(numlines, dtype=np.int32)
        segids = np.zeros(numlines, dtype=object)
        resids = np.zeros(numlines, dtype=np.int32)
        resnames = np.zeros(numlines, dtype=object)
        atomnames = np.zeros(numlines, dtype=object)
        atomtypes = np.zeros(numlines, dtype=object)
        charges = np.zeros(numlines, dtype=np.float32)
        masses = np.zeros(numlines, dtype=np.float64)

        for i in range(numlines):
            try:
                line = lines()
            except StopIteration:
                err = f"{self.filename} is not valid PSF file"
                logger.error(err)
                raise ValueError(err) from None
            try:
                vals = set_type(atom_parser(line))
            except ValueError:
                # last ditch attempt: this *might* be a NAMD/VMD
                # space-separated "PSF" file from VMD version < 1.9.1
                atom_parser = atom_parsers['NAMD']
                vals = set_type(atom_parser(line))
                logger.warning("Guessing that this is actually a"
                               " NAMD-type PSF file..."
                               " continuing with fingers crossed!")
                logger.debug("First NAMD-type line: {0}: {1}"
                             "".format(i, line.rstrip()))

            atomids[i] = vals[0]
            segids[i] = vals[1]
            resids[i] = vals[2]
            resnames[i] = vals[3]
            atomnames[i] = vals[4]
            atomtypes[i] = vals[5]
            charges[i] = vals[6]
            masses[i] = vals[7]

        # Atom
        atomids = Atomids(atomids)
        atomnames = Atomnames(atomnames)
        atomtypes = Atomtypes(atomtypes)
        charges = Charges(charges)
        masses = Masses(masses)

        # Residue
        # resids, resnames
        residx, (new_resids, new_resnames, perres_segids) = change_squash(
            (resids, resnames, segids),
            (resids, resnames, segids))
        # transform from atom:Rid to atom:Rix
        residueids = Resids(new_resids)
        residuenums = Resnums(new_resids.copy())
        residuenames = Resnames(new_resnames)

        # Segment
        segidx, perseg_segids = squash_by(perres_segids)[:2]
        segids = Segids(perseg_segids)

        top = Topology(len(atomids), len(new_resids), len(segids),
                       attrs=[atomids, atomnames, atomtypes,
                              charges, masses,
                              residueids, residuenums, residuenames,
                              segids],
                       atom_resindex=residx,
                       residue_segindex=segidx)

        return top

    def _parsesection(self, lines, atoms_per, numlines):
        section = []

        for i in range(numlines):
            # Subtract 1 from each number to ensure zero-indexing for the atoms
            fields = np.int64(lines().split()) - 1
            for j in range(0, len(fields), atoms_per):
                section.append(tuple(fields[j:j+atoms_per]))
        return section
