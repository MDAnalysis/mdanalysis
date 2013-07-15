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
PSF topology parser
===================

Reads a CHARMM/NAMD/XPLOR PSF_ file to build the system. Currently uses
the list of atoms (including atom types, which can be either integers
or strings, masses and partial charges) and the bond connectivity.

It reads both standard and extended ("EXT") PSF formats and can also parse NAMD
space-separated "PSF" file variants.

.. _PSF: http://www.charmm.org/documentation/c35b1/struct.html
"""

import warnings
from MDAnalysis import FileFormatWarning

import logging
logger = logging.getLogger("MDAnalysis.topology.PSF")


class PSFParseError(Exception):
    """Signifies an error during parsing of a CHARMM PSF file."""
    pass

def parse(filename):
    """Parse CHARMM/NAMD/XPLOR PSF_ file *filename*.

    :Returns: MDAnalysis internal *structure* dict as defined here.
    """
    # Open and check psf validity
    psffile = open(filename,'r')
    next_line = skip_line = psffile.next
    header = next_line()
    if header[:3] != "PSF":
        logger.error("%s is not valid PSF file (header = %r)", psffile.name, header)
        raise PSFParseError("%s is not a valid PSF file" % psffile.name)
    header_flags = header[3:].split()
    if "NAMD" in header_flags:
        format = "NAMD"        # NAMD/VMD
    elif "EXT" in header_flags:
        format = "EXTENDED"    # CHARMM
    else:
        format = "STANDARD"    # CHARMM
    skip_line()
    title = next_line().split()
    if not (title[1] == "!NTITLE"):
        logger.error("%s is not a valid PSF file", psffile.name)
        raise PSFParseError("%s is not a valid PSF file" % psffile.name)
    psfremarks = [next_line() for i in range(int(title[0]))]
    logger.debug("PSF file %r: format %r", psffile.name, format)

    structure = {}

    def parse_sec(section_info, **kwargs):
        desc, atoms_per, per_line, parsefunc, data_struc = section_info
        from math import ceil
        header = next_line()
        while header.strip() == "": header = next_line()
        header = header.split()
        # Get the number
        num = int(header[0])
        sect_type = header[1].strip('!:')
        # Make sure the section type matches the desc
        if not (sect_type == desc):
            logger.error("Expected section %r but found %r", desc, sect_type)
            raise PSFParseError("Expected section {0} but found {1}".format(desc, sect_type))
        # Now figure out how many lines to read
        numlines = int(ceil(float(num)/per_line))
        # Too bad I don't have generator expressions
        #def repeat(func, num):
        #    for i in xrange(num):
        #        yield func()
        #lines = repeat(next_line, numlines)
        parsefunc(next_line, atoms_per, data_struc, structure, numlines, **kwargs)

    sections = [("NATOM", 1, 1, __parseatoms_, "_atoms"),
                ("NBOND", 2, 4, __parsesection_, "_bonds"),
                ("NTHETA", 3, 3, __parsesection_, "_angles"),
                ("NPHI", 4, 2, __parsesection_, "_dihe"),
                ("NIMPHI", 4, 2, __parsesection_, "_impr"),
                ("NDON", 2, 4, __parsesection_,"_donors"),
                ("NACC", 2, 4, __parsesection_,"_acceptors")]

    try:
        for info in sections:
            skip_line()
            parse_sec(info, format=format)
    except StopIteration:
        # Reached the end of the file before we expected
        if not structure.has_key("_atoms"):
            logger.error("The PSF file didn't contain the minimum required section of NATOM")
            raise PSFParseError("The PSF file didn't contain the minimum required section of NATOM")
    # Who cares about the rest
    psffile.close()
    return structure

def __parseatoms_(lines, atoms_per, attr, structure, numlines, **kwargs):
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
    format = kwargs.pop('format', 'STANDARD')

    # how to partition the line into the individual atom components
    atom_parsers = {
        'STANDARD': lambda l:
            (l[:8], l[9:13].strip() or "SYSTEM", l[14:18], l[19:23].strip(), l[24:28].strip(),
             l[29:33].strip(), l[34:48], l[48:62], l[62:70]),   #  l[70:84], l[84:98] ignore ECH and EHA,
        'EXTENDED': lambda l:
            (l[:10], l[11:19].strip() or "SYSTEM", l[20:28], l[29:37].strip(), l[38:46].strip(),
             l[47:51].strip(), l[52:66], l[66:70], l[70:78]),   #  l[78:84], l[84:98] ignore ECH and EHA,
        'NAMD': lambda l: l.split()[:9],  # will fail if SEGID is empty!
        }
    atom_parser = atom_parsers[format]

    atoms = [None,]*numlines
    from MDAnalysis.core.AtomGroup import Atom

    # Oli: I don't think that this is the correct OUTPUT format:
    #   psf_atom_format = "   %5d %4s %4d %4s %-4s %-4s %10.6f      %7.4f%s\n"
    # It should be rather something like:
    #   psf_ATOM_format = '%(iatom)8d %(segid)4s %(resid)-4d %(resname)4s '+\
    #                     '%(name)-4s %(type)4s %(charge)-14.6f%(mass)-14.4f%(imove)8d\n'

    # source/psfres.src (CHEQ and now can be used for CHEQ EXTended), see comments above
    #   II,LSEGID,LRESID,LRES,TYPE(I),IAC(I),CG(I),AMASS(I),IMOVE(I),ECH(I),EHA(I)
    #  (I8,1X,A4, 1X,A4,  1X,A4,  1X,A4,  1X,I4,  1X,2G14.6,     I8,   2G14.6)
    #   0:8   9:13   14:18   19:23   24:28   29:33   34:48 48:62 62:70 70:84 84:98


    for i in xrange(numlines):
        line = lines()
        try:
            iatom, segid, resid, resname, atomname, atomtype, charge, mass, imove = atom_parser(line)
            # Atom(atomno, atomname, type, resname, resid, segid, mass, charge)
            # We want zero-indexing for atom numbers to make it easy
            atom_desc = Atom(int(iatom)-1,atomname,atomtype,resname,int(resid),segid,float(mass),float(charge))
        except ValueError:
            # last ditch attempt: this *might* be a NAMD/VMD space-separated "PSF" file from
            # VMD version < 1.9.1
            atom_parser = atom_parsers['NAMD']
            iatom, segid, resid, resname, atomname, atomtype, charge, mass, imove = atom_parser(line)
            atom_desc = Atom(int(iatom)-1,atomname,atomtype,resname,int(resid),segid,float(mass),float(charge))
            warnings.warn("Guessing that this is actually a NAMD-type PSF file... continuing with fingers crossed!",
                          category=FileFormatWarning)
            logger.warn("Guessing that this is actually a NAMD-type PSF file... continuing with fingers crossed!")
            logger.debug("First NAMD-type line: %d: %r", i, line.rstrip())

        atoms[i] = atom_desc

    structure[attr] = atoms

def __parsesection_(lines, atoms_per, attr, structure, numlines, **kwargs):
    section = [] #[None,]*numlines
    #for l in lines:
    for i in xrange(numlines):
        l = lines()
        # Subtract 1 from each number to ensure zero-indexing for the atoms
        f = map(int, l.split())
        fields = [a-1 for a in f]
        for j in range(0, len(fields), atoms_per):
            section.append(tuple(fields[j:j+atoms_per]))
    structure[attr] = section

