# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2012 Naveen Michaud-Agrawal,
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
Setting up systems for CHARMM
=============================

Helpers for building and setting up simulations for the CHARMM_
simulation program.

.. _CHARMM: http://www.charmm.org

Preparing a structure for CHARMM
--------------------------------

CHARMM works best if each chain (or segment) is provided as an
individual file. Furthermore, a few atom names need to be changed from
the PDB standard to match the CHARMM force field names. This is all
easily accomplished with the :class:`Preprocessor`::

  from MDAnalysis.builder.charmm import Preprocessor
  P = Preprocessor("protein.pdb")
  P.fix_names()
  P.split_chains("protein_charmm.crd")

If the protein contains chains A, B, and C then this will produce
three files "protein_charmm_A.crd", "protein_charmm_B.crd" and
"protein_charmm_C.crd" that should be readable by CHARMM.


Classes
-------

.. autoclass:: Preprocessor
"""

import numpy as np
import logging

from MDAnalysis import NoDataError

logger = logging.getLogger("MDAnalysis.build.charmm")

import MDAnalysis
class Preprocessor(object):
    """Perform actions on input structure to prepare it for CHARMM_.

    Typically, you load the structure from a PDB (or
    :class:`~MDAnalysis.core.AtomGroup.Universe`), apply fixes
    (e.g. :meth:`Preprocessor.fix_names`) and finally write out coordinate
    files, one per segment, with :meth:`Preprocessor.split_chains`.

    .. Note::

       If a :class:`Universe` instance is provided as input then it is
       *modified* through the :meth:`Preprocessor.fix_names` method.

    """

    def __init__(self, *args, **kwargs):
        """Set up universe with :func:`~MDAnalysis.core.AtomGroup.asUniverse`"""
        self.universe = MDAnalysis.asUniverse(*args, **kwargs)
        self._fixes = [fixILE(self.universe),
                       fixHIS(self.universe),
                       fixCterm(self.universe),
                       ]

    def fix_names(self):
        """Apply one or many fixes to the universe.

        Fixes:

        ===========  =====================================================
          name        action
        ===========  =====================================================
          ILE         ILE CD1 --> CD

          HIS         HIS --> HSD [HSE]

          water       HOH --> TIP3, O --> OH2, segid: XWAT;
                      (not implemented)

          cterm       O --> OT1, OXT --> OT2;
                      uses a distance search between C[i] and N[i+1]
                      to determine the breaks. TODO: only fix those at
                      the end of segments/chains and leave internal gaps
                      unfixed.
        ===========  =====================================================

        """
        for fix in self._fixes:
            fix.apply()

    def split_chains(self, basefilename, lowercase=True, **kwargs):
        """Split input structure into segments (typically by chain).

        Filenames are constructed from *basefilename* by inserting
        the segment/chain identifier before the extension. The
        extension determines the output format (unless it is
        explicitly supplied in the *format* keyword )

        *lowercase* = ``True`` instructs the writer to make the output filename
        all lower case (generally better with CHARMM scripts).
        """
        import os.path
        root, ext = os.path.splitext(basefilename)
        if kwargs.get('format', None):
            ext = "." + kwargs.pop('format').lower()
        for s in self.universe.segments:
            filename = root + "_" + s.id + ext
            if lowercase:
                path, basename = os.path.split(filename)
                filename = os.path.join(path, basename.lower())
            s.atoms.write(filename)
            logger.info("Wrote segment %s to %r", s.id, filename)

class Fix(object):
    """Base class for code that modifies the universe/topology"""

    def __init__(self, universe):
        self.universe = universe
        self.name = self.__class__.__name__
    def apply(self, universe):
        """Apply fix to universe"""
        pass
    def rename_atom(self, selection, old, new):
        """Select group defined by *selection* and rename atoms.

        *old* --> *new*

        *selection* can be a string or an :class:`AtomGroup`.

        """
        try:
            if isinstance(selection, basestring):
                group = self.universe.selectAtoms(selection)
            else:
                # let's hope it's an AtomGroup
                group = selection
            a = group.selectAtoms("name %(old)s" % vars())
            a.set_name(new)
            logger.info("%s: %s -> %s for %d residues.",
                        self.name, old, new, a.numberOfResidues())
        except NoDataError:
            logger.warn("%s: no %s atoms found for selection %r",
                        self.name, old, selection)

class fixILE(Fix):
    """ILE CD1 --> CD"""
    def apply(self):
        """Apply fix to universe"""
        self.rename_atom("resname ILE", "CD1", "CD")

class fixHIS(Fix):
    """HIS --> HSD (default) or HSE.

    Hydrogen bonding pattern is not taken into account; this is a pure rename.
    """
    def apply(self, name="HSD"):
        assert name in ("HSD", "HSE")
        try:
            His = self.universe.selectAtoms("resname HIS")
            His.set_resname(name)
            logger.info("fixHIS: HIS -> %s for %d residues.",
                        name, His.numberOfResidues())
        except NoDataError:
            logger.warn("fixHIS: no His residues found")

class fixCterm(Fix):
    """C-terminal oxygens: O --> OT1, OXT --> OT2"""
    #: Two chains are considered broken if the C[i]-N[i+1] distances
    #: is greater than this value (3 A).
    break_distance = 3.0

    def apply(self):
        """Apply fix to universe"""
        try:
            Cterm = self.find_termini()
        except NoDataError:
            logger.error("fixCterm: failed to detect termini, nothing done")
            return
        self.rename_atom(Cterm, "O", "OT1")
        self.rename_atom(Cterm, "OXT", "OT2")

    def find_termini(self):
        """distance search between backbone C[i] and N[i+1]

        Relies on EACH residue having N and C to get the
        sequential ordering right.

        :Returns: ResidueGroup of C-terminal residues
        """

        proteins = self.universe.selectAtoms("protein")

        N = proteins.selectAtoms("backbone and name N")
        C = proteins.selectAtoms("backbone and name C")
        r = N.coordinates()[1:] - C.coordinates()[:-1]
        d = np.sqrt(np.sum(r**2, axis=-1))
        index_broken, = np.where(d > self.break_distance)
        # any leading broken residues and the very last one
        Cterm = (N[index_broken] + N[-1]).residues
        logger.debug("fixCterm.find_termini(): from distance search with cutoff %g A", self.break_distance)
        logger.debug("fixCterm.find_termini(): %r", Cterm)
        return Cterm

# TODO:
# - rename water oxygens:
#        water       HOH --> TIP3, O --> OH2, segid: XWAT
# - special things such as fatty acids
#        PLM: CA --> C10, CB --> C11, ...
