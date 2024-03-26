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
PDB Topology Parser
=========================================================================

This topology parser uses a standard PDB file to build a minimum
internal structure representation (list of atoms).

The topology reader reads a PDB file line by line and ignores atom
numbers but only reads residue numbers up to 9,999 correctly. If you
have systems containing at least 10,000 residues then you need to use
a different file format (e.g. the "extended" PDB, *XPDB* format, see
:mod:`~MDAnalysis.topology.ExtendedPDBParser`) that can handle residue
numbers up to 99,999.



.. Note::

   The parser processes atoms and their names.
   Partial charges are not set. Elements are parsed if they are
   valid. If partially missing or incorrect, empty records are assigned.

See Also
--------
* :mod:`MDAnalysis.topology.ExtendedPDBParser`
* :class:`MDAnalysis.coordinates.PDB.PDBReader`
* :class:`MDAnalysis.core.universe.Universe`


Classes
-------

.. autoclass:: PDBParser
   :members:
   :inherited-members:

"""
import numpy as np
import warnings

from ..guesser.tables import SYMB2Z
from ..lib import util
from .base import TopologyReaderBase, change_squash
from ..core.topology import Topology
from ..core.topologyattrs import (
    Atomnames,
    Atomids,
    AltLocs,
    Bonds,
    ChainIDs,
    Atomtypes,
    Elements,
    ICodes,
    Occupancies,
    RecordTypes,
    Resids,
    Resnames,
    Resnums,
    Segids,
    Tempfactors,
    FormalCharges,
)


def float_or_default(val, default):
    try:
        return float(val)
    except ValueError:
        return default


DIGITS_UPPER = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
DIGITS_LOWER = DIGITS_UPPER.lower()
DIGITS_UPPER_VALUES = dict([pair for pair in zip(DIGITS_UPPER, range(36))])
DIGITS_LOWER_VALUES = dict([pair for pair in zip(DIGITS_LOWER, range(36))])


def decode_pure(digits_values, s):
    """Decodes the string s using the digit, value associations for each
    character.

    Parameters
    ----------
    digits_values: dict
        A dictionary containing the base-10 numbers that each hexadecimal
        number corresponds to.
    s: str
        The contents of the pdb index columns.

    Returns
    -------
    The integer in base-10 corresponding to traditional base-36.
    """
    result = 0
    n = len(digits_values)
    for c in s:
        result *= n
        result += digits_values[c]
    return result


def hy36decode(width, s):
    """
    Decodes base-10/upper-case base-36/lower-case base-36 hybrid.

    Parameters
    ----------
    width: int
        The number of columns in the pdb file store atom index.
    s: str
        The contents of the pdb index columns.

    Returns
    -------
    int
        Base-10 integer corresponding to hybrid36.
    """
    if (len(s) == width):
        f = s[0]
        if (f == "-" or f == " " or f.isdigit()):
            return int(s)
        elif (f in DIGITS_UPPER_VALUES):
            return decode_pure(digits_values=DIGITS_UPPER_VALUES,
                               s=s) - 10 * 36 ** (width - 1) + 10 ** width
        elif (f in DIGITS_LOWER_VALUES):
            return decode_pure(digits_values=DIGITS_LOWER_VALUES,
                               s=s) + 16 * 36 ** (width - 1) + 10 ** width
    raise ValueError("invalid number literal.")


class PDBParser(TopologyReaderBase):
    """Parser that obtains a list of atoms from a standard PDB file.

    Creates the following Attributes (if present):
     - names
     - chainids
     - tempfactors
     - occupancies
     - record_types (ATOM/HETATM)
     - resids
     - resnames
     - segids
     - elements
     - bonds
     - formalcharges

    See Also
    --------
    :class:`MDAnalysis.coordinates.PDB.PDBReader`


    .. versionadded:: 0.8
    .. versionchanged:: 0.18.0
       Added parsing of Record types
    .. versionchanged:: 1.0.0
       Added parsing of valid Elements
    .. versionchanged:: 2.0.0
       Bonds attribute is not added if no bonds are present in PDB file.
       If elements are invalid or partially missing, empty elements records
       are now assigned (Issue #2422).
       Aliased ``bfactors`` topologyattribute to ``tempfactors``.
       ``bfactors`` is deprecated and will be removed in 3.0 (Issue #1901)
    .. versionchanged:: 2.3.0
       Formal charges are now read from PDB files if present. No formalcharge
       attribute is created if no formal charges are present in the PDB file.
       Any formal charges not set are assumed to have a value of 0.
       Raise `UserWarning` instead `RuntimeError`
       when CONECT records are corrupt.
    .. versionchanged:: 2.5.0
       Formal charges will not be populated if an unknown entry is encountered,
       instead a UserWarning is emitted.
    .. versionchanged:: 2.8.0
        Removed type and mass guessing (attributes guessing takes place now
        through universe.guess_TopologyAttrs() API).
    """
    format = ['PDB', 'ENT']

    def parse(self, **kwargs):
        """Parse atom information from PDB file

        Returns
        -------
        MDAnalysis Topology object
        """
        top = self._parseatoms()

        try:
            bonds = self._parsebonds(top.ids.values)
        except AttributeError:
            warnings.warn("Invalid atom serials were present, "
                          "bonds will not be parsed")
        except RuntimeError:
            warnings.warn("CONECT records was corrupt, "
                          "bonds will not be parsed")
        else:
            # Issue 2832: don't append Bonds if there are no bonds
            if bonds:
                top.add_TopologyAttr(bonds)

        return top

    def _parseatoms(self):
        """Create the initial Topology object"""
        resid_prev = 0  # resid looping hack

        record_types = []
        serials = []
        names = []
        altlocs = []
        chainids = []
        icodes = []
        tempfactors = []
        occupancies = []
        resids = []
        resnames = []
        segids = []
        elements = []
        formalcharges = []

        self._wrapped_serials = False  # did serials go over 100k?
        last_wrapped_serial = 100000  # if serials wrap, start from here
        with util.openany(self.filename) as f:
            for line in f:
                line = line.strip()  # Remove extra spaces
                if not line:  # Skip line if empty
                    continue
                if line.startswith('END'):
                    break
                if not line.startswith(('ATOM', 'HETATM')):
                    continue

                record_types.append(line[:6].strip())
                try:
                    serial = int(line[6:11])
                except:
                    try:
                        serial = hy36decode(5, line[6:11])
                    except ValueError:
                        # serial can become '***' when they get too high
                        self._wrapped_serials = True
                        serial = last_wrapped_serial
                        last_wrapped_serial += 1
                finally:
                    serials.append(serial)

                names.append(line[12:16].strip())
                altlocs.append(line[16:17].strip())
                resnames.append(line[17:21].strip())
                chainids.append(line[21:22].strip())
                elements.append(line[76:78].strip())
                formalcharges.append(line[78:80].strip())

                # Resids are optional
                try:
                    if self.format == "XPDB":  # fugly but keeps code DRY
                        # extended non-standard format used by VMD
                        resid = int(line[22:27])
                    else:
                        resid = int(line[22:26])
                        # Wrapping
                        while resid - resid_prev < -5000:
                            resid += 10000
                        resid_prev = resid
                except ValueError:
                    warnings.warn("PDB file is missing resid information.  "
                                  "Defaulted to '1'")
                    resid = 1
                finally:
                    resids.append(resid)
                    icodes.append(line[26:27].strip())

                occupancies.append(float_or_default(line[54:60], 0.0))
                tempfactors.append(float_or_default(line[60:66], 1.0))  # AKA bfactor

                segids.append(line[66:76].strip())

        # Warn about wrapped serials
        if self._wrapped_serials:
            warnings.warn("Serial numbers went over 100,000.  "
                          "Higher serials have been guessed")

        # If segids not present, try to use chainids
        if not any(segids):
            segids = chainids

        n_atoms = len(serials)

        attrs = []
        # Make Atom TopologyAttrs
        for vals, Attr, dtype in (
                (names, Atomnames, object),
                (altlocs, AltLocs, object),
                (chainids, ChainIDs, object),
                (record_types, RecordTypes, object),
                (serials, Atomids, np.int32),
                (tempfactors, Tempfactors, np.float32),
                (occupancies, Occupancies, np.float32),
        ):
            attrs.append(Attr(np.array(vals, dtype=dtype)))
        # OPT: We do this check twice, maybe could refactor to avoid this
        if any(elements):
            # Feed atomtypes as raw element column, but validate elements
            atomtypes = elements
            attrs.append(Atomtypes(np.array(elements, dtype=object)))

            validated_elements = []
            for elem in elements:
                if elem.capitalize() in SYMB2Z:
                    validated_elements.append(elem.capitalize())
                else:
                    wmsg = (f"Unknown element {elem} found for some atoms. "
                            f"These have been given an empty element record. "
                            f"If needed they can be guessed using "
                            f"universe.guess_TopologyAttrs(context='default',"
                            " to_guess=['elements']).")
                    warnings.warn(wmsg)
                    validated_elements.append('')
            attrs.append(Elements(np.array(validated_elements, dtype=object)))
        else:
            warnings.warn("Element information is missing, elements attribute "
                          "will not be populated. If needed these can be"
                          " guessed using universe.guess_TopologyAttrs("
                          "context='default', to_guess=['elements']).")

        if any(formalcharges):
            try:
                for i, entry in enumerate(formalcharges):
                    if not entry == '':
                        if entry == '0':
                            # Technically a lack of charge shouldn't be in the
                            # PDB but MDA has a few files that specifically
                            # have 0 entries, indicating that some folks
                            # interpret 0 as an allowed entry
                            formalcharges[i] = 0
                        elif ('+' in entry) or ('-' in entry):
                            formalcharges[i] = int(entry[::-1])
                        else:
                            raise ValueError
                    else:
                        formalcharges[i] = 0
            except ValueError:
                wmsg = (f"Unknown entry {entry} encountered in formal charge "
                        "field. This likely indicates that the PDB file is "
                        "not fully standard compliant. The formalcharges "
                        "attribute will not be populated.")
                warnings.warn(wmsg)
            else:
                attrs.append(FormalCharges(np.array(formalcharges, dtype=int)))

        # Residue level stuff from here
        resids = np.array(resids, dtype=np.int32)
        resnames = np.array(resnames, dtype=object)
        if self.format == 'XPDB':  # XPDB doesn't have icodes
            icodes = [''] * n_atoms
        icodes = np.array(icodes, dtype=object)
        resnums = resids.copy()
        segids = np.array(segids, dtype=object)

        residx, (resids, resnames, icodes, resnums, segids) = change_squash(
            (resids, resnames, icodes, segids), (resids, resnames, icodes, resnums, segids))
        n_residues = len(resids)
        attrs.append(Resnums(resnums))
        attrs.append(Resids(resids))
        attrs.append(Resnums(resids.copy()))
        attrs.append(ICodes(icodes))
        attrs.append(Resnames(resnames))

        if any(segids) and not any(val is None for val in segids):
            segidx, (segids,) = change_squash((segids,), (segids,))
            n_segments = len(segids)
            attrs.append(Segids(segids))
        else:
            n_segments = 1
            attrs.append(Segids(np.array(['SYSTEM'], dtype=object)))
            segidx = None

        top = Topology(n_atoms, n_residues, n_segments,
                       attrs=attrs,
                       atom_resindex=residx,
                       residue_segindex=segidx)

        return top

    def _parsebonds(self, serials):
        # Could optimise this by saving lines in the main loop
        # then doing post processing after all Atoms have been read
        # ie do one pass through the file only
        # Problem is that in multiframe PDB, the CONECT is at end of file,
        # so the "break" call happens before bonds are reached.

        # If the serials wrapped, this won't work
        if self._wrapped_serials:
            warnings.warn("Invalid atom serials were present, bonds will not"
                          " be parsed")
            raise AttributeError  # gets caught in parse

        # Mapping between the atom array indicies a.index and atom ids
        # (serial) in the original PDB file
        mapping = dict((s, i) for i, s in enumerate(serials))

        bonds = set()
        with util.openany(self.filename) as f:
            lines = (line for line in f if line[:6] == "CONECT")
            for line in lines:
                atom, atoms = _parse_conect(line.strip())
                for a in atoms:
                    try:
                        bond = tuple([mapping[atom], mapping[a]])
                    except KeyError:
                        # Bonds to TER records have no mapping
                        # Ignore these as they are not real atoms
                        warnings.warn(
                            "PDB file contained CONECT record to TER entry. "
                            "These are not included in bonds.")
                    else:
                        bonds.add(bond)

        bonds = tuple(bonds)

        return Bonds(bonds)


def _parse_conect(conect):
    """parse a CONECT record from pdbs

    Parameters
    ----------
    conect : str
        white space striped CONECT record

    Returns
    -------
    atom_id : int
        atom index of bond
    bonds : set
        atom ids of bonded atoms

    Raises
    ------
    RuntimeError
        Raised if ``conect`` is not a valid CONECT record
    """
    atom_id = int(conect[6:11])
    n_bond_atoms = len(conect[11:]) // 5

    try:
        if len(conect[11:]) % n_bond_atoms != 0:
            raise RuntimeError("Bond atoms aren't aligned proberly for CONECT "
                               "record: {}".format(conect))
    except ZeroDivisionError:
        # Conect record with only one entry (CONECT A\n)
        warnings.warn("Found CONECT record with single entry, ignoring this")
        return atom_id, []  # return empty list to allow iteration over nothing

    bond_atoms = (int(conect[11 + i * 5: 16 + i * 5]) for i in
                  range(n_bond_atoms))
    return atom_id, bond_atoms
