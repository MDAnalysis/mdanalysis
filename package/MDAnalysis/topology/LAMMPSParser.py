# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.mdanalysis.org
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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""
LAMMPSParser
============

Parses data files for LAMMPS_.

.. _LAMMPS: http://lammps.sandia.gov/

Classes
-------

.. autoclass:: DATAParser
   :members:
   :inherited-members:

Deprecated classes
------------------

.. autoclass:: LAMMPSDataConverter
   :members:

"""
from __future__ import absolute_import, print_function

from six.moves import range

import numpy as np
import logging
import string
import functools

from . import guessers
from ..lib.util import openany, conv_float
from ..lib.mdamath import triclinic_box
from .base import TopologyReaderBase, squash_by
from ..core.topology import Topology
from ..core.topologyattrs import (
    Atomtypes,
    Atomids,
    Angles,
    Bonds,
    Charges,
    Dihedrals,
    Impropers,
    Masses,
    Resids,
    Resnums,
    Segids,
)

logger = logging.getLogger("MDAnalysis.topology.LAMMPS")


# Sections will all start with one of these words
# and run until the next section title
SECTIONS = set([
    'Atoms',  # Molecular topology sections
    'Velocities',
    'Masses',
    'Ellipsoids',
    'Lines',
    'Triangles',
    'Bodies',
    'Bonds',  # Forcefield sections
    'Angles',
    'Dihedrals',
    'Impropers',
    'Pair',
    'Pair LJCoeffs',
    'Bond Coeffs',
    'Angle Coeffs',
    'Dihedral Coeffs',
    'Improper Coeffs',
    'BondBond Coeffs',  # Class 2 FF sections
    'BondAngle Coeffs',
    'MiddleBondTorsion Coeffs',
    'EndBondTorsion Coeffs',
    'AngleTorsion Coeffs',
    'AngleAngleTorsion Coeffs',
    'BondBond13 Coeffs',
    'AngleAngle Coeffs',
])
# We usually check by splitting around whitespace, so check
# if any SECTION keywords will trip up on this
# and add them
for val in list(SECTIONS):
    if len(val.split()) > 1:
        SECTIONS.add(val.split()[0])


HEADERS = set([
    'atoms',
    'bonds',
    'angles',
    'dihedrals',
    'impropers',
    'atom types',
    'bond types',
    'angle types',
    'dihedral types',
    'improper types',
    'extra bond per atom',
    'extra angle per atom',
    'extra dihedral per atom',
    'extra improper per atom',
    'extra special per atom',
    'ellipsoids',
    'lines',
    'triangles',
    'bodies',
    'xlo xhi',
    'ylo yhi',
    'zlo zhi',
    'xy xz yz',
])


class DATAParser(TopologyReaderBase):
    """Parse a LAMMPS DATA file for topology and coordinates.

    Note that LAMMPS_ DATA files can be used standalone.

    Both topology and coordinate parsing functionality is kept in this
    class as the topology and coordinate reader share many common
    functions

    The parser implements the `LAMMPS DATA file format`_ but only for
    the LAMMPS `atom_style`_ *full* (numeric ids 7, 10) and
    *molecular* (6, 9).

    .. _LAMMPS DATA file format: :http://lammps.sandia.gov/doc/2001/data_format.html
    .. _`atom_style`: http://lammps.sandia.gov/doc/atom_style.html

    .. versionadded:: 0.9.0
    """
    format = 'DATA'

    def iterdata(self):
        with openany(self.filename, 'r') as f:
            for line in f:
                line = line.partition('#')[0].strip()
                if line:
                    yield line

    def grab_datafile(self):
        """Split a data file into dict of header and sections

        Returns
        -------
        header - dict of header section: value
        sections - dict of section name: content
        """
        f = list(self.iterdata())

        starts = [i for i, line in enumerate(f)
                  if line.split()[0] in SECTIONS]
        starts += [None]

        header = {}
        for line in f[:starts[0]]:
            for token in HEADERS:
                if line.endswith(token):
                    header[token] = line.split(token)[0]
                    continue

        sects = {f[l]:f[l+1:starts[i+1]]
                 for i, l in enumerate(starts[:-1])}

        return header, sects

    def parse(self):
        """Parses a LAMMPS_ DATA file.

        Returns
        -------
        MDAnalysis Topology object.
        """
        # Can pass atom_style to help parsing
        atom_style = self.kwargs.get('atom_style', None)

        head, sects = self.grab_datafile()

        try:
            masses = self._parse_masses(sects['Masses'])
        except KeyError:
            masses = None

        try:
            top = self._parse_atoms(sects['Atoms'], masses)
        except KeyError:
            raise ValueError("Data file was missing Atoms section")

        # create mapping of id to index (ie atom id 10 might be the 0th atom)
        mapping = {atom_id: i for i, atom_id in enumerate(top.ids.values)}

        for attr, L, nentries in [
                (Bonds, 'Bonds', 2),
                (Angles, 'Angles', 3),
                (Dihedrals, 'Dihedrals', 4),
                (Impropers, 'Impropers', 4)
        ]:
            try:
                type, sect = self._parse_bond_section(sects[L], nentries, mapping)
            except KeyError:
                pass
            else:
                top.add_TopologyAttr(attr(sect, type))

        return top

    def read_DATA_timestep(self, n_atoms, TS_class, TS_kwargs):
        """Read a DATA file and try and extract x, v, box.

        - positions
        - velocities (optional)
        - box information

        Fills this into the Timestep object and returns it

        .. versionadded:: 0.9.0
        """
        header, sects = self.grab_datafile()

        unitcell = self._parse_box(header)

        try:
            positions, ordering = self._parse_pos(sects['Atoms'])
        except KeyError:
            raise IOError("Position information not found")

        if 'Velocities' in sects:
            velocities = self._parse_vel(sects['Velocities'], ordering)
        else:
            velocities = None

        ts = TS_class.from_coordinates(positions,
                                       velocities=velocities,
                                       **TS_kwargs)
        ts.dimensions = unitcell

        return ts

    def _parse_pos(self, datalines):
        """Strip coordinate info into np array"""
        pos = np.zeros((len(datalines), 3), dtype=np.float32)
        # TODO: could maybe store this from topology parsing?
        # Or try to reach into Universe?
        # but ugly because assumes lots of things, and Reader should be standalone
        ids = np.zeros(len(pos), dtype=np.int32)

        for i, line in enumerate(datalines):
            line = line.split()
            n = len(line)
            ids[i] = line[0]

            if n in (7, 10):
                pos[i] = line[4:7]
            elif n in (6, 9):
                pos[i] = line[3:6]

        order = np.argsort(ids)
        pos = pos[order]

        # return order for velocities
        return pos, order

    def _parse_vel(self, datalines, order):
        """Strip velocity info into np array

        Parameters
        ----------
        datalines : list
          list of strings from file
        order : np.array
          array which rearranges the velocities into correct order
          (from argsort on atom ids)

        Returns
        -------
        velocities : np.ndarray
        """
        vel = np.zeros((len(datalines), 3), dtype=np.float32)

        for i, line in enumerate(datalines):
            line = line.split()
            vel[i] = line[1:4]

        vel = vel[order]

        return vel

    def _parse_bond_section(self, datalines, nentries, mapping):
        """Read lines and strip information

        Arguments
        ---------
        datalines : list
          the raw lines from the data file
        nentries : int
          number of integers per line
        mapping : dict
          converts atom_ids to index within topology
        """
        section = []
        type = []
        for line in datalines:
            line = line.split()
            # map to 0 based int
            section.append(tuple([mapping[int(x)] for x in line[2:2 + nentries]]))
            type.append(line[1])
        return tuple(type), tuple(section)

    def _parse_atoms(self, datalines, massdict=None):
        """Creates a Topology object

        Adds the following attributes
         - resid
         - type
         - masses (optional)
         - charge (optional)

        Lammps atoms can have lots of different formats,
        and even custom formats

        http://lammps.sandia.gov/doc/atom_style.html

        Treated here are
        - atoms with 7 fields (with charge) "full"
        - atoms with 6 fields (no charge) "molecular"

        Arguments
        ---------
        datalines - the relevent lines from the data file
        massdict - dictionary relating type to mass

        Returns
        -------
        top - Topology object
        """
        logger.info("Doing Atoms section")

        n_atoms = len(datalines)

        # Fields per line
        n = len(datalines[0].split())
        has_charge = True if n in [7, 10] else False

        # atom ids aren't necessarily sequential
        atom_ids = np.zeros(n_atoms, dtype=np.int32)
        types = np.zeros(n_atoms, dtype=object)
        resids = np.zeros(n_atoms, dtype=np.int32)
        if has_charge:
            charges = np.zeros(n_atoms, dtype=np.float32)

        for i, line in enumerate(datalines):
            line = line.split()

            # these numpy array are already typed correctly,
            # so just pass the raw strings
            # and let numpy handle the conversion
            atom_ids[i] = line[0]
            resids[i] = line[1]
            types[i] = line[2]
            if has_charge:
                charges[i] = line[3]

        # at this point, we've read the atoms section,
        # but it's still (potentially) unordered
        # TODO: Maybe we can optimise by checking if we need to sort
        # ie `if np.any(np.diff(atom_ids) > 1)`  but we want to search
        # in a generatorish way, np.any() would check everything at once
        order = np.argsort(atom_ids)
        atom_ids = atom_ids[order]
        types = types[order]
        resids = resids[order]
        if has_charge:
            charges = charges[order]

        attrs = []
        attrs.append(Atomtypes(types))
        if has_charge:
            attrs.append(Charges(charges))
        if massdict is not None:
            masses = np.zeros(n_atoms, dtype=np.float64)
            for i, at in enumerate(types):
                masses[i] = massdict[at]
            attrs.append(Masses(masses))
        else:
            # Guess them
            masses = guessers.guess_masses(types)
            attrs.append(Masses(masses, guessed=True))

        residx, resids = squash_by(resids)[:2]
        n_residues = len(resids)

        attrs.append(Atomids(atom_ids))
        attrs.append(Resids(resids))
        attrs.append(Resnums(resids.copy()))
        attrs.append(Segids(np.array(['SYSTEM'], dtype=object)))

        top = Topology(n_atoms, n_residues, 1,
                       attrs=attrs,
                       atom_resindex=residx)

        return top

    def _parse_masses(self, datalines):
        """Lammps defines mass on a per atom type basis.

        This reads mass for each type and stores in dict
        """
        logger.info("Doing Masses section")

        masses = {}
        for line in datalines:
            line = line.split()
            masses[line[0]] = float(line[1])

        return masses

    def _parse_box(self, header):
        x1, x2 = np.float32(header['xlo xhi'].split())
        x = x2 - x1
        y1, y2 = np.float32(header['ylo yhi'].split())
        y = y2 - y1
        z1, z2 = np.float32(header['zlo zhi'].split())
        z = z2 - z1

        if 'xy xz yz' in header:
            # Triclinic
            unitcell = np.zeros((3, 3), dtype=np.float32)

            xy, xz, yz = np.float32(header['xy xz yz'].split())

            unitcell[0][0] = x
            unitcell[1][0] = xy
            unitcell[1][1] = y
            unitcell[2][0] = xz
            unitcell[2][1] = yz
            unitcell[2][2] = z

            unitcell = triclinic_box(*unitcell)
        else:
            # Orthogonal
            unitcell = np.zeros(6, dtype=np.float32)
            unitcell[:3] = x, y, z
            unitcell[3:] = 90., 90., 90.

        return unitcell


@functools.total_ordering
class LAMMPSAtom(object):  # pragma: no cover
    __slots__ = ("index", "name", "type", "chainid", "charge", "mass", "_positions")

    def __init__(self, index, name, type, chain_id, charge=0, mass=1):
        self.index = index
        self.name = repr(type)
        self.type = type
        self.chainid = chain_id
        self.charge = charge
        self.mass = mass

    def __repr__(self):
        return "<LAMMPSAtom " + repr(self.index + 1) + ": name " + repr(self.type) + " of chain " + repr(
            self.chainid) + ">"

    def __lt__(self, other):
        return self.index < other.index

    def __eq__(self, other):
        return self.index == other.index

    def __hash__(self):
        return hash(self.index)

    def __getattr__(self, attr):
        if attr == 'pos':
            return self._positions[self.index]
        else:
            super(LAMMPSAtom, self).__getattribute__(attr)

    def __iter__(self):
        pos = self.pos
        return iter((self.index + 1, self.chainid, self.type, self.charge, self.mass, pos[0], pos[1], pos[2]))


class LAMMPSDataConverter(object):  # pragma: no cover
    """Class to parse a LAMMPS_ DATA file and convert it to PSF/PDB.

    The DATA file contains both topology and coordinate information.

    The :class:`LAMMPSDataConverter` class can extract topology information and
    coordinates from a LAMMPS_ data file. For instance, in order to
    produce a PSF file of the topology and a PDB file of the coordinates
    from a data file "lammps.data" you can use::

      from MDAnalysis.topology.LAMMPSParser import LAMPPSData
      d = LAMMPSDataConverter("lammps.data")
      d.writePSF("lammps.psf")
      d.writePDB("lammps.pdb")

    You can then read a trajectory (e.g. a LAMMPS DCD, see
    :class:`MDAnalysis.coordinates.LAMMPS.DCDReader`) with ::

      u = MDAnalysis.Unverse("lammps.psf", "lammps.dcd", format="LAMMPS")

    .. deprecated:: 0.9.0

    .. versionchanged:: 0.9.0
       Renamed from ``LAMMPSData`` to ``LAMMPSDataConverter``.
    """
    header_keywords = [
        "atoms", "bonds", "angles", "dihedrals", "impropers",
        "atom types", "bond types", "angle types",
        "dihedral types", "improper types",
        "xlo xhi", "ylo yhi", "zlo zhi"]

    connections = dict([
        ["Bonds", ("bonds", 3)],
        ["Angles", ("angles", 3)],
        ["Dihedrals", ("dihedrals", 4)],
        ["Impropers", ("impropers", 2)]])

    coeff = dict([
        ["Masses", ("atom types", 1)],
        ["Velocities", ("atoms", 3)],
        ["Pair Coeffs", ("atom types", 4)],
        ["Bond Coeffs", ("bond types", 2)],
        ["Angle Coeffs", ("angle types", 4)],
        ["Dihedral Coeffs", ("dihedral types", 3)],
        ["Improper Coeffs", ("improper types", 2)]])

    def __init__(self, filename=None):
        self.names = {}
        self.headers = {}
        self.sections = {}
        if filename is None:
            self.title = "LAMMPS data file"
        else:
            # Open and check validity
            with openany(filename, 'r') as file:
                file_iter = file.xreadlines()
                self.title = file_iter.next()
                # Parse headers
                headers = self.headers
                for l in file_iter:
                    line = l.strip()
                    if len(line) == 0:
                        continue
                    found = False
                    for keyword in self.header_keywords:
                        if line.find(keyword) >= 0:
                            found = True
                            values = line.split()
                            if keyword in ("xlo xhi", "ylo yhi", "zlo zhi"):
                                headers[keyword] = (float(values[0]), float(values[1]))
                            else:
                                headers[keyword] = int(values[0])
                    if found is False:
                        break

            # Parse sections
            # XXX This is a crappy way to do it
            with openany(filename, 'r') as file:
                file_iter = file.xreadlines()
                # Create coordinate array
                positions = np.zeros((headers['atoms'], 3), np.float64)
                sections = self.sections
                for l in file_iter:
                    line = l.strip()
                    if len(line) == 0:
                        continue
                    if line in self.coeff:
                        h, numcoeff = self.coeff[line]
                        # skip line
                        file_iter.next()
                        data = []
                        for i in range(headers[h]):
                            fields = file_iter.next().strip().split()
                            data.append(tuple([conv_float(el) for el in fields[1:]]))
                        sections[line] = data
                    elif line in self.connections:
                        h, numfields = self.connections[line]
                        # skip line
                        file_iter.next()
                        data = []
                        for i in range(headers[h]):
                            fields = file_iter.next().strip().split()
                            data.append(tuple(np.int64(fields[1:])))
                        sections[line] = data
                    elif line == "Atoms":
                        file_iter.next()
                        data = []
                        for i in range(headers["atoms"]):
                            fields = file_iter.next().strip().split()
                            index = int(fields[0]) - 1
                            a = LAMMPSAtom(index=index, name=fields[2], type=int(fields[2]), chain_id=int(fields[1]),
                                           charge=float(fields[3]))
                            a._positions = positions
                            data.append(a)
                            positions[index] = np.array([float(fields[4]), float(fields[5]), float(fields[6])])
                        sections[line] = data
                    elif line == "Masses":
                        file_iter.next()
                        data = []
                        for i in range(headers["atom type"]):
                            fields = file_iter.next().strip().split()
                            print("help")
                self.positions = positions

    def writePSF(self, filename, names=None):
        """Export topology information to a simple PSF file."""
        # Naveen formatted -- works with MDAnalysis verison 52
        #psf_atom_format = "   %5d %-4s %-4d %-4s %-4s %-4s %10.6f      %7.4f            %1d\n"
        # Liz formatted -- works with MDAnalysis verison 59
        #psf_atom_format = "%8d %4.4s %-4.4s %-4.4s %-4.4s %-4.4s %16.8e %1s %-7.4f %7.7s %s\n"
        # Oli formatted -- works with MDAnalysis verison 81
        psf_atom_format = "%8d %4s %-4s %4s %-4s% 4s %-14.6f%-14.6f%8s\n"
        with openany(filename, 'w') as file:
            file.write("PSF\n\n")
            file.write(string.rjust('0', 8) + ' !NTITLE\n\n')
            file.write(string.rjust(str(len(self.sections["Atoms"])), 8) + ' !NATOM\n')
            #print self.sections["Masses"]
            for i, atom in enumerate(self.sections["Atoms"]):
                if names is not None:
                    resname, atomname = names[i]
                else:
                    resname, atomname = 'TEMP', 'XXXX'
                for j, liz in enumerate(self.sections["Masses"]):
                    liz = liz[0]
                    #print j+1, atom.type, liz
                    if j + 1 == atom.type:
                        line = [
                            i + 1, 'TEMP',
                            str(atom.chainid), resname, atomname, str(atom.type + 1), atom.charge,
                            float(liz), 0.]
                    else:
                        continue
                #print line
                file.write(psf_atom_format % tuple(line))

            file.write("\n")
            num_bonds = len(self.sections["Bonds"])
            bond_list = self.sections["Bonds"]
            file.write(string.rjust(str(num_bonds), 8) + ' !NBOND\n')
            for index in range(0, num_bonds, 4):
                try:
                    bonds = bond_list[index:index + 4]
                except IndexError:
                    bonds = bond_list[index:-1]
                bond_line = [string.rjust(str(bond[1]), 8) + string.rjust(str(bond[2]), 8) for bond in bonds]
                file.write(''.join(bond_line) + '\n')

    def writePDB(self, filename):
        """Export coordinates to a simple PDB file."""
        atom_format = "%6s%.5s %4s %4s %.4s    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n"
        p = self.positions
        with openany(filename, 'w') as file:
            for i, atom in enumerate(self.sections["Atoms"]):
                line = [
                    "ATOM  ", str(i + 1), 'XXXX', 'TEMP', str(atom.type + 1), p[i, 0], p[i, 1], p[i, 2], 0.0, 0.0,
                    str(atom.type)]
                file.write(atom_format % tuple(line))
