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
LAMMPSParser
============

Parses data_ or dump_ files produced by LAMMPS_.

.. _LAMMPS: http://lammps.sandia.gov/
.. _data: DATA file format: :https://docs.lammps.org/read_data.html
.. _dump: http://lammps.sandia.gov/doc/dump.html

.. versionchanged:: 1.0.0
   Deprecated :class:`LAMMPSDataConverter` has now been removed.


.. _atom_style_kwarg:

Atom styles
-----------

The parsers and readers for LAMMPS DATA files can handle any `atom_style`_
defined in LAMMPS; however, they only parse `id` ('atom-ID' in LAMMPS), `type`
('atom-type' in LAMMPS), `resid` ('molecule-ID' in LAMMPS), `charge` ('q' or
'charge' in LAMMPS), `x`, `y` and `z` attributes. The resid and charge
attributes are optional and any other specified attribute define an atom style
will be ignored.

Valid atom styles as defined in `data`_::

    'angle', 'atomic', 'body', 'bond', 'charge', 'dipole', 'dpd', 'edpd',
    'electron', 'ellipsoid', 'full', 'line', 'mdpd', 'mesont', 'molecular',
    'peri', 'smd', 'sph', 'sphere', 'spin', 'tdpd', 'template', 'tri',
    'wavepacket', 'hybrid'


Examples
--------
To parse a DATA file in 'atomic'::

  Atoms # atomic

  1 1 3.7151744275286681e+01 1.8684434743140471e+01 1.9285127961842125e+01 0 0 0


The following code could be used::

  >>> import MDAnalysis as mda
  >>>
  >>> u = mda.Universe('myfile.data', atom_style='atomic')


.. _`atom_style`: http://lammps.sandia.gov/doc/atom_style.html

Classes
-------

.. autoclass:: DATAParser
   :members:
   :inherited-members:

.. autoclass:: LammpsDumpParser
   :members:


"""
import numpy as np
import logging
import functools
import warnings

from . import guessers
from ..lib.util import openany
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

# NOTE: The SECTIONS, HEADERS, and ATOM_STYLES are only used in DATAParser, so
# it is better to defined them as class attributes in the DATAParser, not.

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
    'Pair Coeffs',
    'PairLJ Coeffs',
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

# List of valid atom styles for Atoms section:
# Each atom style can optionally have three image flags: nx, ny, nz
# This mapping is used between LAMMPS and MDAnalysis attributes:
# LAMMPS -> MDAnalysis: atom-ID -> id, molecule-ID -> segid, atom-type -> type
# q -> charge.
# Some commments about the some atom styles:
# 'smd': molecule is defined differently from molecule-ID.
# 'tdpd': this style can have n+5 columns for n cc species.
# 'template': the order of id and resid is swapped.
# 'hybrid': this style can have n+5 columns for n sub-styles.
# See below for an up-to-date list of valid atom styles and their usage:
# https://docs.lammps.org/atom_style.html
# See below for an up-to-date list of valid atom styles and how they defined in
# a LAMMPS DATA file:
# https://docs.lammps.org/read_data.html

# NOTE: In MDAnalysis, 'id', 'type', 'x', 'y', 'z' are the mandatory
# attributes, and 'segid' and 'charge' are the optional ones.
# What about the other attributes defined in LAMMPS and cannot/are not
# mapped to corresponding attributes in MDAnalysis? Are they ignored as
# mentioned in the doc above?
# An idea is to keep all the atom styles defined below, and just read in
# the attributes that are compatible with MDAnalysis and ignore the rest.
# This way, MDAnalysis can read all types of LAMMPS data file but only read
# in the 5 mendatory and 2 optional attributes, and warns about how it handles
# atom styles that have other attributes than these 7 attributes.

ATOM_STYLES = {
    'angle': ['id', 'segid', 'type', 'x', 'y', 'z'],
    'atomic': ['id', 'type', 'x', 'y', 'z'],
    'body': ['id', 'type', 'bodyflag', 'mass', 'x', 'y', 'z'],
    'bond': ['id', 'resid', 'type', 'x', 'y', 'z'],
    'charge': ['id', 'type', 'charge', 'x', 'y', 'z'],
    'dipole': ['id', 'type', 'charge', 'x', 'y', 'z', 'mux', 'muy', 'muz'],
    'dpd': ['id', 'type', 'theta', 'x', 'y', 'z'],
    'edpd': ['id', 'type', 'edpd_temp', 'edpd_cv', 'x', 'y', 'z'],
    'electron': ['id', 'type', 'charge', 'spin', 'eradius', 'x', 'y', 'z'],
    'ellipsoid': ['id', 'type', 'ellipsoidflag', 'density', 'x', 'y', 'z'],
    'full': ['id', 'resid', 'type', 'charge', 'x', 'y', 'z'],
    'line': ['id', 'resid', 'type', 'lineflag', 'density', 'x', 'y', 'z'],
    'mdpd': ['id', 'type', 'rho', 'x', 'y', 'z'],
    'mesont': ['id', 'segid', 'type', 'bond_nt', 'mass', 'mradius', 'mlength',
               'buckling', 'x', 'y', 'z'],
    'molecular': ['id', 'segid', 'type', 'x', 'y', 'z'],
    'peri': ['id', 'type', 'volume', 'density', 'x', 'y', 'z'],
    'smd': ['id', 'type', 'molecule', 'volume', 'mass', 'kradius', 'cradius',
            'x0', 'y0', 'z0', 'x', 'y', 'z'],
    'sph': ['id', 'type', 'rho', 'esph', 'cv', 'x', 'y', 'z'],
    'sphere': ['id', 'type', 'diameter', 'density', 'x', 'y', 'z'],
    'spin': ['id', 'type', 'x', 'y', 'z', 'spx', 'spy', 'spz', 'sp'],
    'tdpd': ['id', 'type', 'x', 'y', 'z'],
    'template': ['id', 'type', 'resid', 'template-index', 'template-atom', 'x',
                 'y', 'z'],
    'tri': ['id', 'segid', 'type', 'triangleflag', 'density', 'x', 'y', 'z'],
    'wavepacket': ['id', 'type', 'charge', 'spin', 'eradius', 'etag', 'cs_re',
                   'cs_im', 'x', 'y', 'z'],
    'hybrid': ['id', 'type', 'x', 'y', 'z']
}


class DATAParser(TopologyReaderBase):
    """Parse a LAMMPS DATA file for topology and coordinates.

    Note that LAMMPS DATA files can be used standalone.

    Both topology and coordinate parsing functionality is kept in this
    class as the topology and coordinate reader share many common
    functions

    By default the parser expects either *atomic* or *full* `atom_style`
    however this can be by passing an `atom_style` keyword argument,
    see :ref:`atom_style_kwarg`.

    .. versionadded:: 0.9.0
    """
    format = 'DATA'

    def iterdata(self):
        with openany(self.filename) as f:
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

        sects = {f[l]: f[l+1:starts[i+1]]
                 for i, l in enumerate(starts[:-1])}

        return header, sects

    @staticmethod
    def _interpret_atom_style(atom_style):
        """Map field names in a given `atom_style` to their corresponding
        locations.

        Required attributes: id, type, x, y, z
        Optional attributes: resid, charge
        Note: other attributes are ignored.

        Exammples
        ---------
        "full" atom_style has the following attributes:
            "id resid type charge x y z"
        `_interpret_atom_style` returns this dict for "full" atom_style:
            {'id': 0,
             'resid': 1,
             'type': 2,
             'charge': 3,
             'x': 4,
             'y': 5,
             'z': 6,
            }

        Parameters
        ----------
        atom_style: str, default 'full'
            One of the following atom styles:
            {'angle', 'atomic', 'body', 'bond', 'charge', 'dipole', 'dpd',
             'edpd', 'electron', 'ellipsoid', 'full', 'line', 'mdpd', 'mesont',
             'molecular', 'peri', 'smd', 'sph', 'sphere', 'spin', 'tdpd',
             'template', 'tri', 'wavepacket', 'hybrid'}

        Returns
        -------
        dict:
            a map between the attributes' names and attributes' locations.
        """
        style_dict = {}

        for loc, attr in enumerate(ATOM_STYLES[atom_style]):
            style_dict[attr] = loc

        return style_dict

    def parse(self, **kwargs):
        """Parses a LAMMPS_ DATA file.

        Returns
        -------
        MDAnalysis Topology object.
        """
        if kwargs['atom_style'] in ATOM_STYLES.keys():
            self.style_dict = self._interpret_atom_style(kwargs['atom_style'])
        else:
            atom_styles_string = "'" + "', '".join(
                ATOM_STYLES.keys()) + "'"
            raise ValueError(
                f"'{kwargs['atom_style']}' "
                "is not a valid 'atom_style'. Please select one of "
                f"{atom_styles_string} atom_styles.")

        head, sects = self.grab_datafile()

        try:
            masses = self._parse_masses(sects['Masses'])
        except KeyError:
            masses = None

        if 'Atoms' not in sects:
            raise ValueError("Data file was missing Atoms section")

        try:
            top = self._parse_atoms(sects['Atoms'], masses)
        except Exception:
            errmsg = (
                "Failed to parse atoms section.  You can supply a description "
                "of the atom_style as a keyword argument, "
                "eg mda.Universe(..., atom_style='id resid x y z')")
            raise ValueError(errmsg) from None

        # create mapping of id to index (ie atom id 10 might be the 0th atom)
        mapping = {atom_id: i for i, atom_id in enumerate(top.ids.values)}

        for attr, L, nentries in [
                (Bonds, 'Bonds', 2),
                (Angles, 'Angles', 3),
                (Dihedrals, 'Dihedrals', 4),
                (Impropers, 'Impropers', 4)
        ]:
            try:
                type, sect = self._parse_bond_section(
                    sects[L], nentries, mapping)
            except KeyError:
                type, sect = [], []

            top.add_TopologyAttr(attr(sect, type))

        return top

    def read_DATA_timestep(self, n_atoms, TS_class, TS_kwargs,
                           atom_style=None):
        """Read a DATA file and try and extract x, v, box.

        - positions
        - velocities (optional)
        - box information

        Fills this into the Timestep object and returns it

        .. versionadded:: 0.9.0
        .. versionchanged:: 0.18.0
           Added atom_style kwarg
        """
        if atom_style is None:
            self.style_dict = None
        else:
            self.style_dict = self._interpret_atom_style(atom_style)

        header, sects = self.grab_datafile()

        unitcell = self._parse_box(header)

        try:
            positions, ordering = self._parse_pos(sects['Atoms'])
        except KeyError as err:
            errmsg = f"Position information not found: {err}"
            raise IOError(errmsg) from None

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
        # but ugly since assumes lots of things, and Reader should be standalone
        ids = np.zeros(len(pos), dtype=np.int32)

        if self.style_dict is None:
            if len(datalines[0].split()) in (7, 10):
                style_dict = {'id': 0, 'x': 4, 'y': 5, 'z': 6}
            else:
                style_dict = {'id': 0, 'x': 3, 'y': 4, 'z': 5}
        else:
            style_dict = self.style_dict

        for i, line in enumerate(datalines):
            line = line.split()

            ids[i] = line[style_dict['id']]

            pos[i, :] = [line[style_dict['x']],
                         line[style_dict['y']],
                         line[style_dict['z']]]

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

        Returns
        -------
        types : tuple of strings
          type of the bond/angle/dihedral/improper
        indices : tuple of ints
          indices of atoms involved
        """
        section = []
        type = []
        for line in datalines:
            line = line.split()
            # map to 0 based int
            section.append(
                tuple([mapping[int(x)] for x in line[2:2 + nentries]])
                )
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

        if self.style_dict is None:
            sd = {'id': 0,
                  'resid': 1,
                  'type': 2}
            # Fields per line
            n = len(datalines[0].split())
            if n in (7, 10):
                sd['charge'] = 3
        else:
            sd = self.style_dict

        has_charge = 'charge' in sd
        has_resid = 'resid' in sd

        # atom ids aren't necessarily sequential
        atom_ids = np.zeros(n_atoms, dtype=np.int32)
        types = np.zeros(n_atoms, dtype=object)
        if has_resid:
            resids = np.zeros(n_atoms, dtype=np.int32)
        else:
            resids = np.ones(n_atoms, dtype=np.int32)
        if has_charge:
            charges = np.zeros(n_atoms, dtype=np.float32)

        for i, line in enumerate(datalines):
            line = line.split()

            # these numpy array are already typed correctly,
            # so just pass the raw strings
            # and let numpy handle the conversion
            atom_ids[i] = line[sd['id']]
            if has_resid:
                resids[i] = line[sd['resid']]
            types[i] = line[sd['type']]
            if has_charge:
                charges[i] = line[sd['charge']]

        # at this point, we've read the atoms section,
        # but it's still (potentially) unordered
        # TODO: Maybe we can optimise by checking if we need to sort
        # ie `if np.any(np.diff(atom_ids) > 1)`  but we want to search
        # in a generatorish way, np.any() would check everything at once
        order = np.argsort(atom_ids)
        atom_ids = atom_ids[order]
        types = types[order]
        if has_resid:
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


class LammpsDumpParser(TopologyReaderBase):
    """Parses Lammps ascii dump files in 'atom' format

    Only reads atom ids.  Sets all masses to 1.0.

    .. versionadded:: 0.19.0
    """
    format = 'LAMMPSDUMP'

    def parse(self, **kwargs):
        with openany(self.filename) as fin:
            fin.readline()  # ITEM TIMESTEP
            fin.readline()  # 0

            fin.readline()  # ITEM NUMBER OF ATOMS
            natoms = int(fin.readline())

            fin.readline()  # ITEM BOX
            fin.readline()  # x
            fin.readline()  # y
            fin.readline()  # z

            indices = np.zeros(natoms, dtype=int)
            types = np.zeros(natoms, dtype=object)

            fin.readline()  # ITEM ATOMS
            for i in range(natoms):
                idx, atype, _, _, _ = fin.readline().split()

                indices[i] = idx
                types[i] = atype

        order = np.argsort(indices)
        indices = indices[order]
        types = types[order]

        attrs = []
        attrs.append(Atomids(indices))
        attrs.append(Atomtypes(types))
        attrs.append(Masses(np.ones(natoms, dtype=np.float64), guessed=True))
        warnings.warn('Guessed all Masses to 1.0')
        attrs.append(Resids(np.array([1], dtype=int)))
        attrs.append(Resnums(np.array([1], dtype=int)))
        attrs.append(Segids(np.array(['SYSTEM'], dtype=object)))

        return Topology(natoms, 1, 1, attrs=attrs)


@functools.total_ordering
class LAMMPSAtom(object):  # pragma: no cover
    __slots__ = ("index", "name", "type", "chainid", "charge", "mass",
                 "_positions")

    def __init__(self, index, name, type, chain_id, charge=0, mass=1):
        self.index = index
        self.name = repr(type)
        self.type = type
        self.chainid = chain_id
        self.charge = charge
        self.mass = mass

    def __repr__(self):
        return "<LAMMPSAtom " + repr(self.index + 1) + ": name " + repr(self.type) + " of chain " + repr(self.chainid) + ">"

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
        return iter((self.index + 1, self.chainid, self.type, self.charge,
                     self.mass, pos[0], pos[1], pos[2]))
