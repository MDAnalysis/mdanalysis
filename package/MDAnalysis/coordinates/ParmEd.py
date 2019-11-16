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

"""ParmEd structure reader --- :mod:`MDAnalysis.coordinates.ParmEd`
================================================================

Reads coordinates data from a ParmEd Structure.

Classes
-------

.. autoclass:: ParmEdReader
   :members:

"""
from __future__ import absolute_import

import parmed as pmd

from . import base
from ..core import flags
from ..core.universe import Universe



class ParmEdReader(base.SingleFrameReaderBase):
    """Coordinate reader for ParmEd."""
    format = 'PARMED'

    # Structure.coordinates always in Angstrom
    units = {'time': None, 'length': 'Angstrom'}

    def _read_first_frame(self):
        self.n_atoms = len(self.filename.atoms)

        self.ts = ts = self._Timestep(self.n_atoms,
                                      **self._ts_kwargs)
        
        if self.filename.coordinates is not None:
            ts._pos = self.filename.coordinates

        if self.filename.box is not None:
            # optional field
            ts.dimensions = self.filename.box
        else:
            ts._unitcell = None

        ts.frame = 0

        return ts


MDA2PMD = {
    'tempfactor': 'bfactor',
    'gbscreen': 'screen',
    'altLoc': 'altloc',
    'nbindex': 'nb_idx',
    'solventradius': 'solvent_radius',
    'id': 'number'
}

SYMB2Z = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99}


class ParmEdWriter(base.WriterBase):
    """Convert MDAnalysis AtomGroup or Universe to ParmEd Structure
    """

    format = 'PARMED'
    units = {'time': None, 'length': 'Angstrom'}

    def __init__(self, n_atoms=None, convert_units=None, **kwargs):
        """Set up a ParmEdWriter with a precision of 3 decimal places.

        Parameters
        -----------
        filename : str
            output filename

        n_atoms : int (optional)
            number of atoms

        """
        self.n_atoms = n_atoms

        if convert_units is None:
            convert_units = flags['convert_lengths']
        self.convert_units = convert_units  # convert length and time to base units

    def write(self, obj):
        """Write selection at current trajectory frame to ParmEdStructure.

        Parameters
        -----------
        obj : AtomGroup or Universe or :class:`Timestep`
        """
        # write() method that complies with the Trajectory API

        try:

            # make sure to use atoms (Issue 46)
            ag_or_ts = obj.atoms
            # can write from selection == Universe (Issue 49)

        except AttributeError:
            if isinstance(obj, base.Timestep):
                ag_or_ts = obj.copy()
            else:
                raise_from(TypeError("No Timestep found in obj argument"), None)

        # Check for topology information
        missing_topology = []
        try:
            names = ag_or_ts.names
        except (AttributeError, NoDataError):
            names = itertools.cycle(('X',))
            missing_topology.append('names')
        try:
            resnames = ag_or_ts.resnames
        except (AttributeError, NoDataError):
            resnames = itertools.cycle(('UNK',))
            missing_topology.append('resnames')
        
        if missing_topology:
            warnings.warn(
                "Supplied AtomGroup was missing the following attributes: "
                "{miss}. These will be written with default values. "
                "Alternatively these can be supplied as keyword arguments."
                "".format(miss=', '.join(missing_topology)))
        
        try:
            positions = ag_or_ts.positions
        except:
            positions = [None]*ag_or_ts.n_atoms
        
        try:
            velocities = ag_or_ts.velocities
        except:
            velocities = [None]*ag_or_ts.n_atoms

        atom_kwargs = []
        for atom, name, resname, xyz, vel in zip(ag_or_ts, names, resnames, positions, velocities):
            akwargs = {'name': name}
            chain_seg = {'segid': atom.segid}
            for attrname in ('mass', 'charge', 'type',
                             'altLoc', 'tempfactor',
                             'occupancy', 'gbscreen', 'solventradius',
                             'nbindex', 'rmin', 'epsilon', 'rmin14',
                             'epsilon14', 'id'):
                try:
                    akwargs[MDA2PMD.get(attrname, attrname)] = getattr(atom, attrname)
                except AttributeError:
                    pass
                try:
                    el = atom.element.lower().capitalize()
                    akwargs['atomic_number'] = SYMB2Z[el]
                except (KeyError, AttributeError):
                    try:
                        akwargs['atomic_number'] = SYMB2Z[atom.type.lower().capitalize()]
                    except (KeyError, AttributeError):
                        pass
                try:
                    chain_seg['chain'] = atom.chainID
                except AttributeError:
                    pass
                try:
                    chain_seg['inscode'] = atom.icode
                except AttributeError:
                    pass
            atom_kwargs.append((akwargs, resname, atom.resid, chain_seg, xyz, vel))
        

        struct = pmd.Structure()
        
        for akwarg, resname, resid, kw, xyz, vel in atom_kwargs:
            atom = pmd.Atom(**akwarg)
            if xyz is not None:
                atom.xx, atom.xy, atom.xz = xyz

            if vel is not None:
                atom.vx, atom.vy, atom.vz = vel
            struct.add_atom(atom, resname, resid, **kw)
        
        try:
            struct.box = ag_or_ts.dimensions
        except AttributeError:
            struct.box = None
        
        
        
        try:
            params = ag_or_ts.bonds
        except AttributeError:
            pass
        else:
            for p in params:
                atoms = [struct.atoms[i] for i in p.indices]
                try:
                    for obj in p.type:
                        bond = pmd.Bond(*atoms, type=obj.type, order=obj.order)
                        struct.bonds.append(bond)
                except TypeError:
                    order = p.order if p.order is not None else 1
                    bond = pmd.Bond(*atoms, order=order)
                    struct.bonds.append(bond)
        
        for param, pmdtype, trackedlist in (
            ('ureybradleys', pmd.UreyBradley, struct.urey_bradleys),
            ('angles', pmd.Angle, struct.angles),
            ('dihedrals', pmd.Dihedral, struct.dihedrals),
            ('impropers', pmd.Improper, struct.impropers),
            ('cmaps', pmd.Cmap, struct.cmaps)
        ):
            try:
                values = getattr(ag_or_ts, param)
            except AttributeError:
                return
            
            for v in values:
                atoms = [struct.atoms[i] for i in v.indices]

                try:
                    for parmed_obj in v.type:
                        p = pmdtype(*atoms, type=parmed_obj.type)
                        trackedlist.append(p)
                except TypeError:
                    p = pmdtype(*atoms)
                    trackedlist.append(p)
        return struct