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

"""ParmEd structure I/O --- :mod:`MDAnalysis.converters.ParmEd`
================================================================

Read coordinates data from a `ParmEd <https://parmed.github.io/ParmEd/html>`_ :class:`parmed.structure.Structure`
with :class:`ParmEdReader` into a MDAnalysis Universe. Convert it back to a
:class:`parmed.structure.Structure` with :class:`ParmEdConverter`.

Example
-------

ParmEd has some neat functions. One such is `HMassRepartition`_.
This function changes the mass of the hydrogens in your system to your desired
value. It then adjusts the mass of the atom to which it is bonded by the same
amount, so that the total mass is unchanged. ::

    >>> import MDAnalysis as mda
    >>> from MDAnalysis.tests.datafiles import PRM
    >>> u = mda.Universe(PRM)
    >>> u.atoms.masses[:10]
    array([14.01 ,  1.008,  1.008,  1.008, 12.01 ,  1.008, 12.01 ,  1.008,
        1.008,  1.008])

We can convert our universe to a ParmEd structure to change our hydrogen
masses. ::

    >>> import parmed.tools as pmt
    >>> parm = u.atoms.convert_to('PARMED')
    >>> hmass = pmt.HMassRepartition(parm, 5)  # convert to 5 daltons
    >>> hmass.execute()

We can then convert it back to an MDAnalysis Universe for further analysis. ::

    >>> u2 = mda.Universe(parm)
    >>> u2.atoms.masses[:10]
    array([2.03399992, 5.        , 5.        , 5.        , 8.01799965,
       5.        , 0.034     , 5.        , 5.        , 5.        ])

.. _`HMassRepartition`: http://parmed.github.io/ParmEd/html/parmed.html#hmassrepartition



Classes
-------

.. autoclass:: ParmEdReader
   :members:

.. autoclass:: ParmEdConverter
   :members:


.. versionchanged:: 2.0.0
   The ParmEdReader and ParmEdConverter classes were moved from :mod:`~MDAnalysis.coordinates`
   to :mod:`~MDAnalysis.converters`

"""
import functools
import itertools
import warnings

from ..guesser.tables import SYMB2Z
import numpy as np
from numpy.lib import NumpyVersion

from . import base
from ..coordinates.base import SingleFrameReaderBase
from ..core.universe import Universe
from ..exceptions import NoDataError


class ParmEdReader(SingleFrameReaderBase):
    """Coordinate reader for ParmEd."""
    format = 'PARMED'

    # Structure.coordinates always in Angstrom
    units = {'time': None, 'length': 'Angstrom'}

    @staticmethod
    def _format_hint(thing):
        """Can this reader read *thing*?

        .. versionadded:: 1.0.0
        """
        try:
            import parmed as pmd
        except ImportError:
            # if we can't import parmed, it's probably not parmed
            return False
        else:
            return isinstance(thing, pmd.Structure)

    def _read_first_frame(self):
        self.n_atoms = len(self.filename.atoms)

        self.ts = ts = self._Timestep(self.n_atoms,
                                      **self._ts_kwargs)

        if self.filename.coordinates is not None:
            ts._pos = self.filename.coordinates

        # optional field
        ts.dimensions = self.filename.box

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


def get_indices_from_subset(i, atomgroup=None, universe=None):
    return atomgroup[universe.atoms[i]]


class ParmEdConverter(base.ConverterBase):
    """Convert MDAnalysis AtomGroup or Universe to ParmEd :class:`~parmed.structure.Structure`.

    Example
    -------

    .. code-block:: python

        import parmed as pmd
        import MDAnalysis as mda
        from MDAnalysis.tests.datafiles import GRO
        pgro = pmd.load_file(GRO)
        mgro = mda.Universe(pgro)
        parmed_subset = mgro.select_atoms('resname SOL').convert_to('PARMED')


    """

    lib = 'PARMED'
    units = {'time': None, 'length': 'Angstrom'}

    def convert(self, obj):
        """Write selection at current trajectory frame to :class:`~parmed.structure.Structure`.

        Parameters
        -----------
        obj : AtomGroup or Universe or :class:`Timestep`
        """
        try:
            # TODO: remove this guard when parmed has a release
            # that supports NumPy 2
            if NumpyVersion(np.__version__) < "2.0.0":
                import parmed as pmd
            else:
                raise ImportError
        except ImportError:
            if NumpyVersion(np.__version__) >= "2.0.0":
                ermsg = "ParmEd is not compatible with NumPy 2.0+"
            else:
                ermsg = ("ParmEd is required for ParmEdConverter but is not "
                         "installed. Try installing it with \n"
                         "pip install parmed")
            raise ImportError(errmsg)
        try:
            # make sure to use atoms (Issue 46)
            ag_or_ts = obj.atoms
        except AttributeError:
            raise TypeError("No atoms found in obj argument") from None

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
        except (AttributeError, NoDataError):
            positions = [None]*ag_or_ts.n_atoms

        try:
            velocities = ag_or_ts.velocities
        except (AttributeError, NoDataError):
            velocities = [None]*ag_or_ts.n_atoms

        atom_kwargs = []
        for atom, name, resname, xyz, vel in zip(ag_or_ts, names, resnames,
                                                 positions, velocities):
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
                    tp = atom.type.lower().capitalize()
                    akwargs['atomic_number'] = SYMB2Z[tp]
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

            atom.atom_type = pmd.AtomType(akwarg['name'], None,
                                          akwarg['mass'],
                                          atomic_number=akwargs.get('atomic_number'))
            struct.add_atom(atom, resname, resid, **kw)

        try:
            struct.box = ag_or_ts.dimensions
        except AttributeError:
            struct.box = None

        if hasattr(ag_or_ts, 'universe'):
            atomgroup = {atom: index for index,
                         atom in enumerate(list(ag_or_ts))}
            get_atom_indices = functools.partial(get_indices_from_subset,
                                                 atomgroup=atomgroup,
                                                 universe=ag_or_ts.universe)
        else:
            get_atom_indices = lambda x: x

        # bonds
        try:
            params = ag_or_ts.intra_bonds
        except AttributeError:
            pass
        else:
            for p in params:
                atoms = [struct.atoms[i] for i in map(get_atom_indices,
                                                      p.indices)]
                try:
                    for obj in p.type:
                        bond = pmd.Bond(*atoms, type=obj.type, order=obj.order)
                        struct.bonds.append(bond)
                    if isinstance(obj.type, pmd.BondType):
                        struct.bond_types.append(bond.type)
                        bond.type.list = struct.bond_types
                except (TypeError, AttributeError):
                    order = p.order if p.order is not None else 1
                    btype = getattr(p.type, 'type', None)

                    bond = pmd.Bond(*atoms, type=btype, order=order)
                    struct.bonds.append(bond)
                    if isinstance(bond.type, pmd.BondType):
                        struct.bond_types.append(bond.type)
                        bond.type.list = struct.bond_types

        # dihedrals
        try:
            params = ag_or_ts.dihedrals.atomgroup_intersection(ag_or_ts,
                                                               strict=True)
        except AttributeError:
            pass
        else:
            for p in params:
                atoms = [struct.atoms[i] for i in map(get_atom_indices,
                                                      p.indices)]
                try:
                    for obj in p.type:
                        imp = getattr(obj, 'improper', False)
                        ign = getattr(obj, 'ignore_end', False)
                        dih = pmd.Dihedral(*atoms, type=obj.type,
                                           ignore_end=ign, improper=imp)
                        struct.dihedrals.append(dih)
                        if isinstance(dih.type, pmd.DihedralType):
                            struct.dihedral_types.append(dih.type)
                            dih.type.list = struct.dihedral_types
                except (TypeError, AttributeError):
                    btype = getattr(p.type, 'type', None)
                    imp = getattr(p.type, 'improper', False)
                    ign = getattr(p.type, 'ignore_end', False)
                    dih = pmd.Dihedral(*atoms, type=btype,
                                       improper=imp, ignore_end=ign)
                    struct.dihedrals.append(dih)
                    if isinstance(dih.type, pmd.DihedralType):
                        struct.dihedral_types.append(dih.type)
                        dih.type.list = struct.dihedral_types

        for param, pmdtype, trackedlist, typelist, clstype in (
            ('ureybradleys', pmd.UreyBradley, struct.urey_bradleys, struct.urey_bradley_types, pmd.BondType),
            ('angles', pmd.Angle, struct.angles, struct.angle_types, pmd.AngleType),
            ('impropers', pmd.Improper, struct.impropers, struct.improper_types, pmd.ImproperType),
            ('cmaps', pmd.Cmap, struct.cmaps, struct.cmap_types, pmd.CmapType)
        ):
            try:
                params = getattr(ag_or_ts, param)
                values = params.atomgroup_intersection(ag_or_ts, strict=True)
            except AttributeError:
                pass
            else:
                for v in values:
                    atoms = [struct.atoms[i] for i in map(get_atom_indices,
                                                          v.indices)]

                    try:
                        for parmed_obj in v.type:
                            p = pmdtype(*atoms, type=parmed_obj.type)
                            trackedlist.append(p)
                            if isinstance(p.type, clstype):
                                typelist.append(p.type)
                                p.type.list = typelist
                    except (TypeError, AttributeError):
                        vtype = getattr(v.type, 'type', None)

                        p = pmdtype(*atoms, type=vtype)
                        trackedlist.append(p)
                        if isinstance(p.type, clstype):
                            typelist.append(p.type)
                            p.type.list = typelist
        return struct
