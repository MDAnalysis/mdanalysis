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

"""ParmEd structure I/O --- :mod:`MDAnalysis.coordinates.ParmEd`
================================================================

Read coordinates data from a ParmEd Structure with :class:`ParmEdReader` into a MDAnalysis Universe. Convert it back to a ParmEd Structure with :class:`ParmEdConverter`.

Example
-------

::

    import parmed as pmd
    import MDAnalysis as mda
    from MDAnalysis.tests.datafiles import GRO
    universe = mda.Universe(pmd.load_file(GRO))
    structure = universe.atoms.convert_to('PARMED')


Classes
-------

.. autoclass:: ParmEdReader
   :members:

.. autoclass:: ParmEdConverter
   :members:


"""
from __future__ import absolute_import

from . import base
from ..topology.tables import SYMB2Z
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
        parmed_subset = mgro.select_atoms('resname SOL').to_format('PARMED')

        
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
            import parmed as pmd
        except ModuleNotFoundError:
            raise ValueError('Parmed is required for ParmEdConverter but '
                             'is not installed.')
        try:
            # make sure to use atoms (Issue 46)
            ag_or_ts = obj.atoms
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
            struct.add_atom(atom, resname, resid, **kw)
        
        try:
            struct.box = ag_or_ts.dimensions
        except AttributeError:
            struct.box = None
        
        
        try:
            print('trying params')
            params = ag_or_ts.bonds
            print(params)
        except AttributeError:
            pass
        else:
            for p in params:
                atoms = [struct.atoms[i] for i in p.indices]
                try:
                    print('pre trying', p.type)
                    for obj in p.type:
                        bond = pmd.Bond(*atoms, type=obj.type, order=obj.order)
                        struct.bonds.append(bond)
                        print(obj, 'trying', bond, len(struct.bonds))
                except (TypeError, AttributeError):
                    order = p.order if p.order is not None else 1
                    bond = pmd.Bond(*atoms, order=order)
                    struct.bonds.append(bond)
                    print(bond)
                    print(len(struct.bonds))
        
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
                pass
            else:
                for v in values:
                    atoms = [struct.atoms[i] for i in v.indices]

                    try:
                        for parmed_obj in v.type:
                            p = pmdtype(*atoms, type=parmed_obj.type)
                            trackedlist.append(p)
                    except (TypeError, AttributeError):
                        p = pmdtype(*atoms)
                        trackedlist.append(p)
        return struct