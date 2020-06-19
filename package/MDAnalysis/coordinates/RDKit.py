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

"""RDKit molecule --- :mod:`MDAnalysis.coordinates.RDKit`
================================================================

Read coordinates data from an `RDKit <https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html#rdkit.Chem.rdchem.Mol>`_ :class:`rdkit.Chem.rdchem.Mol` with :class:`RDKitReader` 
into a MDAnalysis Universe. Convert it back to a :class:`rdkit.Chem.rdchem.Mol` with 
:class:`RDKitConverter`.


Example
-------

>>> from rdkit import Chem
>>> import MDAnalysis as mda
>>> mol = Chem.MolFromMol2File("docking_poses.mol2", removeHs=False)
>>> u = mda.Universe(mol)
>>> u
<Universe with 42 atoms>
>>> u.trajectory
<RDKitReader with 10 frames of 42 atoms>


Classes
-------

.. autoclass:: RDKitReader
   :members:

.. autoclass:: RDKitConverter
   :members:


"""

import warnings

import numpy as np

from . import memory
from . import base

try:
    from rdkit import Chem
except ImportError:
    pass
else:
    RDBONDTYPE = {
        'AROMATIC': Chem.BondType.AROMATIC,
        'SINGLE': Chem.BondType.SINGLE,
        'DOUBLE': Chem.BondType.DOUBLE,
        'TRIPLE': Chem.BondType.TRIPLE,
    }
    RDBONDORDER = {
        1: Chem.BondType.SINGLE,
        1.5: Chem.BondType.AROMATIC,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE,
    }

class RDKitReader(memory.MemoryReader):
    """Coordinate reader for RDKit.
    
    .. versionadded:: 2.0.0
    """
    format = 'RDKIT'

    # Structure.coordinates always in Angstrom
    units = {'time': None, 'length': 'Angstrom'}

    @staticmethod
    def _format_hint(thing):
        """Can this reader read *thing*?"""
        try:
            from rdkit import Chem
        except ImportError:
            # if we can't import rdkit, it's probably not rdkit
            return False
        else:
            return isinstance(thing, Chem.Mol)

    def __init__(self, filename, **kwargs):
        """Read coordinates from an RDKit molecule.
        Each conformer in the original RDKit molecule will be read as a frame 
        in the resulting universe.

        Parameters
        ----------

        filename : rdkit.Chem.rdchem.Mol
            RDKit molecule
        """
        n_atoms = filename.GetNumAtoms()
        coordinates = np.array([
            conf.GetPositions() for conf in filename.GetConformers()], 
            dtype=np.float32)
        if coordinates.size == 0:
            warnings.warn("No coordinates found in the RDKit molecule")
            coordinates = np.empty((1,n_atoms,3), dtype=np.float32)
            coordinates[:] = np.nan
        super(RDKitReader, self).__init__(coordinates, order='fac', **kwargs)


class RDKitConverter(base.ConverterBase):
    """Convert MDAnalysis AtomGroup or Universe to `RDKit <https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html#rdkit.Chem.rdchem.Mol>`_ :class:`rdkit.Chem.rdchem.Mol`.

    Example
    -------

    .. code-block:: python

        import MDAnalysis as mda
        from MDAnalysis.tests.datafiles import PDB_full
        u = mda.Universe(PDB_full)
        mol = u.select_atoms('resname DMS').convert_to('RDKIT')
    

    .. versionadded:: 2.X.X
    """

    lib = 'RDKIT'
    units = {'time': None, 'length': 'Angstrom'}

    def convert(self, obj):
        """Write selection at current trajectory frame to :class:`~rdkit.Chem.rdchem.Mol`.

        Parameters
        -----------
        obj : AtomGroup or Universe or :class:`Timestep`
        """
        try:
            from rdkit import Chem
        except ImportError:
            raise ImportError('RDKit is required for RDKitConverter but '
                              'is not installed. Try installing it with \n'
                              'conda install -c conda-forge rdkit')
        try:
            # make sure to use atoms (Issue 46)
            ag_or_ts = obj.atoms
        except AttributeError as e:
            if isinstance(obj, base.Timestep):
                ag_or_ts = obj.copy()
            else:
                raise TypeError("No Timestep found in obj argument") from e

        mol = Chem.RWMol()
        atom_mapper = {}
        for atom in ag_or_ts:
            rdatom = Chem.Atom(atom.element)
            index = mol.AddAtom(rdatom)
            atom_mapper[atom.ix] = index

        for bond in ag_or_ts.bonds:
            bond_indices = [atom_mapper[i] for i in bond.indices]
            bond_type = RDBONDTYPE.get(bond.type.upper(), RDBONDORDER.get(
                bond.order, Chem.BondType.SINGLE))
            mol.AddBond(*bond_indices, bond_type)

        Chem.SanitizeMol(mol)
        return mol