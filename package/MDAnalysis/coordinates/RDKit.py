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
import re

import numpy as np

from ..exceptions import NoDataError
from ..topology.guessers import guess_atom_element
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
        "ar": Chem.BondType.AROMATIC,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE,
    }
    # add string version of the key for each bond
    RDBONDORDER.update({str(key):value for key,value in RDBONDORDER.items()})
    RDATTRIBUTES = {
        "altLoc": "AltLoc",
        "chainID": "ChainId",
        "name": "Name",
        "occupancy": "Occupancy",
        "resname": "ResidueName",
        "resid": "ResidueNumber",
        "segid": "SegmentNumber",
        "tempfactor": "TempFactor",
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
        obj : AtomGroup or Universe
        """
        try:
            from rdkit import Chem
        except ImportError:
            raise ImportError('RDKit is required for RDKitConverter but '
                              'is not installed. Try installing it with \n'
                              'conda install -c conda-forge rdkit')
        try:
            # make sure to use atoms (Issue 46)
            ag = obj.atoms
        except AttributeError as e:
            raise TypeError("No `atoms` attribute in object of type {}, "
                            "please use a valid AtomGroup or Universe".format(
                            type(obj))) from e

        mol = Chem.RWMol()
        atom_mapper = {}

        for atom in ag:
            try:
                element = atom.element
            except NoDataError:
                # guess atom element
                # capitalize: transform CL to Cl and so on
                element = guess_atom_element(atom.name).capitalize()
            rdatom = Chem.Atom(element)
            # add properties
            mi = Chem.AtomPDBResidueInfo()
            for attr, rdattr in RDATTRIBUTES.items():
                try: # get value in MDA atom
                    value = getattr(atom, attr)
                except AttributeError:
                    pass
                else:
                    if isinstance(value, np.generic):
                        # convert numpy types to python standard types
                        value = value.item()
                    if attr == "segid":
                        # RDKit needs segid to be an int
                        try:
                            value = int(value)
                        except ValueError:
                            # convert any string to int
                            # can be mapped back with np.base_repr(x, 36)
                            value = int(value, 36)
                    elif attr == "name":
                        # RDKit needs the name to be properly formated for a
                        # PDB file (1 letter elements start at col 14)
                        name = re.findall('(\D+|\d+)', value)
                        if len(name) == 2:
                            symbol, number = name
                        else:
                            symbol, number = name[0], ""
                        value = "{:>2}".format(symbol) + "{:<2}".format(number)
                    # set attribute value in RDKit MonomerInfo
                    getattr(mi, "Set%s" % rdattr)(value)
            rdatom.SetMonomerInfo(mi)
            # TODO other properties (charges)
            index = mol.AddAtom(rdatom)
            # map index in universe to index in mol
            atom_mapper[atom.ix] = index

        try:
            bonds = ag.bonds
        except NoDataError:
            ag.guess_bonds()
            bonds = ag.bonds

        for bond in bonds:
            bond_indices = [atom_mapper[i] for i in bond.indices]
            try:
                bond_type = bond.type.upper()
            except AttributeError:
                # bond type can be a tuple for PDB files
                bond_type = None
            bond_type = RDBONDTYPE.get(bond_type, RDBONDORDER.get(
                bond.order, Chem.BondType.SINGLE))
            mol.AddBond(*bond_indices, bond_type)

        Chem.SanitizeMol(mol)
        return mol