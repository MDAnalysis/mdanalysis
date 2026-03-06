""" This package contains classes responsible for loading rdkit objects """
from io import StringIO
from typing import TYPE_CHECKING

from ..periodic_table import Element
from ..topologyobjects import QualitativeBondType

class RDKit:

    @staticmethod
    def load(rmol):
        """
        Load a :class:`Mol` object and return a populated :class:`Structure`
        instance

        Parameters
        ----------
        rmol: :class:`Mol`
            RDKit :class:`Mol` object to convert

        Examples
        --------
        >>> from rdkit import Chem
        >>> import parmed as pmd
        >>> mol = Chem.MolFromSmiles('Cc1ccccc1')
        >>> struct = pmd.load_rdkit(mol)
        """
        # TODO - We can convert to a Structure programmatically without having to go through
        # a PDB file first
        from ..formats.pdb import PDBFile
        from rdkit import Chem
        fh = StringIO(Chem.MolToPDBBlock(rmol))
        return PDBFile.parse(fh)

    @staticmethod
    def from_smiles(smiles, coordinates=True, hydrogens=True):
        """
        Load smiles string to :class:`Structure`

        Parameters
        ----------
        smiles : str, smiles
        coordinates : bool, default True
            if True, use `rdkit.Chem.AllChem.EmbedMultipleConfs to assign coordinates
        hydrogens : bool, default True
            if True, use `rdkit.Chem.AddHs` to generate explicit hydrogens

        Returns
        -------
        parm : :class:`Structure`
        """
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles(smiles)

        if hydrogens:
            mol = Chem.AddHs(mol)

        if coordinates:
            AllChem.EmbedMultipleConfs(mol, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)

        parm = RDKit.load(mol)
        if not coordinates:
            parm.coordinates = None
            parm._coordinates = None
        return parm

    @staticmethod
    def from_sdf(filename, structure=False):
        """
        Load SDF file to :class:`Structure`

        Parameters
        ----------
        filename: str
        structure : bool, default False
            if True, return a :class:`Structure`
            if False, return a list of :class:`Structure`
        """
        from rdkit import Chem
        sdf_collection = Chem.SDMolSupplier(filename, removeHs=False)
        if structure:
            mol = next(sdf_collection)
            return RDKit.load(mol)
        else:
            return [RDKit.load(mol) for mol in sdf_collection]

    @classmethod
    def to_mol(cls, structure: "Structure"):
        """ Instantiates an RDKit Mol object from a ParmEd Structure """
        from rdkit.Chem import Atom, RWMol, Conformer, HybridizationType
        from rdkit.Chem.rdchem import BondType

        mol = RWMol()
        for atom in structure.atoms:
            rdatom = Atom(atom.atomic_number)
            if atom.aromatic is not None:
                rdatom.SetIsAromatic(atom.aromatic)
            if atom.formal_charge is not None:
                rdatom.SetFormalCharge(atom.formal_charge)
            if atom.hybridization is not None:
                rdatom.SetHybridization(getattr(HybridizationType, atom.hybridization.name))
            pdb_info = cls._get_pdb_info(atom)
            rdatom.SetMonomerInfo(pdb_info)
            mol.AddAtom(rdatom)

        added_bonds = set()
        for bond in structure.bonds:
            key = frozenset({bond.atom1.idx, bond.atom2.idx})
            if key in added_bonds:
                continue
            added_bonds.add(key)
            bond_type = (
                getattr(BondType, bond.qualitative_type.name)
                if isinstance(bond.qualitative_type, QualitativeBondType)
                else BondType.UNSPECIFIED
            )
            mol.AddBond(bond.atom1.idx, bond.atom2.idx, bond_type)

        coordinates = structure.get_coordinates("all")
        if coordinates is None:
            return mol.GetMol()

        for i in range(coordinates.shape[0]):
            conformer = Conformer(len(structure.atoms))
            for j in range(len(structure.atoms)):
                conformer.SetAtomPosition(j, coordinates[i, j, :])
            mol.AddConformer(conformer)

        return mol.GetMol()

    @classmethod
    def _get_pdb_info(cls, atom: "Atom"):
        from rdkit.Chem.rdchem import AtomPDBResidueInfo
        try:
            residue_number = atom.residue.number if atom.residue.number != -1 else atom.residue.idx + 1
            chain = atom.residue.chain
            insertion_code = atom.residue.insertion_code
        except AttributeError:
            residue_number = 1
            chain = ""
            insertion_code = ""
        return AtomPDBResidueInfo(
            cls._to_pdb_atom_name(atom_name=atom.name, symbol=Element[atom.atomic_number]),
            atom.number if atom.number != -1 else atom.idx + 1,
            atom.altloc,
            atom.residue.name,
            residue_number,
            chain,
            insertion_code,
            atom.occupancy if atom.occupancy != 0.0 else 1.0,
            atom.bfactor,
        )

    @staticmethod
    def _to_pdb_atom_name(*, atom_name: str, symbol: str) -> str:
        """Pad atom_name according to PDB specification.

        For the 13-16 columns relating to the atom name:
        'Alignment of one-letter atom name such as C starts at column 14, while two-letter atom name such as FE starts at column 13.'
        See, http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
        """
        pad_left = " " * (2 - len(symbol)) * int(len(atom_name) != 4)
        return (pad_left + atom_name).ljust(4)


if TYPE_CHECKING:
    from ..structure import Structure
    from ..topologyobjects import Atom
