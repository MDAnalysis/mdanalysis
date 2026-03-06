""" Factory functions that will generate Entos Molecule and QMMMInput from ParmEd Structure objects """
from typing import Optional, Sequence, Union

import numpy as np

from ..structure import Structure
from ..utils.decorators import needs_openmm
from .imports import Molecule, constants, QMMMSystem, HAS_ENTOS, HAS_MISTRAL

def to_entos_molecule(
    struct: Structure,
    charge: Optional[float] = None,
    multiplicity: Optional[int] = None,
) -> Molecule:
    """
    Create an Entos Molecule object from a Structure

    Parameters
    ----------
    struct: :class:`Structure`
        The structure from which a Molecule should be created
    charge: float = None
        The net charge on the molecule. If it is None, the charge is determined as the sum
        of partial atomic charges on each atom in the structure
    multiplicity: int = None
        The spin-multiplicity on the structure. If not set, the lowest multiplicity with the
        given number of electrons will be used

    Returns
    -------
    Molecule
        The molecule given the input structure with coordinates
    """
    if not HAS_ENTOS:
        raise ImportError("to_entos_molecule will not work without installing Entos Sierra")

    if struct.coordinates is None:
        raise ValueError("Cannot create an Entos Molecule without coordinates")

    if charge is None:
        charge = sum(a.charge for a in struct.atoms)

    return Molecule(
        atomic_numbers=[a.atomic_number for a in struct.atoms],
        geometry=struct.coordinates * constants.cf("angstrom", "bohr"),
        multiplicity=multiplicity,
        charge=charge,
    )

@needs_openmm
def to_entos_qmmm_system(
    struct: Structure,
    qm_selection: Union[str, Sequence[int]],
    charge: Optional[float] = None,
    multiplicity: Optional[int] = None,
    **kwargs,
) -> QMMMSystem:
    """
    Create a multiscale input for QM/MM simulations with Qcore from a Structure object

    Parameters
    ----------
    struct: :class:`Structure`
        The input structure to create the multiscale input from
    qm_selection: :class:`Sequence` or str
        The selection array (starting from 0), mask array (if length matches
        struct atom count), or atom mask syntax (if selection is a string) for the
        QM atoms in the system
    charge: float = None
        The charge on the QM region. If not set, the value is taken from the sum of
        the charges on the qm_selection
    multiplicity: int = None
        The multiplicity on the QM region. The default is taken as the minimum multiplicity
        given the number of electrons in the QM region.
    **kwargs
        Additional keyword arguments are passed to the `createSystem` call when creating an OpenMM
        System object for doing QM/MM calculations

    Returns
    -------
    multiscale_input: QMMMInput
        The QM/MM input
    """
    from openmm.app import LJPME, NoCutoff
    if not HAS_MISTRAL:
        raise ImportError("You must install Entos mistral to create a QMMMInput object")
    if struct.coordinates is None:
        raise ValueError("Cannot create a QMMMSystem with no coordinates")
    selection_array = _selection_to_array(struct, qm_selection)
    if selection_array[0] < 0 or selection_array[-1] >= len(struct.atoms):
        raise ValueError("QM selection array is out of bounds")
    create_system_args = {"nonbondedMethod": NoCutoff if struct.box is None else LJPME}
    create_system_args.update(kwargs)
    system = struct.createSystem(**create_system_args)
    if charge is None:
        charge = sum(struct.atoms[i].charge for i in selection_array)
    return QMMMSystem(
        geometry=struct.coordinates * constants.cf("angstrom", "bohr"),
        atomic_numbers=[a.atomic_number for a in struct.atoms],
        forcefield=system,
        topology=struct.topology,
        qm_indices=selection_array,
        qm_charge=charge,
        qm_multiplicity=multiplicity,
    )

def _selection_to_array(struct: Structure, qm_selection: Union[str, Sequence[int]]) -> np.ndarray:
    """ Converts a selection into an array of indices """
    if isinstance(qm_selection, str):
        from ..amber import AmberMask
        mask = AmberMask(struct, qm_selection)
        return np.array(list(mask.Selected()))
    elif len(qm_selection) == len(struct.atoms):
        return np.array([i for i, mask in enumerate(qm_selection) if mask])
    else:
        return np.array(sorted(list(set(qm_selection))))
