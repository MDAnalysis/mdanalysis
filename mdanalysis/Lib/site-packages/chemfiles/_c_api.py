# -* coding: utf-8 -*
# Chemfiles, a modern library for chemistry file reading and writing
# Copyright (C) Guillaume Fraux and contributors -- BSD license
#
# =========================================================================== #
# !!!! AUTO-GENERATED FILE !!!! Do not edit. See the bindgen repository for
# the generation code (https://github.com/chemfiles/bindgen).
# This file contains Python ctype interface to the C API
# =========================================================================== #

# flake8: noqa
'''
Foreign function interface declaration for the Python interface to chemfiles
'''
from numpy.ctypeslib import ndpointer
import numpy as np
from ctypes import c_int, c_int64, c_uint64, c_double, c_char, c_char_p, c_void_p, c_bool
from ctypes import CFUNCTYPE, ARRAY, POINTER, Structure

from .utils import _check_return_code


class chfl_status(c_int):
    CHFL_SUCCESS = 0
    CHFL_MEMORY_ERROR = 1
    CHFL_FILE_ERROR = 2
    CHFL_FORMAT_ERROR = 3
    CHFL_SELECTION_ERROR = 4
    CHFL_CONFIGURATION_ERROR = 5
    CHFL_OUT_OF_BOUNDS = 6
    CHFL_PROPERTY_ERROR = 7
    CHFL_GENERIC_ERROR = 254
    CHFL_CXX_ERROR = 255


class chfl_bond_order(c_int):
    CHFL_BOND_UNKNOWN = 0
    CHFL_BOND_SINGLE = 1
    CHFL_BOND_DOUBLE = 2
    CHFL_BOND_TRIPLE = 3
    CHFL_BOND_QUADRUPLE = 4
    CHFL_BOND_QUINTUPLET = 5
    CHFL_BOND_AMIDE = 254
    CHFL_BOND_AROMATIC = 255


class chfl_property_kind(c_int):
    CHFL_PROPERTY_BOOL = 0
    CHFL_PROPERTY_DOUBLE = 1
    CHFL_PROPERTY_STRING = 2
    CHFL_PROPERTY_VECTOR3D = 3


class chfl_cellshape(c_int):
    CHFL_CELL_ORTHORHOMBIC = 0
    CHFL_CELL_TRICLINIC = 1
    CHFL_CELL_INFINITE = 2


class CHFL_TRAJECTORY(Structure):
    pass


class CHFL_CELL(Structure):
    pass


class CHFL_ATOM(Structure):
    pass


class CHFL_FRAME(Structure):
    pass


class CHFL_TOPOLOGY(Structure):
    pass


class CHFL_SELECTION(Structure):
    pass


class CHFL_RESIDUE(Structure):
    pass


class CHFL_PROPERTY(Structure):
    pass

# Some hand-defined type. Make sure to edit the bindgen code to make this
# correspond to the current chemfiles.h header
chfl_vector3d = ARRAY(c_double, 3)

chfl_warning_callback = CFUNCTYPE(None, c_char_p)

class chfl_match(Structure):
    _fields_ = [
        ('size', c_uint64),
        ('atoms', ARRAY(c_uint64, 4))
    ]

class chfl_format_metadata(Structure):
    _fields_ = [
        ('name', c_char_p),
        ('extension', c_char_p),
        ('description', c_char_p),
        ('reference', c_char_p),

        ('read', c_bool),
        ('write', c_bool),
        ('memory', c_bool),
        ('positions', c_bool),
        ('velocities', c_bool),
        ('unit_cell', c_bool),
        ('atoms', c_bool),
        ('bonds', c_bool),
        ('residues', c_bool),
    ]

# end of hand-defined types


def set_interface(c_lib):
    # Manually defined functions
    c_lib.chfl_free.argtypes = [c_void_p]
    c_lib.chfl_trajectory_close.argtypes = [POINTER(CHFL_TRAJECTORY)]
    # End of manually defined functions

    # Function "chfl_version", at types.h:150
    c_lib.chfl_version.argtypes = []
    c_lib.chfl_version.restype = c_char_p

    # Function "chfl_last_error", at misc.h:23
    c_lib.chfl_last_error.argtypes = []
    c_lib.chfl_last_error.restype = c_char_p

    # Function "chfl_clear_errors", at misc.h:33
    c_lib.chfl_clear_errors.argtypes = []
    c_lib.chfl_clear_errors.restype = chfl_status
    c_lib.chfl_clear_errors.errcheck = _check_return_code

    # Function "chfl_set_warning_callback", at misc.h:42
    c_lib.chfl_set_warning_callback.argtypes = [chfl_warning_callback]
    c_lib.chfl_set_warning_callback.restype = chfl_status
    c_lib.chfl_set_warning_callback.errcheck = _check_return_code

    # Function "chfl_add_configuration", at misc.h:58
    c_lib.chfl_add_configuration.argtypes = [c_char_p]
    c_lib.chfl_add_configuration.restype = chfl_status
    c_lib.chfl_add_configuration.errcheck = _check_return_code

    # Function "chfl_formats_list", at misc.h:71
    c_lib.chfl_formats_list.argtypes = [POINTER(POINTER(chfl_format_metadata)), POINTER(c_uint64)]
    c_lib.chfl_formats_list.restype = chfl_status
    c_lib.chfl_formats_list.errcheck = _check_return_code

    # Function "chfl_guess_format", at misc.h:93
    c_lib.chfl_guess_format.argtypes = [c_char_p, c_char_p, c_uint64]
    c_lib.chfl_guess_format.restype = chfl_status
    c_lib.chfl_guess_format.errcheck = _check_return_code

    # Function "chfl_property_bool", at property.h:37
    c_lib.chfl_property_bool.argtypes = [c_bool]
    c_lib.chfl_property_bool.restype = POINTER(CHFL_PROPERTY)

    # Function "chfl_property_double", at property.h:47
    c_lib.chfl_property_double.argtypes = [c_double]
    c_lib.chfl_property_double.restype = POINTER(CHFL_PROPERTY)

    # Function "chfl_property_string", at property.h:57
    c_lib.chfl_property_string.argtypes = [c_char_p]
    c_lib.chfl_property_string.restype = POINTER(CHFL_PROPERTY)

    # Function "chfl_property_vector3d", at property.h:67
    c_lib.chfl_property_vector3d.argtypes = [chfl_vector3d]
    c_lib.chfl_property_vector3d.restype = POINTER(CHFL_PROPERTY)

    # Function "chfl_property_get_kind", at property.h:74
    c_lib.chfl_property_get_kind.argtypes = [POINTER(CHFL_PROPERTY), POINTER(chfl_property_kind)]
    c_lib.chfl_property_get_kind.restype = chfl_status
    c_lib.chfl_property_get_kind.errcheck = _check_return_code

    # Function "chfl_property_get_bool", at property.h:87
    c_lib.chfl_property_get_bool.argtypes = [POINTER(CHFL_PROPERTY), POINTER(c_bool)]
    c_lib.chfl_property_get_bool.restype = chfl_status
    c_lib.chfl_property_get_bool.errcheck = _check_return_code

    # Function "chfl_property_get_double", at property.h:100
    c_lib.chfl_property_get_double.argtypes = [POINTER(CHFL_PROPERTY), POINTER(c_double)]
    c_lib.chfl_property_get_double.restype = chfl_status
    c_lib.chfl_property_get_double.errcheck = _check_return_code

    # Function "chfl_property_get_string", at property.h:115
    c_lib.chfl_property_get_string.argtypes = [POINTER(CHFL_PROPERTY), c_char_p, c_uint64]
    c_lib.chfl_property_get_string.restype = chfl_status
    c_lib.chfl_property_get_string.errcheck = _check_return_code

    # Function "chfl_property_get_vector3d", at property.h:128
    c_lib.chfl_property_get_vector3d.argtypes = [POINTER(CHFL_PROPERTY), chfl_vector3d]
    c_lib.chfl_property_get_vector3d.restype = chfl_status
    c_lib.chfl_property_get_vector3d.errcheck = _check_return_code

    # Function "chfl_atom", at atom.h:24
    c_lib.chfl_atom.argtypes = [c_char_p]
    c_lib.chfl_atom.restype = POINTER(CHFL_ATOM)

    # Function "chfl_atom_copy", at atom.h:34
    c_lib.chfl_atom_copy.argtypes = [POINTER(CHFL_ATOM)]
    c_lib.chfl_atom_copy.restype = POINTER(CHFL_ATOM)

    # Function "chfl_atom_from_frame", at atom.h:60
    c_lib.chfl_atom_from_frame.argtypes = [POINTER(CHFL_FRAME), c_uint64]
    c_lib.chfl_atom_from_frame.restype = POINTER(CHFL_ATOM)

    # Function "chfl_atom_from_topology", at atom.h:83
    c_lib.chfl_atom_from_topology.argtypes = [POINTER(CHFL_TOPOLOGY), c_uint64]
    c_lib.chfl_atom_from_topology.restype = POINTER(CHFL_ATOM)

    # Function "chfl_atom_mass", at atom.h:94
    c_lib.chfl_atom_mass.argtypes = [POINTER(CHFL_ATOM), POINTER(c_double)]
    c_lib.chfl_atom_mass.restype = chfl_status
    c_lib.chfl_atom_mass.errcheck = _check_return_code

    # Function "chfl_atom_set_mass", at atom.h:103
    c_lib.chfl_atom_set_mass.argtypes = [POINTER(CHFL_ATOM), c_double]
    c_lib.chfl_atom_set_mass.restype = chfl_status
    c_lib.chfl_atom_set_mass.errcheck = _check_return_code

    # Function "chfl_atom_charge", at atom.h:112
    c_lib.chfl_atom_charge.argtypes = [POINTER(CHFL_ATOM), POINTER(c_double)]
    c_lib.chfl_atom_charge.restype = chfl_status
    c_lib.chfl_atom_charge.errcheck = _check_return_code

    # Function "chfl_atom_set_charge", at atom.h:121
    c_lib.chfl_atom_set_charge.argtypes = [POINTER(CHFL_ATOM), c_double]
    c_lib.chfl_atom_set_charge.restype = chfl_status
    c_lib.chfl_atom_set_charge.errcheck = _check_return_code

    # Function "chfl_atom_type", at atom.h:131
    c_lib.chfl_atom_type.argtypes = [POINTER(CHFL_ATOM), c_char_p, c_uint64]
    c_lib.chfl_atom_type.restype = chfl_status
    c_lib.chfl_atom_type.errcheck = _check_return_code

    # Function "chfl_atom_set_type", at atom.h:142
    c_lib.chfl_atom_set_type.argtypes = [POINTER(CHFL_ATOM), c_char_p]
    c_lib.chfl_atom_set_type.restype = chfl_status
    c_lib.chfl_atom_set_type.errcheck = _check_return_code

    # Function "chfl_atom_name", at atom.h:152
    c_lib.chfl_atom_name.argtypes = [POINTER(CHFL_ATOM), c_char_p, c_uint64]
    c_lib.chfl_atom_name.restype = chfl_status
    c_lib.chfl_atom_name.errcheck = _check_return_code

    # Function "chfl_atom_set_name", at atom.h:163
    c_lib.chfl_atom_set_name.argtypes = [POINTER(CHFL_ATOM), c_char_p]
    c_lib.chfl_atom_set_name.restype = chfl_status
    c_lib.chfl_atom_set_name.errcheck = _check_return_code

    # Function "chfl_atom_full_name", at atom.h:173
    c_lib.chfl_atom_full_name.argtypes = [POINTER(CHFL_ATOM), c_char_p, c_uint64]
    c_lib.chfl_atom_full_name.restype = chfl_status
    c_lib.chfl_atom_full_name.errcheck = _check_return_code

    # Function "chfl_atom_vdw_radius", at atom.h:185
    c_lib.chfl_atom_vdw_radius.argtypes = [POINTER(CHFL_ATOM), POINTER(c_double)]
    c_lib.chfl_atom_vdw_radius.restype = chfl_status
    c_lib.chfl_atom_vdw_radius.errcheck = _check_return_code

    # Function "chfl_atom_covalent_radius", at atom.h:195
    c_lib.chfl_atom_covalent_radius.argtypes = [POINTER(CHFL_ATOM), POINTER(c_double)]
    c_lib.chfl_atom_covalent_radius.restype = chfl_status
    c_lib.chfl_atom_covalent_radius.errcheck = _check_return_code

    # Function "chfl_atom_atomic_number", at atom.h:205
    c_lib.chfl_atom_atomic_number.argtypes = [POINTER(CHFL_ATOM), POINTER(c_uint64)]
    c_lib.chfl_atom_atomic_number.restype = chfl_status
    c_lib.chfl_atom_atomic_number.errcheck = _check_return_code

    # Function "chfl_atom_properties_count", at atom.h:212
    c_lib.chfl_atom_properties_count.argtypes = [POINTER(CHFL_ATOM), POINTER(c_uint64)]
    c_lib.chfl_atom_properties_count.restype = chfl_status
    c_lib.chfl_atom_properties_count.errcheck = _check_return_code

    # Function "chfl_atom_list_properties", at atom.h:228
    c_lib.chfl_atom_list_properties.argtypes = [POINTER(CHFL_ATOM), POINTER(c_char_p), c_uint64]
    c_lib.chfl_atom_list_properties.restype = chfl_status
    c_lib.chfl_atom_list_properties.errcheck = _check_return_code

    # Function "chfl_atom_set_property", at atom.h:240
    c_lib.chfl_atom_set_property.argtypes = [POINTER(CHFL_ATOM), c_char_p, POINTER(CHFL_PROPERTY)]
    c_lib.chfl_atom_set_property.restype = chfl_status
    c_lib.chfl_atom_set_property.errcheck = _check_return_code

    # Function "chfl_atom_get_property", at atom.h:254
    c_lib.chfl_atom_get_property.argtypes = [POINTER(CHFL_ATOM), c_char_p]
    c_lib.chfl_atom_get_property.restype = POINTER(CHFL_PROPERTY)

    # Function "chfl_residue", at residue.h:25
    c_lib.chfl_residue.argtypes = [c_char_p]
    c_lib.chfl_residue.restype = POINTER(CHFL_RESIDUE)

    # Function "chfl_residue_with_id", at residue.h:35
    c_lib.chfl_residue_with_id.argtypes = [c_char_p, c_int64]
    c_lib.chfl_residue_with_id.restype = POINTER(CHFL_RESIDUE)

    # Function "chfl_residue_from_topology", at residue.h:63
    c_lib.chfl_residue_from_topology.argtypes = [POINTER(CHFL_TOPOLOGY), c_uint64]
    c_lib.chfl_residue_from_topology.restype = POINTER(CHFL_RESIDUE)

    # Function "chfl_residue_for_atom", at residue.h:90
    c_lib.chfl_residue_for_atom.argtypes = [POINTER(CHFL_TOPOLOGY), c_uint64]
    c_lib.chfl_residue_for_atom.restype = POINTER(CHFL_RESIDUE)

    # Function "chfl_residue_copy", at residue.h:102
    c_lib.chfl_residue_copy.argtypes = [POINTER(CHFL_RESIDUE)]
    c_lib.chfl_residue_copy.restype = POINTER(CHFL_RESIDUE)

    # Function "chfl_residue_atoms_count", at residue.h:109
    c_lib.chfl_residue_atoms_count.argtypes = [POINTER(CHFL_RESIDUE), POINTER(c_uint64)]
    c_lib.chfl_residue_atoms_count.restype = chfl_status
    c_lib.chfl_residue_atoms_count.errcheck = _check_return_code

    # Function "chfl_residue_atoms", at residue.h:123
    c_lib.chfl_residue_atoms.argtypes = [POINTER(CHFL_RESIDUE), ndpointer(np.uint64, flags="C_CONTIGUOUS", ndim=1), c_uint64]
    c_lib.chfl_residue_atoms.restype = chfl_status
    c_lib.chfl_residue_atoms.errcheck = _check_return_code

    # Function "chfl_residue_id", at residue.h:136
    c_lib.chfl_residue_id.argtypes = [POINTER(CHFL_RESIDUE), POINTER(c_int64)]
    c_lib.chfl_residue_id.restype = chfl_status
    c_lib.chfl_residue_id.errcheck = _check_return_code

    # Function "chfl_residue_name", at residue.h:148
    c_lib.chfl_residue_name.argtypes = [POINTER(CHFL_RESIDUE), c_char_p, c_uint64]
    c_lib.chfl_residue_name.restype = chfl_status
    c_lib.chfl_residue_name.errcheck = _check_return_code

    # Function "chfl_residue_add_atom", at residue.h:157
    c_lib.chfl_residue_add_atom.argtypes = [POINTER(CHFL_RESIDUE), c_uint64]
    c_lib.chfl_residue_add_atom.restype = chfl_status
    c_lib.chfl_residue_add_atom.errcheck = _check_return_code

    # Function "chfl_residue_contains", at residue.h:167
    c_lib.chfl_residue_contains.argtypes = [POINTER(CHFL_RESIDUE), c_uint64, POINTER(c_bool)]
    c_lib.chfl_residue_contains.restype = chfl_status
    c_lib.chfl_residue_contains.errcheck = _check_return_code

    # Function "chfl_residue_properties_count", at residue.h:176
    c_lib.chfl_residue_properties_count.argtypes = [POINTER(CHFL_RESIDUE), POINTER(c_uint64)]
    c_lib.chfl_residue_properties_count.restype = chfl_status
    c_lib.chfl_residue_properties_count.errcheck = _check_return_code

    # Function "chfl_residue_list_properties", at residue.h:192
    c_lib.chfl_residue_list_properties.argtypes = [POINTER(CHFL_RESIDUE), POINTER(c_char_p), c_uint64]
    c_lib.chfl_residue_list_properties.restype = chfl_status
    c_lib.chfl_residue_list_properties.errcheck = _check_return_code

    # Function "chfl_residue_set_property", at residue.h:204
    c_lib.chfl_residue_set_property.argtypes = [POINTER(CHFL_RESIDUE), c_char_p, POINTER(CHFL_PROPERTY)]
    c_lib.chfl_residue_set_property.restype = chfl_status
    c_lib.chfl_residue_set_property.errcheck = _check_return_code

    # Function "chfl_residue_get_property", at residue.h:218
    c_lib.chfl_residue_get_property.argtypes = [POINTER(CHFL_RESIDUE), c_char_p]
    c_lib.chfl_residue_get_property.restype = POINTER(CHFL_PROPERTY)

    # Function "chfl_topology", at topology.h:25
    c_lib.chfl_topology.argtypes = []
    c_lib.chfl_topology.restype = POINTER(CHFL_TOPOLOGY)

    # Function "chfl_topology_from_frame", at topology.h:38
    c_lib.chfl_topology_from_frame.argtypes = [POINTER(CHFL_FRAME)]
    c_lib.chfl_topology_from_frame.restype = POINTER(CHFL_TOPOLOGY)

    # Function "chfl_topology_copy", at topology.h:48
    c_lib.chfl_topology_copy.argtypes = [POINTER(CHFL_TOPOLOGY)]
    c_lib.chfl_topology_copy.restype = POINTER(CHFL_TOPOLOGY)

    # Function "chfl_topology_atoms_count", at topology.h:56
    c_lib.chfl_topology_atoms_count.argtypes = [POINTER(CHFL_TOPOLOGY), POINTER(c_uint64)]
    c_lib.chfl_topology_atoms_count.restype = chfl_status
    c_lib.chfl_topology_atoms_count.errcheck = _check_return_code

    # Function "chfl_topology_resize", at topology.h:68
    c_lib.chfl_topology_resize.argtypes = [POINTER(CHFL_TOPOLOGY), c_uint64]
    c_lib.chfl_topology_resize.restype = chfl_status
    c_lib.chfl_topology_resize.errcheck = _check_return_code

    # Function "chfl_topology_add_atom", at topology.h:75
    c_lib.chfl_topology_add_atom.argtypes = [POINTER(CHFL_TOPOLOGY), POINTER(CHFL_ATOM)]
    c_lib.chfl_topology_add_atom.restype = chfl_status
    c_lib.chfl_topology_add_atom.errcheck = _check_return_code

    # Function "chfl_topology_remove", at topology.h:86
    c_lib.chfl_topology_remove.argtypes = [POINTER(CHFL_TOPOLOGY), c_uint64]
    c_lib.chfl_topology_remove.restype = chfl_status
    c_lib.chfl_topology_remove.errcheck = _check_return_code

    # Function "chfl_topology_bonds_count", at topology.h:95
    c_lib.chfl_topology_bonds_count.argtypes = [POINTER(CHFL_TOPOLOGY), POINTER(c_uint64)]
    c_lib.chfl_topology_bonds_count.restype = chfl_status
    c_lib.chfl_topology_bonds_count.errcheck = _check_return_code

    # Function "chfl_topology_angles_count", at topology.h:104
    c_lib.chfl_topology_angles_count.argtypes = [POINTER(CHFL_TOPOLOGY), POINTER(c_uint64)]
    c_lib.chfl_topology_angles_count.restype = chfl_status
    c_lib.chfl_topology_angles_count.errcheck = _check_return_code

    # Function "chfl_topology_dihedrals_count", at topology.h:113
    c_lib.chfl_topology_dihedrals_count.argtypes = [POINTER(CHFL_TOPOLOGY), POINTER(c_uint64)]
    c_lib.chfl_topology_dihedrals_count.restype = chfl_status
    c_lib.chfl_topology_dihedrals_count.errcheck = _check_return_code

    # Function "chfl_topology_impropers_count", at topology.h:122
    c_lib.chfl_topology_impropers_count.argtypes = [POINTER(CHFL_TOPOLOGY), POINTER(c_uint64)]
    c_lib.chfl_topology_impropers_count.restype = chfl_status
    c_lib.chfl_topology_impropers_count.errcheck = _check_return_code

    # Function "chfl_topology_bonds", at topology.h:135
    c_lib.chfl_topology_bonds.argtypes = [POINTER(CHFL_TOPOLOGY), ndpointer(np.uint64, flags="C_CONTIGUOUS", ndim=2), c_uint64]
    c_lib.chfl_topology_bonds.restype = chfl_status
    c_lib.chfl_topology_bonds.errcheck = _check_return_code

    # Function "chfl_topology_angles", at topology.h:148
    c_lib.chfl_topology_angles.argtypes = [POINTER(CHFL_TOPOLOGY), ndpointer(np.uint64, flags="C_CONTIGUOUS", ndim=2), c_uint64]
    c_lib.chfl_topology_angles.restype = chfl_status
    c_lib.chfl_topology_angles.errcheck = _check_return_code

    # Function "chfl_topology_dihedrals", at topology.h:162
    c_lib.chfl_topology_dihedrals.argtypes = [POINTER(CHFL_TOPOLOGY), ndpointer(np.uint64, flags="C_CONTIGUOUS", ndim=2), c_uint64]
    c_lib.chfl_topology_dihedrals.restype = chfl_status
    c_lib.chfl_topology_dihedrals.errcheck = _check_return_code

    # Function "chfl_topology_impropers", at topology.h:176
    c_lib.chfl_topology_impropers.argtypes = [POINTER(CHFL_TOPOLOGY), ndpointer(np.uint64, flags="C_CONTIGUOUS", ndim=2), c_uint64]
    c_lib.chfl_topology_impropers.restype = chfl_status
    c_lib.chfl_topology_impropers.errcheck = _check_return_code

    # Function "chfl_topology_add_bond", at topology.h:185
    c_lib.chfl_topology_add_bond.argtypes = [POINTER(CHFL_TOPOLOGY), c_uint64, c_uint64]
    c_lib.chfl_topology_add_bond.restype = chfl_status
    c_lib.chfl_topology_add_bond.errcheck = _check_return_code

    # Function "chfl_topology_remove_bond", at topology.h:197
    c_lib.chfl_topology_remove_bond.argtypes = [POINTER(CHFL_TOPOLOGY), c_uint64, c_uint64]
    c_lib.chfl_topology_remove_bond.restype = chfl_status
    c_lib.chfl_topology_remove_bond.errcheck = _check_return_code

    # Function "chfl_topology_clear_bonds", at topology.h:207
    c_lib.chfl_topology_clear_bonds.argtypes = [POINTER(CHFL_TOPOLOGY)]
    c_lib.chfl_topology_clear_bonds.restype = chfl_status
    c_lib.chfl_topology_clear_bonds.errcheck = _check_return_code

    # Function "chfl_topology_residues_count", at topology.h:215
    c_lib.chfl_topology_residues_count.argtypes = [POINTER(CHFL_TOPOLOGY), POINTER(c_uint64)]
    c_lib.chfl_topology_residues_count.restype = chfl_status
    c_lib.chfl_topology_residues_count.errcheck = _check_return_code

    # Function "chfl_topology_add_residue", at topology.h:227
    c_lib.chfl_topology_add_residue.argtypes = [POINTER(CHFL_TOPOLOGY), POINTER(CHFL_RESIDUE)]
    c_lib.chfl_topology_add_residue.restype = chfl_status
    c_lib.chfl_topology_add_residue.errcheck = _check_return_code

    # Function "chfl_topology_residues_linked", at topology.h:238
    c_lib.chfl_topology_residues_linked.argtypes = [POINTER(CHFL_TOPOLOGY), POINTER(CHFL_RESIDUE), POINTER(CHFL_RESIDUE), POINTER(c_bool)]
    c_lib.chfl_topology_residues_linked.restype = chfl_status
    c_lib.chfl_topology_residues_linked.errcheck = _check_return_code

    # Function "chfl_topology_bond_with_order", at topology.h:251
    c_lib.chfl_topology_bond_with_order.argtypes = [POINTER(CHFL_TOPOLOGY), c_uint64, c_uint64, chfl_bond_order]
    c_lib.chfl_topology_bond_with_order.restype = chfl_status
    c_lib.chfl_topology_bond_with_order.errcheck = _check_return_code

    # Function "chfl_topology_bond_orders", at topology.h:265
    c_lib.chfl_topology_bond_orders.argtypes = [POINTER(CHFL_TOPOLOGY), ndpointer(chfl_bond_order, flags="C_CONTIGUOUS", ndim=1), c_uint64]
    c_lib.chfl_topology_bond_orders.restype = chfl_status
    c_lib.chfl_topology_bond_orders.errcheck = _check_return_code

    # Function "chfl_topology_bond_order", at topology.h:278
    c_lib.chfl_topology_bond_order.argtypes = [POINTER(CHFL_TOPOLOGY), c_uint64, c_uint64, POINTER(chfl_bond_order)]
    c_lib.chfl_topology_bond_order.restype = chfl_status
    c_lib.chfl_topology_bond_order.errcheck = _check_return_code

    # Function "chfl_cell", at cell.h:40
    c_lib.chfl_cell.argtypes = [chfl_vector3d, chfl_vector3d]
    c_lib.chfl_cell.restype = POINTER(CHFL_CELL)

    # Function "chfl_cell_from_matrix", at cell.h:55
    c_lib.chfl_cell_from_matrix.argtypes = [ARRAY(chfl_vector3d, (3))]
    c_lib.chfl_cell_from_matrix.restype = POINTER(CHFL_CELL)

    # Function "chfl_cell_from_frame", at cell.h:68
    c_lib.chfl_cell_from_frame.argtypes = [POINTER(CHFL_FRAME)]
    c_lib.chfl_cell_from_frame.restype = POINTER(CHFL_CELL)

    # Function "chfl_cell_copy", at cell.h:78
    c_lib.chfl_cell_copy.argtypes = [POINTER(CHFL_CELL)]
    c_lib.chfl_cell_copy.restype = POINTER(CHFL_CELL)

    # Function "chfl_cell_volume", at cell.h:85
    c_lib.chfl_cell_volume.argtypes = [POINTER(CHFL_CELL), POINTER(c_double)]
    c_lib.chfl_cell_volume.restype = chfl_status
    c_lib.chfl_cell_volume.errcheck = _check_return_code

    # Function "chfl_cell_lengths", at cell.h:94
    c_lib.chfl_cell_lengths.argtypes = [POINTER(CHFL_CELL), chfl_vector3d]
    c_lib.chfl_cell_lengths.restype = chfl_status
    c_lib.chfl_cell_lengths.errcheck = _check_return_code

    # Function "chfl_cell_set_lengths", at cell.h:110
    c_lib.chfl_cell_set_lengths.argtypes = [POINTER(CHFL_CELL), chfl_vector3d]
    c_lib.chfl_cell_set_lengths.restype = chfl_status
    c_lib.chfl_cell_set_lengths.errcheck = _check_return_code

    # Function "chfl_cell_angles", at cell.h:119
    c_lib.chfl_cell_angles.argtypes = [POINTER(CHFL_CELL), chfl_vector3d]
    c_lib.chfl_cell_angles.restype = chfl_status
    c_lib.chfl_cell_angles.errcheck = _check_return_code

    # Function "chfl_cell_set_angles", at cell.h:137
    c_lib.chfl_cell_set_angles.argtypes = [POINTER(CHFL_CELL), chfl_vector3d]
    c_lib.chfl_cell_set_angles.restype = chfl_status
    c_lib.chfl_cell_set_angles.errcheck = _check_return_code

    # Function "chfl_cell_matrix", at cell.h:146
    c_lib.chfl_cell_matrix.argtypes = [POINTER(CHFL_CELL), ARRAY(chfl_vector3d, (3))]
    c_lib.chfl_cell_matrix.restype = chfl_status
    c_lib.chfl_cell_matrix.errcheck = _check_return_code

    # Function "chfl_cell_shape", at cell.h:155
    c_lib.chfl_cell_shape.argtypes = [POINTER(CHFL_CELL), POINTER(chfl_cellshape)]
    c_lib.chfl_cell_shape.restype = chfl_status
    c_lib.chfl_cell_shape.errcheck = _check_return_code

    # Function "chfl_cell_set_shape", at cell.h:164
    c_lib.chfl_cell_set_shape.argtypes = [POINTER(CHFL_CELL), chfl_cellshape]
    c_lib.chfl_cell_set_shape.restype = chfl_status
    c_lib.chfl_cell_set_shape.errcheck = _check_return_code

    # Function "chfl_cell_wrap", at cell.h:173
    c_lib.chfl_cell_wrap.argtypes = [POINTER(CHFL_CELL), chfl_vector3d]
    c_lib.chfl_cell_wrap.restype = chfl_status
    c_lib.chfl_cell_wrap.errcheck = _check_return_code

    # Function "chfl_frame", at frame.h:25
    c_lib.chfl_frame.argtypes = []
    c_lib.chfl_frame.restype = POINTER(CHFL_FRAME)

    # Function "chfl_frame_copy", at frame.h:35
    c_lib.chfl_frame_copy.argtypes = [POINTER(CHFL_FRAME)]
    c_lib.chfl_frame_copy.restype = POINTER(CHFL_FRAME)

    # Function "chfl_frame_atoms_count", at frame.h:43
    c_lib.chfl_frame_atoms_count.argtypes = [POINTER(CHFL_FRAME), POINTER(c_uint64)]
    c_lib.chfl_frame_atoms_count.restype = chfl_status
    c_lib.chfl_frame_atoms_count.errcheck = _check_return_code

    # Function "chfl_frame_positions", at frame.h:62
    c_lib.chfl_frame_positions.argtypes = [POINTER(CHFL_FRAME), POINTER(POINTER(chfl_vector3d)), POINTER(c_uint64)]
    c_lib.chfl_frame_positions.restype = chfl_status
    c_lib.chfl_frame_positions.errcheck = _check_return_code

    # Function "chfl_frame_velocities", at frame.h:85
    c_lib.chfl_frame_velocities.argtypes = [POINTER(CHFL_FRAME), POINTER(POINTER(chfl_vector3d)), POINTER(c_uint64)]
    c_lib.chfl_frame_velocities.restype = chfl_status
    c_lib.chfl_frame_velocities.errcheck = _check_return_code

    # Function "chfl_frame_add_atom", at frame.h:97
    c_lib.chfl_frame_add_atom.argtypes = [POINTER(CHFL_FRAME), POINTER(CHFL_ATOM), chfl_vector3d, chfl_vector3d]
    c_lib.chfl_frame_add_atom.restype = chfl_status
    c_lib.chfl_frame_add_atom.errcheck = _check_return_code

    # Function "chfl_frame_remove", at frame.h:110
    c_lib.chfl_frame_remove.argtypes = [POINTER(CHFL_FRAME), c_uint64]
    c_lib.chfl_frame_remove.restype = chfl_status
    c_lib.chfl_frame_remove.errcheck = _check_return_code

    # Function "chfl_frame_resize", at frame.h:122
    c_lib.chfl_frame_resize.argtypes = [POINTER(CHFL_FRAME), c_uint64]
    c_lib.chfl_frame_resize.restype = chfl_status
    c_lib.chfl_frame_resize.errcheck = _check_return_code

    # Function "chfl_frame_add_velocities", at frame.h:134
    c_lib.chfl_frame_add_velocities.argtypes = [POINTER(CHFL_FRAME)]
    c_lib.chfl_frame_add_velocities.restype = chfl_status
    c_lib.chfl_frame_add_velocities.errcheck = _check_return_code

    # Function "chfl_frame_has_velocities", at frame.h:142
    c_lib.chfl_frame_has_velocities.argtypes = [POINTER(CHFL_FRAME), POINTER(c_bool)]
    c_lib.chfl_frame_has_velocities.restype = chfl_status
    c_lib.chfl_frame_has_velocities.errcheck = _check_return_code

    # Function "chfl_frame_set_cell", at frame.h:151
    c_lib.chfl_frame_set_cell.argtypes = [POINTER(CHFL_FRAME), POINTER(CHFL_CELL)]
    c_lib.chfl_frame_set_cell.restype = chfl_status
    c_lib.chfl_frame_set_cell.errcheck = _check_return_code

    # Function "chfl_frame_set_topology", at frame.h:163
    c_lib.chfl_frame_set_topology.argtypes = [POINTER(CHFL_FRAME), POINTER(CHFL_TOPOLOGY)]
    c_lib.chfl_frame_set_topology.restype = chfl_status
    c_lib.chfl_frame_set_topology.errcheck = _check_return_code

    # Function "chfl_frame_step", at frame.h:173
    c_lib.chfl_frame_step.argtypes = [POINTER(CHFL_FRAME), POINTER(c_uint64)]
    c_lib.chfl_frame_step.restype = chfl_status
    c_lib.chfl_frame_step.errcheck = _check_return_code

    # Function "chfl_frame_set_step", at frame.h:182
    c_lib.chfl_frame_set_step.argtypes = [POINTER(CHFL_FRAME), c_uint64]
    c_lib.chfl_frame_set_step.restype = chfl_status
    c_lib.chfl_frame_set_step.errcheck = _check_return_code

    # Function "chfl_frame_guess_bonds", at frame.h:194
    c_lib.chfl_frame_guess_bonds.argtypes = [POINTER(CHFL_FRAME)]
    c_lib.chfl_frame_guess_bonds.restype = chfl_status
    c_lib.chfl_frame_guess_bonds.errcheck = _check_return_code

    # Function "chfl_frame_distance", at frame.h:203
    c_lib.chfl_frame_distance.argtypes = [POINTER(CHFL_FRAME), c_uint64, c_uint64, POINTER(c_double)]
    c_lib.chfl_frame_distance.restype = chfl_status
    c_lib.chfl_frame_distance.errcheck = _check_return_code

    # Function "chfl_frame_angle", at frame.h:214
    c_lib.chfl_frame_angle.argtypes = [POINTER(CHFL_FRAME), c_uint64, c_uint64, c_uint64, POINTER(c_double)]
    c_lib.chfl_frame_angle.restype = chfl_status
    c_lib.chfl_frame_angle.errcheck = _check_return_code

    # Function "chfl_frame_dihedral", at frame.h:225
    c_lib.chfl_frame_dihedral.argtypes = [POINTER(CHFL_FRAME), c_uint64, c_uint64, c_uint64, c_uint64, POINTER(c_double)]
    c_lib.chfl_frame_dihedral.restype = chfl_status
    c_lib.chfl_frame_dihedral.errcheck = _check_return_code

    # Function "chfl_frame_out_of_plane", at frame.h:239
    c_lib.chfl_frame_out_of_plane.argtypes = [POINTER(CHFL_FRAME), c_uint64, c_uint64, c_uint64, c_uint64, POINTER(c_double)]
    c_lib.chfl_frame_out_of_plane.restype = chfl_status
    c_lib.chfl_frame_out_of_plane.errcheck = _check_return_code

    # Function "chfl_frame_properties_count", at frame.h:248
    c_lib.chfl_frame_properties_count.argtypes = [POINTER(CHFL_FRAME), POINTER(c_uint64)]
    c_lib.chfl_frame_properties_count.restype = chfl_status
    c_lib.chfl_frame_properties_count.errcheck = _check_return_code

    # Function "chfl_frame_list_properties", at frame.h:264
    c_lib.chfl_frame_list_properties.argtypes = [POINTER(CHFL_FRAME), POINTER(c_char_p), c_uint64]
    c_lib.chfl_frame_list_properties.restype = chfl_status
    c_lib.chfl_frame_list_properties.errcheck = _check_return_code

    # Function "chfl_frame_set_property", at frame.h:276
    c_lib.chfl_frame_set_property.argtypes = [POINTER(CHFL_FRAME), c_char_p, POINTER(CHFL_PROPERTY)]
    c_lib.chfl_frame_set_property.restype = chfl_status
    c_lib.chfl_frame_set_property.errcheck = _check_return_code

    # Function "chfl_frame_get_property", at frame.h:290
    c_lib.chfl_frame_get_property.argtypes = [POINTER(CHFL_FRAME), c_char_p]
    c_lib.chfl_frame_get_property.restype = POINTER(CHFL_PROPERTY)

    # Function "chfl_frame_add_bond", at frame.h:299
    c_lib.chfl_frame_add_bond.argtypes = [POINTER(CHFL_FRAME), c_uint64, c_uint64]
    c_lib.chfl_frame_add_bond.restype = chfl_status
    c_lib.chfl_frame_add_bond.errcheck = _check_return_code

    # Function "chfl_frame_bond_with_order", at frame.h:309
    c_lib.chfl_frame_bond_with_order.argtypes = [POINTER(CHFL_FRAME), c_uint64, c_uint64, chfl_bond_order]
    c_lib.chfl_frame_bond_with_order.restype = chfl_status
    c_lib.chfl_frame_bond_with_order.errcheck = _check_return_code

    # Function "chfl_frame_remove_bond", at frame.h:321
    c_lib.chfl_frame_remove_bond.argtypes = [POINTER(CHFL_FRAME), c_uint64, c_uint64]
    c_lib.chfl_frame_remove_bond.restype = chfl_status
    c_lib.chfl_frame_remove_bond.errcheck = _check_return_code

    # Function "chfl_frame_clear_bonds", at frame.h:331
    c_lib.chfl_frame_clear_bonds.argtypes = [POINTER(CHFL_FRAME)]
    c_lib.chfl_frame_clear_bonds.restype = chfl_status
    c_lib.chfl_frame_clear_bonds.errcheck = _check_return_code

    # Function "chfl_frame_add_residue", at frame.h:341
    c_lib.chfl_frame_add_residue.argtypes = [POINTER(CHFL_FRAME), POINTER(CHFL_RESIDUE)]
    c_lib.chfl_frame_add_residue.restype = chfl_status
    c_lib.chfl_frame_add_residue.errcheck = _check_return_code

    # Function "chfl_trajectory_open", at trajectory.h:26
    c_lib.chfl_trajectory_open.argtypes = [c_char_p, c_char]
    c_lib.chfl_trajectory_open.restype = POINTER(CHFL_TRAJECTORY)

    # Function "chfl_trajectory_with_format", at trajectory.h:43
    c_lib.chfl_trajectory_with_format.argtypes = [c_char_p, c_char, c_char_p]
    c_lib.chfl_trajectory_with_format.restype = POINTER(CHFL_TRAJECTORY)

    # Function "chfl_trajectory_memory_reader", at trajectory.h:59
    c_lib.chfl_trajectory_memory_reader.argtypes = [c_char_p, c_uint64, c_char_p]
    c_lib.chfl_trajectory_memory_reader.restype = POINTER(CHFL_TRAJECTORY)

    # Function "chfl_trajectory_memory_writer", at trajectory.h:74
    c_lib.chfl_trajectory_memory_writer.argtypes = [c_char_p]
    c_lib.chfl_trajectory_memory_writer.restype = POINTER(CHFL_TRAJECTORY)

    # Function "chfl_trajectory_path", at trajectory.h:86
    c_lib.chfl_trajectory_path.argtypes = [POINTER(CHFL_TRAJECTORY), c_char_p, c_uint64]
    c_lib.chfl_trajectory_path.restype = chfl_status
    c_lib.chfl_trajectory_path.errcheck = _check_return_code

    # Function "chfl_trajectory_read", at trajectory.h:98
    c_lib.chfl_trajectory_read.argtypes = [POINTER(CHFL_TRAJECTORY), POINTER(CHFL_FRAME)]
    c_lib.chfl_trajectory_read.restype = chfl_status
    c_lib.chfl_trajectory_read.errcheck = _check_return_code

    # Function "chfl_trajectory_read_step", at trajectory.h:110
    c_lib.chfl_trajectory_read_step.argtypes = [POINTER(CHFL_TRAJECTORY), c_uint64, POINTER(CHFL_FRAME)]
    c_lib.chfl_trajectory_read_step.restype = chfl_status
    c_lib.chfl_trajectory_read_step.errcheck = _check_return_code

    # Function "chfl_trajectory_write", at trajectory.h:119
    c_lib.chfl_trajectory_write.argtypes = [POINTER(CHFL_TRAJECTORY), POINTER(CHFL_FRAME)]
    c_lib.chfl_trajectory_write.restype = chfl_status
    c_lib.chfl_trajectory_write.errcheck = _check_return_code

    # Function "chfl_trajectory_set_topology", at trajectory.h:130
    c_lib.chfl_trajectory_set_topology.argtypes = [POINTER(CHFL_TRAJECTORY), POINTER(CHFL_TOPOLOGY)]
    c_lib.chfl_trajectory_set_topology.restype = chfl_status
    c_lib.chfl_trajectory_set_topology.errcheck = _check_return_code

    # Function "chfl_trajectory_topology_file", at trajectory.h:144
    c_lib.chfl_trajectory_topology_file.argtypes = [POINTER(CHFL_TRAJECTORY), c_char_p, c_char_p]
    c_lib.chfl_trajectory_topology_file.restype = chfl_status
    c_lib.chfl_trajectory_topology_file.errcheck = _check_return_code

    # Function "chfl_trajectory_set_cell", at trajectory.h:154
    c_lib.chfl_trajectory_set_cell.argtypes = [POINTER(CHFL_TRAJECTORY), POINTER(CHFL_CELL)]
    c_lib.chfl_trajectory_set_cell.restype = chfl_status
    c_lib.chfl_trajectory_set_cell.errcheck = _check_return_code

    # Function "chfl_trajectory_nsteps", at trajectory.h:164
    c_lib.chfl_trajectory_nsteps.argtypes = [POINTER(CHFL_TRAJECTORY), POINTER(c_uint64)]
    c_lib.chfl_trajectory_nsteps.restype = chfl_status
    c_lib.chfl_trajectory_nsteps.errcheck = _check_return_code

    # Function "chfl_trajectory_memory_buffer", at trajectory.h:178
    c_lib.chfl_trajectory_memory_buffer.argtypes = [POINTER(CHFL_TRAJECTORY), POINTER(c_char_p), POINTER(c_uint64)]
    c_lib.chfl_trajectory_memory_buffer.restype = chfl_status
    c_lib.chfl_trajectory_memory_buffer.errcheck = _check_return_code

    # Function "chfl_selection", at selection.h:24
    c_lib.chfl_selection.argtypes = [c_char_p]
    c_lib.chfl_selection.restype = POINTER(CHFL_SELECTION)

    # Function "chfl_selection_copy", at selection.h:37
    c_lib.chfl_selection_copy.argtypes = [POINTER(CHFL_SELECTION)]
    c_lib.chfl_selection_copy.restype = POINTER(CHFL_SELECTION)

    # Function "chfl_selection_size", at selection.h:49
    c_lib.chfl_selection_size.argtypes = [POINTER(CHFL_SELECTION), POINTER(c_uint64)]
    c_lib.chfl_selection_size.restype = chfl_status
    c_lib.chfl_selection_size.errcheck = _check_return_code

    # Function "chfl_selection_string", at selection.h:62
    c_lib.chfl_selection_string.argtypes = [POINTER(CHFL_SELECTION), c_char_p, c_uint64]
    c_lib.chfl_selection_string.restype = chfl_status
    c_lib.chfl_selection_string.errcheck = _check_return_code

    # Function "chfl_selection_evaluate", at selection.h:75
    c_lib.chfl_selection_evaluate.argtypes = [POINTER(CHFL_SELECTION), POINTER(CHFL_FRAME), POINTER(c_uint64)]
    c_lib.chfl_selection_evaluate.restype = chfl_status
    c_lib.chfl_selection_evaluate.errcheck = _check_return_code

    # Function "chfl_selection_matches", at selection.h:87
    c_lib.chfl_selection_matches.argtypes = [POINTER(CHFL_SELECTION), ndpointer(chfl_match, flags="C_CONTIGUOUS", ndim=1), c_uint64]
    c_lib.chfl_selection_matches.restype = chfl_status
    c_lib.chfl_selection_matches.errcheck = _check_return_code
