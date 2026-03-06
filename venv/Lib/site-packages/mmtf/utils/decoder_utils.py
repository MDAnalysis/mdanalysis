import sys

from mmtf.utils import constants
def add_atom_data(data_api, data_setters, atom_names, element_names, atom_charges, group_atom_ind):
    """Add the atomic data to the DataTransferInterface.
    :param data_api the data api from where to get the data
    :param data_setters the class to push the data to
    :param atom_nams the list of atom names for the group
    :param element_names the list of element names for this group
    :param atom_charges the list formal atomic charges for this group
    :param group_atom_ind the index of this atom in the group"""
    atom_name = atom_names[group_atom_ind]
    element = element_names[group_atom_ind]
    charge = atom_charges[group_atom_ind]
    alternative_location_id = data_api.alt_loc_list[data_api.atom_counter]
    serial_number = data_api.atom_id_list[data_api.atom_counter]
    x = data_api.x_coord_list[data_api.atom_counter]
    y = data_api.y_coord_list[data_api.atom_counter]
    z = data_api.z_coord_list[data_api.atom_counter]
    occupancy = data_api.occupancy_list[data_api.atom_counter]
    temperature_factor = data_api.b_factor_list[data_api.atom_counter]
    data_setters.set_atom_info(atom_name, serial_number, alternative_location_id,
                               x, y, z, occupancy, temperature_factor, element, charge)


def add_group_bonds(data_setters, bond_indices, bond_orders):
    """Add the bonds for this group.
    :param data_setters the class to push the data to
    :param bond_indices the indices of the atoms in the group that
    are bonded (in pairs)
    :param bond_orders the orders of the bonds"""
    for bond_index in range(len(bond_orders)):
        data_setters.set_group_bond(bond_indices[bond_index*2],bond_indices[bond_index*2+1],bond_orders[bond_index])


def add_group(data_api, data_setters, group_index):
    """Add the data for a whole group.
    :param data_api the data api from where to get the data
    :param data_setters the class to push the data to
    :param group_index the index for this group"""
    group_type_ind = data_api.group_type_list[group_index]
    atom_count = len(data_api.group_list[group_type_ind]["atomNameList"])
    insertion_code = data_api.ins_code_list[group_index]
    data_setters.set_group_info(data_api.group_list[group_type_ind]["groupName"],
                                data_api.group_id_list[group_index], insertion_code,
                                data_api.group_list[group_type_ind]["chemCompType"],
                                atom_count, data_api.num_bonds,
                                data_api.group_list[group_type_ind]["singleLetterCode"],
                                data_api.sequence_index_list[group_index],
                                data_api.sec_struct_list[group_index])
    for group_atom_ind in range(atom_count):
        add_atom_data(data_api, data_setters,
                      data_api.group_list[group_type_ind]["atomNameList"],
                      data_api.group_list[group_type_ind]["elementList"],
                      data_api.group_list[group_type_ind]["formalChargeList"],
                      group_atom_ind)
        data_api.atom_counter +=1
    add_group_bonds(data_setters,
                    data_api.group_list[group_type_ind]["bondAtomList"],
                    data_api.group_list[group_type_ind]["bondOrderList"])
    return atom_count


def add_chain_info(data_api, data_setters, chain_index):
    """Add the data for a whole chain.
    :param data_api the data api from where to get the data
    :param data_setters the class to push the data to
    :param chain_index the index for this chain"""
    chain_id = data_api.chain_id_list[chain_index]
    chain_name = data_api.chain_name_list[chain_index]
    num_groups = data_api.groups_per_chain[chain_index]
    data_setters.set_chain_info(chain_id, chain_name, num_groups)
    next_ind = data_api.group_counter + num_groups
    last_ind = data_api.group_counter
    for group_ind in range(last_ind, next_ind):
        add_group(data_api, data_setters, group_ind)
        data_api.group_counter +=1
    data_api.chain_counter+=1


def add_atomic_information(data_api, data_setters):
    """Add all the structural information.
    :param data_api the data api from where to get the data
    :param data_setters the class to push the data to"""
    for model_chains in data_api.chains_per_model:
        data_setters.set_model_info(data_api.model_counter, model_chains)
        tot_chains_this_model = data_api.chain_counter + model_chains
        last_chain_counter = data_api.chain_counter
        for chain_index in range(last_chain_counter, tot_chains_this_model):
            add_chain_info(data_api, data_setters, chain_index)
        data_api.model_counter+=1


def generate_bio_assembly(data_api, struct_inflator):
    """Generate the bioassembly data.
    :param data_api the interface to the decoded data
    :param struct_inflator the interface to put the data into the client object"""
    bioassembly_count = 0
    for bioassembly in data_api.bio_assembly:
        bioassembly_count += 1
        for transform in bioassembly["transformList"]:
            struct_inflator.set_bio_assembly_trans(bioassembly_count,
                                                   transform["chainIndexList"],
                                                   transform["matrix"])

def add_inter_group_bonds(data_api, struct_inflator):
    """	 Generate inter group bonds.
	 Bond indices are specified within the whole structure and start at 0.
	 :param data_api the interface to the decoded data
	 :param struct_inflator the interface to put the data into the client object"""
    for i in range(len(data_api.bond_order_list)):
        struct_inflator.set_inter_group_bond(data_api.bond_atom_list[i * 2],
                                             data_api.bond_atom_list[i * 2 + 1],
                                             data_api.bond_order_list[i])



def add_header_info(data_api, struct_inflator):
    """ Add ancilliary header information to the structure.
	 :param data_api the interface to the decoded data
	 :param struct_inflator the interface to put the data into the client object
	 """
    struct_inflator.set_header_info(data_api.r_free,
                                    data_api.r_work,
                                    data_api.resolution,
                                    data_api.title,
                                    data_api.deposition_date,
                                    data_api.release_date,
                                    data_api.experimental_methods)



def add_xtalographic_info(data_api, struct_inflator):
    """	 Add the crystallographic data to the structure.
	 :param data_api the interface to the decoded data
	 :param struct_inflator the interface to put the data into the client object"""
    if data_api.unit_cell == None and data_api.space_group is not None:
        struct_inflator.set_xtal_info(data_api.space_group,
                                      constants.UNKNOWN_UNIT_CELL)
    elif data_api.unit_cell is not None and data_api.space_group is None:
        struct_inflator.set_xtal_info(constants.UNKNOWN_SPACE_GROUP,
                                      data_api.unit_cell)
    elif data_api.unit_cell is None and data_api.space_group is None:
        struct_inflator.set_xtal_info(constants.UNKNOWN_SPACE_GROUP,
                                      constants.UNKNOWN_UNIT_CELL)
    else:
        struct_inflator.set_xtal_info(data_api.space_group,
                                      data_api.unit_cell)

def add_entity_info( data_api, struct_inflator):
    """Add the entity info to the structure.
    :param data_api the interface to the decoded data
    :param struct_inflator the interface to put the data into the client object
    """
    for entity in data_api.entity_list:
        struct_inflator.set_entity_info(entity["chainIndexList"],
                                        entity["sequence"],
                                        entity.get("description", ""),
                                        entity["type"])


def get_bonds(input_group):
    """Utility function to get indices (in pairs) of the bonds."""
    out_list = []
    for i in range(len(input_group.bond_order_list)):
        out_list.append((input_group.bond_atom_list[i * 2], input_group.bond_atom_list[i * 2 + 1],))
    return out_list
