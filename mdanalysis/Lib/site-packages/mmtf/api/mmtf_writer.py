from mmtf.codecs import encode_array
import msgpack
from mmtf.utils import constants

def make_entity_dict(chain_indices,sequence,description,entity_type):
    out_d = {}
    out_d["description"] = description
    out_d["type"] = entity_type
    out_d["chainIndexList"] = chain_indices
    out_d["sequence"] = sequence
    return out_d

class Group(object):

    def __eq__(self, other):
        """Function to define equality"""
        if self.atom_name_list != other.atom_name_list:
            return False
        if self.charge_list != other.charge_list:
            return False
        if self.element_list != other.element_list:
            return False
        if self.group_type != other.group_type:
            return False
        if self.group_name != other.group_name:
            return False
        if self.single_letter_code != other.single_letter_code:
            return False
        if self.bond_atom_list != other.bond_atom_list:
            return False
        if self.bond_order_list != other.bond_order_list:
            return False
        return True

    def __init__(self):
        self.atom_name_list = []
        self.bond_order_list = []
        self.bond_atom_list = []
        self.charge_list = []
        self.element_list = []
        self.group_name = constants.UNKOWN_GROUP_NAME
        self.group_type = constants.UNKOWN_GROUP_TYPE
        self.single_letter_code = constants.UNKNOWN_SL

    def convert_to_dict(self):
        """Convert the group object to an appropriate DICT"""
        out_dict = {}
        out_dict["groupName"] = self.group_name
        out_dict["atomNameList"] = self.atom_name_list
        out_dict["elementList"] = self.element_list
        out_dict["bondOrderList"] = self.bond_order_list
        out_dict["bondAtomList"] = self.bond_atom_list
        out_dict["formalChargeList"] = self.charge_list
        out_dict["singleLetterCode"] = self.single_letter_code
        out_dict["chemCompType"] = self.group_type
        return out_dict


def get_unique_groups(input_list):
    """Function to get a unique list of groups."""
    out_list = []
    for item in input_list:
        if item not in out_list:
            out_list.append(item)
    return out_list


class TemplateEncoder(object):
    """Template class to be used by third parties to pass data into other data structures."""

    def init_structure(self, total_num_bonds, total_num_atoms,
                       total_num_groups, total_num_chains, total_num_models,
                       structure_id):
        """Initialise the structure object.
        :param total_num_bonds: the number of bonds in the structure
        :param total_num_atoms: the number of atoms in the structure
        :param total_num_groups: the number of groups in the structure
        :param total_num_chains: the number of chains in the structure
        :param total_num_models: the number of models in the structure
        :param structure_id the: id of the structure (e.g. PDB id)
        """
        raise NotImplementedError

    def set_atom_info(self, atom_name, serial_number, alternative_location_id,
                      x, y, z, occupancy, temperature_factor, element, charge):
        """Create an atom object an set the information.
        :param atom_name: the atom name, e.g. CA for this atom
        :param serial_number: the serial id of the atom (e.g. 1)
        :param alternative_location_id: the alternative location id for the atom, if present
        :param x: the x coordiante of the atom
        :param y: the y coordinate of the atom
        :param z: the z coordinate of the atom
        :param occupancy: the occupancy of the atom
        :param temperature_factor: the temperature factor of the atom
        :param element: the element of the atom, e.g. C for carbon. According to IUPAC. Calcium  is Ca
        :param charge: the formal atomic charge of the atom
        """
        raise NotImplementedError


    def set_chain_info(self, chain_id, chain_name, num_groups):
        """Set the chain information.
        :param chain_id: the asym chain id from mmCIF
        :param chain_name: the auth chain id from mmCIF
        :param num_groups: the number of groups this chain has
        """
        raise NotImplementedError



    def set_entity_info(self, chain_indices, sequence, description, entity_type):
        """Set the entity level information for the structure.
        :param chain_indices: the indices of the chains for this entity
        :param sequence: the one letter code sequence for this entity
        :param description: the description for this entity
        :param entity_type: the entity type (polymer,non-polymer,water)
        """
        raise NotImplementedError


    def set_group_info(self, group_name, group_number, insertion_code,
                       group_type, atom_count, bond_count, single_letter_code,
                       sequence_index, secondary_structure_type):
        """Set the information for a group
        :param group_name: the name of this group,e.g. LYS
        :param group_number: the residue number of this group
        :param insertion_code: the insertion code for this group
        :param group_type: a string indicating the type of group (as found in the chemcomp dictionary.
        Empty string if none available.
        :param atom_count: the number of atoms in the group
        :param bond_count: the number of unique bonds in the group
        :param single_letter_code: the single letter code of the group
        :param sequence_index: the index of this group in the sequence defined by the enttiy
        :param secondary_structure_type: the type of secondary structure used (types are according to DSSP and
        number to type mappings are defined in the specification)
        """
        raise NotImplementedError



    def set_model_info(self, model_id, chain_count):
        # FIXME model_id here is meaningles and potentially misleading.
        """Set the information for a model.
        :param model_id: the index for the model
        :param chain_count: the number of chains in the model
        """
        raise NotImplementedError


    def set_xtal_info(self, space_group, unit_cell):
        """Set the crystallographic information for the structure
        :param space_group: the space group name, e.g. "P 21 21 21"
        :param unit_cell: an array of length 6 with the unit cell parameters in order: a, b, c, alpha, beta, gamma
        """
        raise NotImplementedError


    def set_header_info(self, r_free, r_work, resolution, title,
                        deposition_date, release_date, experimental_methods):
        """Sets the header information.
        :param r_free: the measured R-Free for the structure
        :param r_work: the measure R-Work for the structure
        :param resolution: the resolution of the structure
        :param title: the title of the structure
        :param deposition_date: the deposition date of the structure
        :param release_date: the release date of the structure
        :param experimnetal_methods: the list of experimental methods in the structure
        """
        raise NotImplementedError


    def set_bio_assembly_trans(self, bio_assembly_index, input_chain_indices, input_transform):
        """Set the Bioassembly transformation information. A single bioassembly can have multiple transforms,
        :param bio_assembly_index: the integer index of the bioassembly
        :param input_chain_indices: the list of integer indices for the chains of this bioassembly
        :param input_transformation: the list of doubles for  the transform of this bioassmbly transform"""
        raise NotImplementedError



    def finalize_structure(self):
        """Any functions needed to cleanup the structure."""
        raise NotImplementedError



    def set_group_bond(self, atom_index_one, atom_index_two, bond_order):
        """Add bonds within a group.
        :param atom_index_one: the integer atom index (in the group) of the first partner in the bond
        :param atom_index_two: the integer atom index (in the group) of the second partner in the bond
        :param bond_order: the integer bond order
        """
        raise NotImplementedError



    def set_inter_group_bond(self, atom_index_one, atom_index_two, bond_order):
        """Add bonds between groups.
        :param atom_index_one: the integer atom index (in the structure) of the first partner in the bond
        :param atom_index_two: the integer atom index (in the structure) of the second partner in the bond
        :param bond_order the bond order
        """
        raise NotImplementedError



class MMTFEncoder(TemplateEncoder):

    def encode_data(self):
        """Encode the data back into a dict."""
        output_data = {}
        output_data["groupTypeList"] = encode_array(self.group_type_list, 4, 0)
        output_data["xCoordList"] = encode_array(self.x_coord_list, 10, 1000)
        output_data["yCoordList"] = encode_array(self.y_coord_list, 10, 1000)
        output_data["zCoordList"] = encode_array(self.z_coord_list, 10, 1000)
        output_data["bFactorList"] = encode_array(self.b_factor_list, 10, 100)
        output_data["occupancyList"] = encode_array(self.occupancy_list, 9, 100)
        output_data["atomIdList"] = encode_array(self.atom_id_list, 8, 0)
        output_data["altLocList"] = encode_array(self.alt_loc_list, 6, 0)
        output_data["insCodeList"] = encode_array(self.ins_code_list, 6, 0)
        output_data["groupIdList"] = encode_array(self.group_id_list, 8, 0)
        output_data["groupList"] = self.group_list
        output_data["sequenceIndexList"] = encode_array(self.sequence_index_list, 8, 0)
        output_data["chainNameList"] = encode_array(self.chain_name_list, 5, 4)
        output_data["chainIdList"] = encode_array(self.chain_id_list, 5, 4)
        output_data["bondAtomList"] = encode_array(self.bond_atom_list, 4, 0)
        output_data["bondOrderList"] = encode_array(self.bond_order_list, 2, 0)
        output_data["secStructList"] = encode_array(self.sec_struct_list, 2, 0)
        output_data["chainsPerModel"] = self.chains_per_model
        output_data["groupsPerChain"] = self.groups_per_chain
        output_data["spaceGroup"] = self.space_group
        output_data["mmtfVersion"] = self.mmtf_version
        output_data["mmtfProducer"] = self.mmtf_producer
        output_data["structureId"] = self.structure_id
        output_data["entityList"] = self.entity_list
        output_data["bioAssemblyList"] = self.bio_assembly
        output_data["rFree"] = self.r_free
        output_data["rWork"] = self.r_work
        output_data["resolution"] = self.resolution
        output_data["title"] = self.title
        output_data["experimentalMethods"] = self.experimental_methods
        output_data["depositionDate"] = self.deposition_date
        output_data["releaseDate"] = self.release_date
        output_data["unitCell"] = self.unit_cell
        output_data["numBonds"] = self.num_bonds
        output_data["numChains"] = self.num_chains
        output_data["numModels"] = self.num_models
        output_data["numAtoms"] = self.num_atoms
        output_data["numGroups"] = self.num_groups
        return output_data


    def get_msgpack(self):
        """Get the msgpack of the encoded data."""
        return msgpack.packb(self.encode_data(), use_bin_type=True)


    def write_file(self, file_path):
        with open(file_path, "wb") as out_f:
            out_f.write(self.get_msgpack())


    def init_structure(self, total_num_bonds, total_num_atoms,
                       total_num_groups, total_num_chains, total_num_models,
                       structure_id):
        """Initialise the structure object.
        :param total_num_bonds: the number of bonds in the structure
        :param total_num_atoms: the number of atoms in the structure
        :param total_num_groups: the number of groups in the structure
        :param total_num_chains: the number of chains in the structure
        :param total_num_models: the number of models in the structure
        :param structure_id the: id of the structure (e.g. PDB id)
        """
        self.mmtf_version = constants.MMTF_VERSION
        self.mmtf_producer = constants.PRODUCER
        self.num_atoms = total_num_atoms
        self.num_bonds = total_num_bonds
        self.num_groups = total_num_groups
        self.num_chains = total_num_chains
        self.num_models = total_num_models
        self.structure_id = structure_id
        # initialise the arrays
        self.x_coord_list = []
        self.y_coord_list = []
        self.z_coord_list = []
        self.group_type_list = []
        self.entity_list = []
        self.b_factor_list = []
        self.occupancy_list = []
        self.atom_id_list = []
        self.alt_loc_list = []
        self.ins_code_list = []
        self.group_id_list = []
        self.sequence_index_list = []
        self.group_list = []
        self.chain_name_list = []
        self.chain_id_list = []
        self.bond_atom_list = []
        self.bond_order_list = []
        self.sec_struct_list = []
        self.chains_per_model = []
        self.groups_per_chain = []
        self.current_group = None
        self.bio_assembly = []


    def set_atom_info(self, atom_name, serial_number, alternative_location_id,
                      x, y, z, occupancy, temperature_factor, element, charge):
        """Create an atom object an set the information.
        :param atom_name: the atom name, e.g. CA for this atom
        :param serial_number: the serial id of the atom (e.g. 1)
        :param alternative_location_id: the alternative location id for the atom, if present
        :param x: the x coordiante of the atom
        :param y: the y coordinate of the atom
        :param z: the z coordinate of the atom
        :param occupancy: the occupancy of the atom
        :param temperature_factor: the temperature factor of the atom
        :param element: the element of the atom, e.g. C for carbon. According to IUPAC. Calcium  is Ca
        :param charge: the formal atomic charge of the atom
        """
        self.x_coord_list.append(x)
        self.y_coord_list.append(y)
        self.z_coord_list.append(z)
        self.atom_id_list.append(serial_number)
        self.alt_loc_list.append(alternative_location_id)
        self.occupancy_list.append(occupancy)
        self.b_factor_list.append(temperature_factor)
        ## Now add the group level data
        self.current_group.atom_name_list.append(atom_name)
        self.current_group.charge_list.append(charge)
        self.current_group.element_list.append(element)


    def set_chain_info(self, chain_id, chain_name, num_groups):
        """Set the chain information.
        :param chain_id: the asym chain id from mmCIF
        :param chain_name: the auth chain id from mmCIF
        :param num_groups: the number of groups this chain has
        """
        self.chain_id_list.append(chain_id)
        self.chain_name_list.append(chain_name)
        self.groups_per_chain.append(num_groups)


    def set_entity_info(self, chain_indices, sequence, description, entity_type):
        """Set the entity level information for the structure.
        :param chain_indices: the indices of the chains for this entity
        :param sequence: the one letter code sequence for this entity
        :param description: the description for this entity
        :param entity_type: the entity type (polymer,non-polymer,water)
        """
        self.entity_list.append(make_entity_dict(chain_indices,sequence,description,entity_type))


    def set_group_info(self, group_name, group_number, insertion_code,
                       group_type, atom_count, bond_count, single_letter_code,
                       sequence_index, secondary_structure_type):
        """Set the information for a group
        :param group_name: the name of this group,e.g. LYS
        :param group_number: the residue number of this group
        :param insertion_code: the insertion code for this group
        :param group_type: a string indicating the type of group (as found in the chemcomp dictionary.
        Empty string if none available.
        :param atom_count: the number of atoms in the group
        :param bond_count: the number of unique bonds in the group
        :param single_letter_code: the single letter code of the group
        :param sequence_index: the index of this group in the sequence defined by the enttiy
        :param secondary_structure_type: the type of secondary structure used (types are according to DSSP and
        number to type mappings are defined in the specification)
        """
        # Add the group to the overall list - unless it's the first time round
        if self.current_group is not None:
            self.group_list.append(self.current_group)

        # Add the group level information
        self.group_id_list.append(group_number)
        self.ins_code_list.append(insertion_code)
        self.sequence_index_list.append(sequence_index)
        self.sec_struct_list.append(secondary_structure_type)
        self.current_group = Group()
        self.current_group.group_name = group_name
        self.current_group.group_type = group_type
        self.current_group.single_letter_code = single_letter_code

    def set_model_info(self, model_id, chain_count):
        # FIXME model_id here is meaningles and potentially misleading.
        """Set the information for a model.
        :param model_id: the index for the model
        :param chain_count: the number of chains in the model
        """
        self.chains_per_model.append(chain_count)


    def set_xtal_info(self, space_group, unit_cell):
        """Set the crystallographic information for the structure
        :param space_group: the space group name, e.g. "P 21 21 21"
        :param unit_cell: an array of length 6 with the unit cell parameters in order: a, b, c, alpha, beta, gamma
        """
        self.space_group = space_group
        self.unit_cell = unit_cell

    def set_header_info(self, r_free, r_work, resolution, title,
                        deposition_date, release_date, experimental_methods):
        """Sets the header information.
        :param r_free: the measured R-Free for the structure
        :param r_work: the measure R-Work for the structure
        :param resolution: the resolution of the structure
        :param title: the title of the structure
        :param deposition_date: the deposition date of the structure
        :param release_date: the release date of the structure
        :param experimnetal_methods: the list of experimental methods in the structure
        """
        self.r_free = r_free
        self.r_work = r_work
        self.resolution = resolution
        self.title = title
        self.deposition_date = deposition_date
        self.release_date = release_date
        self.experimental_methods = experimental_methods


    def set_bio_assembly_trans(self, bio_assembly_index, input_chain_indices, input_transform):
        """Set the Bioassembly transformation information. A single bioassembly can have multiple transforms,
        :param bio_assembly_index: the integer index of the bioassembly
        :param input_chain_indices: the list of integer indices for the chains of this bioassembly
        :param input_transformation: the list of doubles for  the transform of this bioassmbly transform"""
        this_bioass = None
        for bioass in self.bio_assembly:
            if bioass['name'] == str(bio_assembly_index):
                this_bioass = bioass
                break
        if not this_bioass:
            this_bioass = {"name": str(bio_assembly_index), 'transformList': []}
        else:
            self.bio_assembly.remove(this_bioass)
        this_bioass['transformList'].append({'chainIndexList':input_chain_indices,'matrix': input_transform})
        self.bio_assembly.append(this_bioass)


    def finalize_structure(self):
        """Any functions needed to cleanup the structure."""
        self.group_list.append(self.current_group)
        group_set = get_unique_groups(self.group_list)
        for item in self.group_list:
            self.group_type_list.append(group_set.index(item))
        self.group_list = [x.convert_to_dict() for x in group_set]


    def set_group_bond(self, atom_index_one, atom_index_two, bond_order):
        """Add bonds within a group.
        :param atom_index_one: the integer atom index (in the group) of the first partner in the bond
        :param atom_index_two: the integer atom index (in the group) of the second partner in the bond
        :param bond_order: the integer bond order
        """
        self.current_group.bond_atom_list.append(atom_index_one)
        self.current_group.bond_atom_list.append(atom_index_two)
        self.current_group.bond_order_list.append(bond_order)


    def set_inter_group_bond(self, atom_index_one, atom_index_two, bond_order):
        """Add bonds between groups.
        :param atom_index_one: the integer atom index (in the structure) of the first partner in the bond
        :param atom_index_two: the integer atom index (in the structure) of the second partner in the bond
        :param bond_order the bond order
        """
        self.bond_atom_list.append(atom_index_one)
        self.bond_atom_list.append(atom_index_two)
        self.bond_order_list.append(bond_order)
