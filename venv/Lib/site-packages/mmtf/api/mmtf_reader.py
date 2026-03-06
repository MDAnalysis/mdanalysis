from mmtf.codecs import decode_array
from mmtf.utils import decoder_utils
import sys


class MMTFDecoder(object):
    """Class to decode raw mmtf data into a parsed data model that can be fed into other data model"""
    model_counter = 0
    chain_counter = 0
    group_counter = 0
    atom_counter = 0


    def get_coords(self):
        """Utility function to get the coordinates as a single list of tuples."""
        out_list = []
        for i in range(len(self.x_coord_list)):
            out_list.append((self.x_coord_list[i],self.y_coord_list[i],self.z_coord_list[i],))
        return out_list

    def get_bonds(self):
        """Utility function to get all the inter group bonds for the structure in pairs."""
        return decoder_utils.get_bonds(self)

    def decode_data(self, input_data):
        """Function to decode the input data and place it onto the class.
        :param input_data: the input data as a dict"""
        self.group_type_list = decode_array(input_data["groupTypeList"])
        self.x_coord_list = decode_array(input_data["xCoordList"])
        self.y_coord_list = decode_array(input_data["yCoordList"])
        self.z_coord_list = decode_array(input_data["zCoordList"])
        if "bFactorList" in input_data:
            self.b_factor_list = decode_array(input_data["bFactorList"])
        else:
            self.b_factor_list = []
        if "occupancyList" in input_data:
            self.occupancy_list = decode_array(input_data["occupancyList"])
        else:
            self.occupancy_list = []
        if "atomIdList" in input_data:
            self.atom_id_list = decode_array(input_data["atomIdList"])
        else:
            self.atom_id_list = []
        if "altLocList" in input_data:
            self.alt_loc_list = decode_array(input_data["altLocList"])
        else:
            self.alt_loc_list = []
        if "insCodeList" in input_data:
            self.ins_code_list = decode_array(input_data["insCodeList"])
        else:
            self.ins_code_list = []
        self.group_id_list = decode_array(input_data["groupIdList"])
        self.group_list = input_data["groupList"]
        if "sequenceIndexList" in input_data:
            self.sequence_index_list = decode_array(input_data["sequenceIndexList"])
        else:
            self.sequence_index_list = []
        self.chains_per_model = input_data["chainsPerModel"]
        self.groups_per_chain = input_data["groupsPerChain"]
        if "chainNameList" in input_data:
            self.chain_name_list = decode_array(input_data["chainNameList"])
        else:
            self.chain_name_list = []
        self.chain_id_list = decode_array(input_data["chainIdList"])
        if "spaceGroup" in input_data:
            self.space_group = input_data["spaceGroup"]
        else:
            self.space_group = None
        if "bondAtomList" in input_data:
            self.bond_atom_list = decode_array(input_data["bondAtomList"])
        else:
            self.bond_atom_list = None
        if "bondOrderList" in input_data:
            self.bond_order_list = decode_array(input_data["bondOrderList"])
        else:
            self.bond_order_list = None
        if "mmtfVersion" in input_data:
            self.mmtf_version = input_data["mmtfVersion"]
        else:
            self.mmtf_version = None
        if "mmtfProducer" in input_data:
            self.mmtf_producer = input_data["mmtfProducer"]
        else:
            self.mmtf_producer = None
        if "structureId" in input_data:
            self.structure_id = input_data["structureId"]
        else:
            self.structure_id = None
        if "title" in input_data:
            self.title = input_data["title"]
        else:
            self.title = None 
        if "experimentalMethods" in input_data:
            self.experimental_methods = input_data["experimentalMethods"]
        else:
            self.experimental_methods = None
        if "depositionDate" in input_data:
            self.deposition_date = input_data["depositionDate"]
        else:
            self.deposition_date = None
        if "releaseDate" in input_data:
            self.release_date = input_data["releaseDate"]
        else:
            self.release_date = None
        if "entityList" in input_data:
            self.entity_list = input_data["entityList"]
        else:
            self.entity_list = []
        if "bioAssemblyList" in input_data:
            self.bio_assembly = input_data["bioAssemblyList"]
        else:
            self.bio_assembly = []
        if "rFree" in input_data:
            self.r_free = input_data["rFree"]
        else:
            self.r_free = None
        if "rWork" in input_data:
            self.r_work = input_data["rWork"]
        else:
            self.r_work = None
        if "resolution" in input_data:
            self.resolution = input_data["resolution"]
        else:
            self.resolution = None
        if "unitCell" in input_data:
            self.unit_cell = input_data["unitCell"]
        else:
            self.unit_cell = None
        if "secStructList" in input_data:
            self.sec_struct_list = decode_array(input_data["secStructList"])
        # Now all the numbers to defien the
        self.num_bonds = int(input_data["numBonds"])
        self.num_chains = int(input_data["numChains"])
        self.num_models = int(input_data["numModels"])
        self.num_atoms = int(input_data["numAtoms"])
        self.num_groups = int(input_data["numGroups"])


    def pass_data_on(self, data_setters):
        """Write the data from the getters to the setters.

        :param data_setters: a series of functions that can fill a chemical
        data structure
        :type data_setters: DataTransferInterface
        """
        data_setters.init_structure(self.num_bonds, len(self.x_coord_list), len(self.group_type_list),
                                    len(self.chain_id_list), len(self.chains_per_model), self.structure_id)
        decoder_utils.add_entity_info(self, data_setters)
        decoder_utils.add_atomic_information(self, data_setters)
        decoder_utils.add_header_info(self, data_setters)
        decoder_utils.add_xtalographic_info(self, data_setters)
        decoder_utils.generate_bio_assembly(self, data_setters)
        decoder_utils.add_inter_group_bonds(self, data_setters)
        data_setters.finalize_structure()
