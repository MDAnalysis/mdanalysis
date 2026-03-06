import unittest

import msgpack
import numpy

from mmtf import fetch,parse,parse_gzip, converters
from mmtf.api.default_api import ungzip_data,write_mmtf,MMTFDecoder,_internet_on
from mmtf.codecs import encoders
from mmtf.utils.codec_utils import parse_header
from mmtf.utils.constants import BASE_URL
from mmtf.codecs.default_codec import codec_dict
from mmtf.codecs.decoders import numpy_decoders as decoders


def run_all(unit_test, encoded_data, decoded_data, param, codec_id):
    """Test that a given codec can work in the forward backward and round trip both ways."""
    try:
        unit_test.assertEqual(codec_dict[codec_id].decode(encoded_data, param).tolist(), decoded_data.tolist())
    except:
        unit_test.assertEqual(codec_dict[codec_id].decode(encoded_data, param), decoded_data.tolist())
    try:
        unit_test.assertEqual(
            codec_dict[codec_id].decode(codec_dict[codec_id].encode(decoded_data, param),
                                               param).tolist(), decoded_data.tolist())
    except:
        unit_test.assertEqual(codec_dict[codec_id].decode(codec_dict[codec_id].encode(decoded_data, param),
                                               param), decoded_data.tolist())
    unit_test.assertEqual(codec_dict[codec_id].encode(decoded_data, param), encoded_data)
    unit_test.assertEqual(codec_dict[codec_id].encode(codec_dict[codec_id].decode(encoded_data, param),
                                                             param), encoded_data)

class CodecTest(unittest.TestCase):
    def test_delt_rec_float(self):
        test_data = b'\x7f\xffD\xab\x01\x8f\xff\xca'
        output_data = numpy.array([50.346, 50.745, 50.691])
        run_all(self, test_data, output_data, 1000, 10)

    def test_run_len_float(self):
        test_data = b'\x00\x00\x00d\x00\x00\x00\x03'
        output_data = numpy.array([1.00,1.00,1.00])
        run_all(self, test_data, output_data, 100, 9)

    def test_run_len_delta_int(self):
        test_data = b'\x00\x00\x00\x01\x00\x00\x00\x07'
        output_data = numpy.array([1,2,3,4,5,6,7])
        run_all(self, test_data, output_data, 0, 8)

    def test_run_len_char(self):
        test_data = b'\x00\x00\x00\x41\x00\x00\x00\x04'
        output_data = numpy.array(["A","A","A","A"])
        run_all(self, test_data, output_data, 0, 6)

    def test_enc_str(self):
        test_data = b'B\x00\x00\x00A\x00\x00\x00C\x00\x00\x00A\x00\x00\x00A\x00\x00\x00A\x00\x00\x00'
        output_data =  numpy.array(["B","A","C","A","A","A"])
        run_all(self, test_data, output_data, 0, 5)

    def test_byte_to_int(self):
        test_data =  b'\x07\x06\x06\x07\x07'
        output_data = numpy.array([7,6,6,7,7])
        run_all(self, test_data, output_data, 0, 2)

    def test_four_byte_int(self):
        test_data = b'\x00\x00\x00\x01\x00\x02\x00\x01\x00\x00\x00\x00\x00\x00\x00\x02'
        output_data = numpy.array([1, 131073, 0, 2])
        run_all(self, test_data, output_data, 0, 4)

class DecoderTests(unittest.TestCase):
    def test_run_length_decode(self):
        input_data = numpy.array([15,3,100,2,111,4,10000,6])
        output_data_test = [15,15,15,100,100,111,111,111,111,10000,10000,10000,10000,10000,10000]
        output_data = decoders.run_length_decode(input_data).tolist()
        self.assertEqual(output_data, output_data_test)

    def test_empty_run_length_decode(self):
        input_data = numpy.array([])
        output_data_test = []
        output_data = decoders.run_length_decode(input_data).tolist()
        self.assertEqual(output_data, output_data_test)


    def test_delta_decode(self):
        input_data = numpy.asarray([15,3,100,-1,11,4],dtype=numpy.int32)
        output_data_test = [15,18,118,117,128,132]
        output_data = decoders.delta_decode(input_data).tolist()
        self.assertEqual(output_data, output_data_test)

    def test_empty_delta_decode(self):
        input_data = numpy.asarray([],dtype=numpy.int32)
        output_data_test = []
        output_data = decoders.delta_decode(input_data).tolist()
        self.assertEqual(output_data, output_data_test)

class EncoderTests(unittest.TestCase):
    def test_run_length_encode(self):
        output_data_test = [15, 3, 100, 2, 111, 4, 10000, 6]
        input_data = [15, 15, 15, 100, 100, 111, 111, 111, 111, 10000, 10000, 10000, 10000, 10000, 10000]
        output_data = encoders.run_length_encode(input_data)
        self.assertEqual(output_data, output_data_test)

    def test_empty_run_length_encode(self):
        input_data = []
        output_data_test = []
        output_data = encoders.run_length_encode(input_data)
        self.assertEqual(output_data, output_data_test)

    def test_delta_encode(self):
        output_data_test = [15, 3, 100, -1, 11, 4]
        input_data = [15, 18, 118, 117, 128, 132]
        output_data = encoders.delta_encode(input_data)
        self.assertEqual(output_data, output_data_test)

    def test_empty_delta_encode(self):
        input_data = []
        output_data_test = []
        output_data = encoders.delta_encode(input_data)
        self.assertEqual(output_data, output_data_test)


class ConverterTests(unittest.TestCase):

    def test_convert_chain_list(self):
        in_bytes = b'A\x00\x00\x00A\x00\x00\x00A\x00\x00\x00A\x00\x00\x00A\x00\x00\x00A\x00\x00\x00'
        out_strings_test =  ["A", "A","A","A","A","A"]
        self.assertEqual(out_strings_test, converters.decode_chain_list(in_bytes))
        self.assertEqual(in_bytes, converters.encode_chain_list(out_strings_test))

    def test_convert_int_to_float(self):
        in_array = numpy.asarray([10001,100203,124542])
        out_array_test = [10.001,100.203,124.542]
        converted = converters.convert_ints_to_floats(in_array, 1000.0).tolist()
        for i in range(len(out_array_test)):
            self.assertAlmostEqual(out_array_test[i], converted[i],places=3)
        self.assertEqual(in_array.tolist(), converters.convert_floats_to_ints(out_array_test, 1000.0))

    def test_recursive_enc(self):
        in_arr = [1,420,32767,120,-32768,34767]
        out_array_test = [1,420,32767,0,120,-32768,0,32767,2000]
        self.assertEqual(out_array_test, converters.recursive_index_encode(in_arr))

    def test_recursive_dec(self):
        in_arr = numpy.asarray([1,420,32767,0,120,-32768,0,32767,2000],dtype=numpy.int32)
        out_array_test = [1,420,32767,120,-32768,34767]
        self.assertEqual(out_array_test, converters.recursive_index_decode(in_arr).tolist())

    def test_convert_one_byte_int(self):
        in_bytes = b'\x07\x06\x06\x07\x07'
        out_array_test = [7,6,6,7,7]
        self.assertEqual(out_array_test, converters.convert_bytes_to_ints(in_bytes,1).tolist())
        self.assertEqual(in_bytes, converters.convert_ints_to_bytes(out_array_test,1))
        self.assertEqual(in_bytes,converters.convert_ints_to_bytes(converters.convert_bytes_to_ints(in_bytes,1),1))

    def test_convert_two_byte_int(self):
        in_bytes = b'\x00\x00\x00\x01\x00\x02\x00\x01\x00\x00\x00\x00\x00\x00\x00\x02'
        out_array_test = [0,1,2,1,0,0,0,2]
        self.assertEqual(out_array_test, converters.convert_bytes_to_ints(in_bytes,2).tolist())
        self.assertEqual(in_bytes, converters.convert_ints_to_bytes(out_array_test, 2))
        self.assertEqual(in_bytes,converters.convert_ints_to_bytes(converters.convert_bytes_to_ints(in_bytes,2),2))


    def test_convert_four_byte_int(self):
        in_bytes = b'\x00\x00\x00\x01\x00\x02\x00\x01\x00\x00\x00\x00\x00\x00\x00\x02'
        out_array_test = [1, 131073, 0, 2]
        self.assertEqual(out_array_test, converters.convert_bytes_to_ints(in_bytes,4).tolist())
        self.assertEqual(in_bytes, converters.convert_ints_to_bytes(out_array_test,4))
        self.assertEqual(in_bytes,converters.convert_ints_to_bytes(converters.convert_bytes_to_ints(in_bytes,4),4))

    def test_parse_header(self):
        in_bytes = b'\x00\x00\x00\x01\x00\x02\x00\x01\x00\x00\x00\x00\x00\x00\x00\x02'
        codec,length,param, bytearray = parse_header(in_bytes)
        self.assertEqual(length,131073)
        self.assertEqual(param,0)
        self.assertEqual(len(bytearray),4)

    def test_convert_int_to_char(self):
        int_array =  [66,63,67]
        out_array_test = ["B", "?","C"]
        self.assertEqual(out_array_test, converters.convert_ints_to_chars(int_array))
        self.assertEqual(int_array, converters.convert_chars_to_ints(out_array_test))

    def test_decoder(self):
        decoded = parse("mmtf/tests/testdatastore/4CUP.mmtf")

    def test_gz_decoder(self):
        decoded = parse_gzip("mmtf/tests/testdatastore/4CUP.mmtf.gz")

    def test_round_trip(self):
        decoded = parse("mmtf/tests/testdatastore/4CUP.mmtf")
        packed = decoded.get_msgpack()
        decoded.decode_data(msgpack.unpackb(packed))

    def test_gzip_open(self):
        with open("mmtf/tests/testdatastore/4CUP.mmtf.gz","rb") as fh:
            ungzip_data(fh.read())

    def test_fetch(self):
        if _internet_on(BASE_URL):
            decoded = fetch("4CUP")
        else:
            print("Warning - cannot connect to "+BASE_URL)


    def array_eq(self,array_one, array_two):
        import numpy as np
        if [x for x in np.isclose(array_one,array_two) if x]:
            return True
        else:
            try:
                if not array_one and not array_two:
                    return True
            except ValueError:
                pass
            print(array_one)
            print(array_two)
            print("Arrays not equal")
            return False

    def char_arr_eq(self,array_one, array_two):
        import numpy as np
        return np.array_equal(array_one,array_two)

    def dict_list_equal(self,list_one,list_two):
        list_one = sorted(list_one, key=lambda x:sorted(x.keys()))
        list_two = sorted(list_two, key=lambda x:sorted(x.keys()))
        len_one = len(list_one)
        if len_one != len(list_two):
            self.assertTrue(False,"Lists of different lengths")
        for i in range(len_one):
            if list_one[i]!=list_two[i]:
                print(list_one[i])
                print(list_two[i])
            self.assertTrue(list_one[i]==list_two[i])

    def iterate(self, data_one, data_two):
        chain_ind = 0
        group_ind = 0
        atom_ind_one = 0
        atom_ind_two = 0
        for model in data_one.chains_per_model:
            for chain in range(model):
                for group in range(data_one.groups_per_chain[chain_ind]):
                    self.char_arr_eq(data_one.group_list[data_one.group_type_list[group_ind]]["atomNameList"],
                                  data_two.group_list[data_two.group_type_list[group_ind]]["atomNameList"])
                    self.char_arr_eq(data_one.group_list[data_one.group_type_list[group_ind]]["elementList"],
                                     data_two.group_list[data_two.group_type_list[group_ind]]["elementList"])
                    self.array_eq(data_one.group_list[data_one.group_type_list[group_ind]]["bondOrderList"],
                                     data_two.group_list[data_two.group_type_list[group_ind]]["bondOrderList"])
                    self.array_eq(data_one.group_list[data_one.group_type_list[group_ind]]["bondAtomList"],
                                     data_two.group_list[data_two.group_type_list[group_ind]]["bondAtomList"])
                    self.array_eq(data_one.group_list[data_one.group_type_list[group_ind]]["formalChargeList"],
                                     data_two.group_list[data_two.group_type_list[group_ind]]["formalChargeList"])
                    self.assertEqual(data_one.group_list[data_one.group_type_list[group_ind]]["groupName"],
                                     data_two.group_list[data_two.group_type_list[group_ind]]["groupName"])
                    self.assertEqual(data_one.group_list[data_one.group_type_list[group_ind]]["singleLetterCode"],
                                     data_two.group_list[data_two.group_type_list[group_ind]]["singleLetterCode"])
                    self.assertEqual(data_one.group_list[data_one.group_type_list[group_ind]]["chemCompType"],
                                     data_two.group_list[data_two.group_type_list[group_ind]]["chemCompType"])
                    group_ind+=1
                chain_ind+=1
        return True

    def check_equal(self, data_one, data_two):
        self.assertTrue(self.array_eq(data_one.x_coord_list,data_two.x_coord_list))
        self.assertTrue(self.array_eq(data_one.y_coord_list,data_two.y_coord_list))
        self.assertTrue(self.array_eq(data_one.z_coord_list,data_two.z_coord_list))
        self.assertTrue(self.array_eq(data_one.b_factor_list,data_two.b_factor_list))
        self.assertTrue(self.array_eq(data_one.occupancy_list,data_two.occupancy_list))
        self.assertTrue(self.array_eq(data_one.atom_id_list,data_two.atom_id_list))
        self.assertTrue(self.char_arr_eq(data_one.alt_loc_list,data_two.alt_loc_list))
        self.assertTrue(self.char_arr_eq(data_one.ins_code_list,data_two.ins_code_list))
        self.assertTrue(self.array_eq(data_one.group_id_list,data_two.group_id_list))
        self.dict_list_equal(data_one.entity_list,data_two.entity_list)
        self.dict_list_equal(data_one.bio_assembly,data_two.bio_assembly)
        self.assertTrue(self.array_eq(data_one.sequence_index_list,data_two.sequence_index_list))
        self.assertEqual(data_one.chains_per_model, data_two.chains_per_model)
        self.assertEqual(data_one.groups_per_chain, data_two.groups_per_chain)
        self.assertEqual(data_one.chain_name_list, data_two.chain_name_list)
        self.assertEqual(data_one.chain_id_list, data_two.chain_id_list)
        self.assertEqual(data_one.space_group,data_two.space_group)
        self.assertTrue(self.array_eq(data_one.bond_atom_list,data_two.bond_atom_list))
        self.assertTrue(self.array_eq(data_one.bond_order_list,data_two.bond_order_list))
        self.assertEqual(data_one.structure_id,data_two.structure_id)
        self.assertEqual(data_one.title,data_two.title)
        self.assertTrue(self.char_arr_eq(data_one.experimental_methods,data_two.experimental_methods))
        self.assertEqual(data_one.deposition_date,data_two.deposition_date)
        self.assertEqual(data_one.release_date,data_two.release_date)
        self.assertTrue(self.array_eq(data_one.sec_struct_list,data_two.sec_struct_list))
        self.assertEqual(data_one.r_free,data_two.r_free)
        self.assertEqual(data_one.r_work,data_two.r_work)
        self.assertEqual(data_one.resolution,data_two.resolution)
        self.assertEqual(data_one.unit_cell,data_two.unit_cell)
        self.assertEqual(data_one.num_bonds, data_two.num_bonds)
        self.assertEqual(data_one.num_chains, data_two.num_chains)
        self.assertEqual(data_one.num_models, data_two.num_models)
        self.assertEqual(data_one.num_atoms, data_two.num_atoms)
        self.assertEqual(data_one.num_groups, data_two.num_groups)
        self.assertTrue(self.iterate(data_one, data_two))

    def test_round_trip(self):
        data_in = parse_gzip("mmtf/tests/testdatastore/4CUP.mmtf.gz")
        write_mmtf("test.mmtf", data_in, MMTFDecoder.pass_data_on)
        data_rt = parse("test.mmtf")
        self.check_equal(data_in, data_rt)

    def round_trip(self,pdb_id):
        if _internet_on(BASE_URL):
            data_in = fetch(pdb_id)
            write_mmtf(pdb_id+".mmtf", data_in, MMTFDecoder.pass_data_on)
            data_rt = parse(pdb_id+".mmtf")
            self.check_equal(data_in, data_rt)
        else:
            print("Warning - cannot connect to "+BASE_URL)

    def test_round_trip_list(self):
        id_list = [
            #
            "1a1q",
            # // Just added to check
            "9pti",
            # // An entity that has no chain
            "2ja5",
            # // A couple of examples of multiple disulpgide bonds being formed.
            "3zxw",
            "1nty",
            # // A weird residue case
            "2eax",
            # // A Deuterated Structure
            "4pdj",
            # // Weird bioassembly
            "4a1i",
            # // Multi model structure
            "1cdr",
            # // Another weird structure (jose's suggestion)
            "3zyb",
            # //Standard structure
            "4cup",
            # // Weird NMR structure
            "1o2f",
            # // B-DNA structure
            "1bna",
            # // DNA structure
            "4y60",
            # // Sugar structure
            "1skm",
            # // Calpha atom is missing (not marked as calpha)
            "1lpv",
            # // NMR structure with multiple models - one of which has chain missing
            "1msh",
            # // No ATOM records just HETATM records (in PDB). Opposite true for MMCif. It's a D-Peptide.
            "1r9v",
            # // Biosynthetic protein
            "5emg",
            # // Micro heterogenity
            "4ck4",
            # // Ribosome
            "4v5a",
            # // Negative residue numbers
            "5esw",
            # // A tiny example case
            "3njw",
            # // A GFP example with weird seqres records
            "1ema"]
        for pdb_id in id_list:
            self.round_trip(pdb_id)

if __name__ == '__main__':
    unittest.main()
