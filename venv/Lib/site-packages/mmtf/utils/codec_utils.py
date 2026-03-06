import struct

import mmtf
import mmtf.utils.constants


def parse_header(input_array):
    """Parse the header and return it along with the input array minus the header.
    :param input_array the array to parse
    :return the codec, the length of the decoded array, the parameter and the remainder
    of the array"""
    codec = struct.unpack(mmtf.utils.constants.NUM_DICT[4], input_array[0:4])[0]
    length = struct.unpack(mmtf.utils.constants.NUM_DICT[4], input_array[4:8])[0]
    param = struct.unpack(mmtf.utils.constants.NUM_DICT[4], input_array[8:12])[0]
    return codec,length,param,input_array[12:]


def add_header(input_array, codec, length, param):
    """Add the header to the appropriate array.
    :param the encoded array to add the header to
    :param the codec being used
    :param the length of the decoded array
    :param the parameter to add to the header
    :return the prepended encoded byte array"""
    return struct.pack(mmtf.utils.constants.NUM_DICT[4], codec) + \
           struct.pack(mmtf.utils.constants.NUM_DICT[4], length) + \
           struct.pack(mmtf.utils.constants.NUM_DICT[4], param) + input_array