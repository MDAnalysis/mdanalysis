from mmtf.utils.codec_utils import parse_header,add_header
from .codecs import DeltaRecursiveFloat,RunLengthFloat,RunLengthDeltaInt,RunLengthChar,EncodeString,ByteToInt,FourByteToInt

codec_dict = {10: DeltaRecursiveFloat,
              9: RunLengthFloat,
              8: RunLengthDeltaInt,
              6: RunLengthChar,
              5: EncodeString,
              2: ByteToInt,
              4: FourByteToInt}

def decode_array(input_array):
    """Parse the header of an input byte array and then decode using the input array,
    the codec and the appropirate parameter.

    :param input_array: the array to be decoded
    :return the decoded array"""
    codec, length, param, input_array = parse_header(input_array)
    return codec_dict[codec].decode(input_array, param)


def encode_array(input_array, codec, param):
    """Encode the array using the method and then add the header to this array.

    :param input_array: the array to be encoded
    :param codec: the integer index of the codec to use
    :param param: the integer parameter to use in the function
    :return an array with the header added to the fornt"""
    return add_header(codec_dict[codec].encode(input_array, param), codec, len(input_array), param)