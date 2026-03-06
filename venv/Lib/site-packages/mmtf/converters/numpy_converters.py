from __future__ import division
import numpy
import mmtf
import mmtf.utils.constants


def convert_bytes_to_ints(in_bytes, num):
    """Convert a byte array into an integer array. The number of bytes forming an integer
    is defined by num

    :param in_bytes: the input bytes
    :param num: the number of bytes per int
    :return the integer array"""
    dt = numpy.dtype('>i' + str(num))
    return numpy.frombuffer(in_bytes, dt)

def decode_chain_list(in_bytes):
    """Convert a list of bytes to a list of strings. Each string is of length mmtf.CHAIN_LEN

    :param in_bytes: the input bytes
    :return the decoded list of strings"""
    bstrings = numpy.frombuffer(in_bytes, numpy.dtype('S' + str(mmtf.utils.constants.CHAIN_LEN)))
    return [s.decode("ascii").strip(mmtf.utils.constants.NULL_BYTE) for s in bstrings]

def convert_ints_to_floats(in_ints, divider):
    """Convert integers to floats by division.
    :param in_ints: the integer array
    :param divider: the divider
    :return the array of floats produced"""
    return (in_ints.astype(numpy.float64) / divider)

def recursive_index_decode(int_array, max=32767, min=-32768):
    """Unpack an array of integers using recursive indexing.

    :param int_array: the input array of integers
    :param max: the maximum integer size
    :param min: the minimum integer size
    :return the array of integers after recursive index decoding"""
    out_arr = []
    decoded_val = 0
    for item in int_array.tolist():
        if item==max or item==min:
            decoded_val += item
        else:
            decoded_val += item
            out_arr.append(decoded_val)
            decoded_val = 0
    return numpy.asarray(out_arr,dtype=numpy.int32)

