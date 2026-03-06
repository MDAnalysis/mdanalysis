from mmtf import converters
from .decoders import delta_decode,run_length_decode
from .encoders import delta_encode,run_length_encode

class DeltaRecursiveFloat():
    """Covert an array of floats to integers, perform delta
    encoding and then use recursive indexing to store as 2
    byte integers in a byte array."""
    @staticmethod
    def decode(in_array, param):
        return converters.convert_ints_to_floats(
            delta_decode(
                converters.recursive_index_decode(
                    converters.convert_bytes_to_ints(in_array,2))),param)
    @staticmethod
    def encode(in_array, param):
        return converters.convert_ints_to_bytes(
            converters.recursive_index_encode(
                delta_encode(
                    converters.convert_floats_to_ints(in_array,param))),2)

class RunLengthFloat():
    """Covert an array of floats to integers, perform run-length
    encoding and then store as four byte integers in a byte array."""
    @staticmethod
    def decode(in_array, param):
        return converters.convert_ints_to_floats(
                run_length_decode(
                    converters.convert_bytes_to_ints(in_array,4)),param)
    @staticmethod
    def encode(in_array,param):
        return converters.convert_ints_to_bytes(
            run_length_encode(
                converters.convert_floats_to_ints(in_array,param)),4)

class RunLengthDeltaInt():
    """Delta encode an array of integers and then perform run-length
    encoding on this and then store as four byte integers in a byte array."""
    @staticmethod
    def decode(in_array, param):
        return delta_decode(
            run_length_decode(
                converters.convert_bytes_to_ints(in_array, 4)))
    @staticmethod
    def encode(in_array,param):
        return converters.convert_ints_to_bytes(
            run_length_encode(
                delta_encode(in_array)),4)

class RunLengthChar():
    """Convert chars to integers and run-length encoode these and then store as
    four byte integers in a byte array."""
    @staticmethod
    def decode(in_array, param):
        return converters.convert_ints_to_chars(
            run_length_decode(
                converters.convert_bytes_to_ints(in_array, 4)))
    @staticmethod
    def encode(in_array, param):
        return converters.convert_ints_to_bytes(
            run_length_encode(
                converters.convert_chars_to_ints(in_array)),4)

class EncodeString():
    """Convert strings to set length byte arrays (in this case four). If
    a string is of lenght less than four a null byte is used instead."""
    @staticmethod
    def decode(in_array,param):
        return converters.decode_chain_list(in_array)
    @staticmethod
    def encode(in_array,param):
        return converters.encode_chain_list(in_array)


class ByteToInt():
    """Convert integers to single bytes and store in byte array."""
    @staticmethod
    def decode(in_array,param):
        return converters.convert_bytes_to_ints(in_array, 1)
    @staticmethod
    def encode(in_array,param):
        return converters.convert_ints_to_bytes(in_array,1)


class FourByteToInt():
    """Convert integers to four bytes and store in byte array."""
    @staticmethod
    def decode(in_array,param):
        return converters.convert_bytes_to_ints(in_array, 4)
    @staticmethod
    def encode(in_array, param):
        return converters.convert_ints_to_bytes(in_array, 4)

