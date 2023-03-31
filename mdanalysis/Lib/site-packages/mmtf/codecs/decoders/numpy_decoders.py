import numpy

def delta_decode(in_array):
    """A function to delta decode an int array.

    :param in_array: the input array of integers
    :return the decoded array"""
    return in_array.cumsum()


def run_length_decode(in_array):
    """A function to run length decode an int array.

    :param in_array: the input array of integers
    :return the decoded array"""
    switch=False
    out_array=[]
    in_array = in_array.tolist()
    for item in in_array:
        if switch==False:
            this_item = item
            switch=True
        else:
            switch=False
            out_array.extend([this_item]*int(item))
    return numpy.asarray(out_array, dtype=numpy.int32)