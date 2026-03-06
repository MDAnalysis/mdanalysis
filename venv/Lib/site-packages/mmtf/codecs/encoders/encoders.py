def run_length_encode(in_array):
    """A function to run length decode an int array.

    :param in_array: the inptut array of integers
    :return the encoded integer array"""
    if(len(in_array)==0):
        return []
    curr_ans = in_array[0]
    out_array = [curr_ans]
    counter = 1
    for in_int in in_array[1:]:
        if in_int == curr_ans:
            counter+=1
        else:
            out_array.append(counter)
            out_array.append(in_int)
            curr_ans = in_int
            counter = 1
    # Add the final counter
    out_array.append(counter)
    return out_array

def delta_encode(in_array):
    """A function to delta decode an int array.

    :param in_array: the inut array to be delta encoded
    :return the encoded integer array"""
    if(len(in_array)==0):
        return []
    curr_ans = in_array[0]
    out_array = [curr_ans]
    for in_int in in_array[1:]:
        out_array.append(in_int-curr_ans)
        curr_ans = in_int
    return out_array