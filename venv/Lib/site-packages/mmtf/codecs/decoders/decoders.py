def run_length_decode(in_array):
    """A function to run length decode an int array.
    :param in_array: the input array of integers
    :return the decoded array"""
    switch=False
    out_array=[]
    for item in in_array:
        if switch==False:
            this_item = item
            switch=True
        else:
            switch=False
            out_array.extend([this_item]*int(item))
    return out_array

def delta_decode(in_array):
    """A function to delta decode an int array.

    :param in_array: the input array of integers
    :return the decoded array"""
    if len(in_array) == 0:
        return []
    this_ans = in_array[0]
    out_array = [this_ans]
    for i in range(1, len(in_array)):
        this_ans += in_array[i]
        out_array.append(this_ans)
    return out_array