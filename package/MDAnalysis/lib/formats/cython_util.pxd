cimport numpy as np
cdef np.ndarray ptr_to_ndarray(void* data_ptr, np.int64_t[:] dim, int data_type)
