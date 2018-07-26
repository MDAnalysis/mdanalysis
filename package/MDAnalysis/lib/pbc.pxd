from cymem.cymem cimport Pool

cdef class PBC:
    # Cymem Memory pool, takes care of heap memory allocation
    cdef Pool mem
    # Triclinic vector representation of box
    cdef float** box
    # Full box diagonal
    cdef float* fbox_diag
    # Half box diagonal
    cdef float* hbox_diag
    # Minus half box diagonal
    cdef float* mhbox_diag
    # Maximum cutoff squared
    cdef float max_cutoff2
    # What type of periodic boundaries [ORTHO, TRICLINIC, NONE]
    cdef public int pbc_type
    # triclinic vectors to search for minimum image
    cdef int ntric_vec
    cdef float **tric_vec

    cdef void define_box(self, float[:] box)
    cdef void ortho_define(self, float[:] box)
    cdef void triclinic_define(self, float[:] box)
    cdef void find_triclinic_shifts(self)
    cdef float calc_max_cutoff2(self)

cdef void minimum_image(float* dx, PBC pbc)
