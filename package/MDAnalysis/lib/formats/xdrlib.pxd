from libc.stdint cimport int64_t
from libc.stdio cimport SEEK_SET

cdef extern from 'include/xdrfile.h':
    ctypedef struct XDRFILE:
        pass

    XDRFILE* xdrfile_open (char * path, char * mode)
    int xdrfile_close (XDRFILE * xfp)
    int xdr_seek(XDRFILE *xfp, int64_t pos, int whence)
    ctypedef float matrix[3][3]
    ctypedef float rvec[3]


cdef extern from 'include/xdrfile_xtc.h':
    int read_xtc_natoms(char * fname, int * natoms)
    int read_xtc(XDRFILE * xfp, int natoms, int * step, float * time, matrix box,
                 rvec * x, float * prec)
    int write_xtc(XDRFILE * xfp, int natoms, int step, float time, matrix box,
                  rvec * x, float prec)
    int read_xtc_n_frames(char *fn, int *n_frames, int *est_nframes, int64_t **offsets);


cdef extern from 'include/xdrfile_trr.h':
    int read_trr_natoms(char *fname, int *natoms)
    int read_trr(XDRFILE *xfp, int natoms, int *step, float *time, float *_lambda,
                 matrix box, rvec *x, rvec *v, rvec *f, int *has_prop)
    int write_trr(XDRFILE *xfp, int natoms, int step, float time, float _lambda,
                  matrix box, rvec *x, rvec *v, rvec *f)


cdef enum:
    EOK = 0
    EHEADER = 1
    ESTRING = 2
    EDOUBLE = 3
    EINTEGER = 4
    EFLOAT = 5
    EUNSIGNEDINTEGER = 6
    ECOMPRESSION = 7
    ECLOSE = 8
    EMAGIC = 9
    EMEMORY = 10
    EENDOFFILE = 11
    ENOTFOUND = 11
