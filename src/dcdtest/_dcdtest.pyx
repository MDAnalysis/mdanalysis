
ctypedef int size_t

cdef extern from "Numeric/arrayobject.h":
    struct PyArray_Descr:
        int type_num, elsize
        char type

    ctypedef class Numeric.ArrayType [object PyArrayObject]:
        cdef char *data
        cdef int nd
        cdef int *dimensions, *strides
        cdef object base
        cdef PyArray_Descr *descr
        cdef int flags

cdef extern from *:
    ctypedef int fio_fd
    ctypedef int fio_size_t
    fio_size_t fio_fseek(fio_fd fd, fio_size_t offset, int whence)

cdef extern from "readdcd.h":
    int read_dcdsubset(fio_fd fd, int natoms, int lowerb, int upperb, float *x, float *y, float *z,
                       float *unitcell, int nfixed, int first, int *freeind,
                       float *fixedcoords, int reverse, int charmm)
    int skip_dcdstep(fio_fd fd, int natoms, int nfixed, int charmm, int numstep)

ctypedef struct dcdhandle:
    fio_fd fd
    fio_size_t header_size
    int natoms
    int nsets
    int setsread
    int istart
    int nsavc
    double delta
    int nfixed
    int *freeind
    float *fixedcoords
    int reverse
    int charmm
    int first
    int with_unitcell

cdef extern from "Python.h":
    void* PyCObject_AsVoidPtr(object o)
    char* PyString_AsString(object o)

cdef extern from "correl.h":
    void copyseries(int frame, char *data, int *strides, float *tempX, float *tempY, float *tempZ, char* datacode, int numdata, int* atomlist, int* atomcounts, int lowerb)

import Numeric

def __read_timecorrel(object self, object atoms, object atomcounts, object format, int sizedata, int skip, int lowerb, int upperb):
    cdef dcdhandle* dcd
    cdef ArrayType atomlist, atomcountslist
    cdef ArrayType data, temp
    cdef float *tempX, *tempY, *tempZ
    cdef int rc
    cdef char* fmtstr

    dcd = <dcdhandle*>PyCObject_AsVoidPtr(self._dcd_C_ptr)
    cdef int numframes
    numframes = dcd.nsets / skip
    cdef int numdata
    numdata = len(format)
    if numdata==0:
        raise Exception("No data requested")
    fmtstr = PyString_AsString(format)
    atomlist = Numeric.array(atoms)
    atomcountslist = Numeric.array(atomcounts)
    cdef int range
    range = upperb - lowerb + 1
    # Create data list
    data = Numeric.zeros((sizedata, numframes), Numeric.Float64)
    temp = Numeric.zeros((3, range), Numeric.Float32)
    tempX = <float*>(temp.data+0*temp.strides[0])
    tempY = <float*>(temp.data+1*temp.strides[0])
    tempZ = <float*>(temp.data+2*temp.strides[0])
    # Jump to the beginning of the trajectory file
    rc = fio_fseek(dcd.fd, dcd.header_size, 0) #FIO_SEEK_SET
    dcd.setsread = 0
    dcd.first = 1

    cdef int index, numskip
    cdef int i, j
    cdef float unitcell[6]
    for i from 0 <= i < numframes:
        if (skip > 1):
            # Check if we have fixed atoms
            # XXX not done
            numskip = skip - (dcd.setsread % skip) - 1
            rc = skip_dcdstep(dcd.fd, dcd.natoms, dcd.nfixed, dcd.charmm, numskip)
            if (rc < 0):
                raise IOError("Error skipping frame from DCD file")
            dcd.setsread = dcd.setsread + numskip
        rc = read_dcdsubset(dcd.fd, dcd.natoms, lowerb, upperb, tempX, tempY, tempZ, unitcell, dcd.nfixed, dcd.first, dcd.freeind, dcd.fixedcoords, dcd.reverse, dcd.charmm)
        dcd.first=0
        dcd.setsread = dcd.setsread + 1
        if (rc < 0):
            raise IOError("Error reading frame from DCD file")
        # Copy into data array based on format
        copyseries(i, data.data, data.strides, tempX, tempY, tempZ, fmtstr, numdata, <int*>atomlist.data, <int*>atomcountslist.data, lowerb);
    return data

def __read_timeseries(object self, object atoms, int skip):
    cdef dcdhandle* dcd
    cdef ArrayType atomlist
    cdef ArrayType coord, temp
    cdef float *tempX, *tempY, *tempZ
    cdef int rc

    dcd = <dcdhandle*>PyCObject_AsVoidPtr(self._dcd_C_ptr)
    cdef int numframes
    numframes = dcd.nsets / skip
    cdef int numatoms 
    numatoms = len(atoms)
    if numatoms==0:
        raise Exception("No atoms passed into __read_timeseries function")
    atomlist = Numeric.array(atoms)
    cdef int lowerb, upperb, range
    lowerb = atoms[0]
    upperb = atoms[-1]
    range = upperb - lowerb + 1
    # Create atom list
    coord = Numeric.zeros((numatoms, numframes, 3), Numeric.Float64)
    temp = Numeric.zeros((3, range), Numeric.Float32)
    tempX = <float*>(temp.data+0*temp.strides[0])
    tempY = <float*>(temp.data+1*temp.strides[0])
    tempZ = <float*>(temp.data+2*temp.strides[0])
    # Jump to the beginning of the trajectory file
    rc = fio_fseek(dcd.fd, dcd.header_size, 0) #FIO_SEEK_SET
    dcd.setsread = 0
    dcd.first = 1

    cdef int index, numskip
    cdef int i, j
    cdef float unitcell[6]
    for i from 0 <= i < numframes:
        if (skip > 1):
            # Check if we have fixed atoms
            # XXX not done
            numskip = skip - (dcd.setsread % skip) - 1
            rc = skip_dcdstep(dcd.fd, dcd.natoms, dcd.nfixed, dcd.charmm, numskip)
            if (rc < 0):
                raise IOError("Error skipping frame from DCD file")
            dcd.setsread = dcd.setsread + numskip
        rc = read_dcdsubset(dcd.fd, dcd.natoms, lowerb, upperb, tempX, tempY, tempZ, unitcell, dcd.nfixed, dcd.first, dcd.freeind, dcd.fixedcoords, dcd.reverse, dcd.charmm)
        dcd.first=0
        dcd.setsread = dcd.setsread + 1
        if (rc < 0):
            raise IOError("Error reading frame from DCD file")
        # Copy into numeric array
        for j from 0 <= j < numatoms:
            index = (<int*>atomlist.data)[j]-lowerb
            (<double*> (coord.data+j*coord.strides[0]+i*coord.strides[1]+0*coord.strides[2]))[0] = tempX[index]
            (<double*> (coord.data+j*coord.strides[0]+i*coord.strides[1]+1*coord.strides[2]))[0] = tempY[index]
            (<double*> (coord.data+j*coord.strides[0]+i*coord.strides[1]+2*coord.strides[2]))[0] = tempZ[index]
    return coord
