cimport numpy as np
cimport cython
cimport xdrlib
from cython_util cimport ptr_to_ndarray

import cython
import numpy as np
from os.path import exists
from collections import namedtuple


np.import_array()

ctypedef np.float32_t DTYPE_T
DTYPE = np.float32
cdef int DIMS = 3
cdef int HASX = 1
cdef int HASV = 2
cdef int HASF = 4

TrrFrame = namedtuple('TrrFrame', 'x v f box step time lmbda hasx hasv hasf')

cdef class TRRFile:
    """Python wrapper for gromacs-trr format

    This class can be similar to the normal file objects in python. The read()
    function will return a frame and all information in it instead of a single
    line. Additionally the context-manager protocoll is supported as well.

    Parameters
    ----------
    fname: str
        The filename to open.
    mode: ('r', 'w')
        The mode in which to open the file, either 'r' read or 'w' write

    Examples
    --------
    >>> import trr
    >>> with trr.TRRFile('foo.trr') as f:
    >>>     for frame in f:
    >>>         xyz, box, step, time, prec = frame
    >>>         print(xyz)
    """
    cdef readonly int n_atoms
    cdef int is_open
    cdef int reached_eof
    cdef xdrlib.XDRFILE *xfp
    cdef readonly fname
    cdef int current_frame
    cdef str mode
    cdef np.ndarray box
    cdef np.ndarray _offsets
    cdef int _has_offsets

    def __cinit__(self, fname, mode='r'):
        self.fname = fname
        self.is_open = False
        self.open(self.fname, mode)

    def __dealloc__(self):
        self.close()

    def close(self):
        """Close the open TRR file

        Raises
        ------
        IOError
            If the TRR file can't be closed for some reason
        """
        cdef int res = 1
        if self.is_open:
            res = xdrlib.xdrfile_close(self.xfp)
            self.is_open = False
            if res != 0:
                raise IOError('Couldn\'t close file: {}'.format(self.fname))
        return True

    def open(self, fname, mode):
        """Open a TRR file

        If another TRR file is currently opened it will be closed

        Parameters
        ----------
        fname: str
            The filename to open.
        mode: ('r', 'w')
            The mode in which to open the file, either 'r' read or 'w' write

        Raises
        ------
        ValueError
            Wrong mode given
        IOError
            Couldn't read the file
        """
        if self.is_open:
            self.close()
        cdef int return_code = 1
        self.fname = fname
        self.n_atoms = 0
        self.reached_eof = False
        self.current_frame = 0

        if mode == 'r':
            opening_mode = b'r'
        elif mode == 'w':
            opening_mode = b'w'
        else:
            raise ValueError('mode must be one of "r" or "w", you '
                             'supplied {}'.format(mode))
        self.mode = mode

        if self.mode == 'r':
            if not exists(self.fname):
                raise IOError('File does not exists: {}'.format(self.fname))

            return_code = xdrlib.read_trr_natoms(fname, &self.n_atoms)

            if return_code != 0:
                raise IOError('TRR read error: {}'.format(return_code))
            if self.n_atoms <= 0:
                raise IOError('Couldn\'t read number of atoms: {}'.format(
                    self.fname))

        self.xfp = xdrlib.xdrfile_open(fname, opening_mode)
        if self.xfp is NULL:
            raise IOError('error opening xtf file: {}'.format(self.fname))
        self.is_open = True

    def __enter__(self):
        """Support context manager"""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Support context manager"""
        self.close()
        # always propagate exceptions forward
        return False

    def __iter__(self):
        self.close()
        self.open(self.fname, self.mode)
        return self

    def __next__(self):
        if self.reached_eof:
            raise StopIteration
        return self.read()

    def __len__(self):
        if not self.is_open:
            raise RuntimeError('No file currently opened')
        return self.offsets.size

    def seek(self, frame):
        """Seek to Frame

        Parameters
        ----------
        frame: int
            wind the file to given frame

        Raises
        ------
        RuntimeError
            If you seek for more frames then are available
        """
        if frame >= self.offsets.size:
            raise RuntimeError('Trying to seek over max number of frames')
        self.reached_eof = False
        ok = xdrlib.xdr_seek(self.xfp, self.offsets[frame], xdrlib.SEEK_SET)
        self.current_frame = frame

    def calc_offsets(self):
        if not self.is_open:
            return np.array([])
        cdef int n_frames = 0
        cdef int est_nframes = 0
        cdef xdrlib.int64_t* offsets = NULL
        ok = xdrlib.read_trr_n_frames(self.fname, &n_frames, &est_nframes, &offsets);
        # the read_xtc_n_frames allocates memory for the offsets with an
        # overestimation. This number is saved in est_nframes and we need to
        # tell the new numpy array about the whole allocated memory to avoid
        # memory leaks.
        cdef np.ndarray dims = np.array([est_nframes], dtype=np.int64)
        # this handles freeing the allocated memory correctly.
        cdef np.ndarray nd_offsets = ptr_to_ndarray(<void*> offsets, dims, np.NPY_INT64)
        return nd_offsets[:n_frames]

    @property
    def offsets(self):
        if not self._has_offsets:
            self._offsets = self.calc_offsets()
            self._has_offsets = True
        return self._offsets

    @offsets.setter
    def offsets(self, offsets):
        self._offsets = offsets
        self._has_offsets = True

    def tell(self):
        """Get current frame"""
        return self.current_frame

    def read(self):
        """Read one frame in the TRR file

        Returns
        -------
        xyz: ndarray, shape=(n_atoms, 3)
            cartesion coordinates
        velocity: ndarray, shape=(n_atoms, 3)
            velocities
        forces: ndarray, shape=(n_atoms, 3)
            forces
        box: ndarray, shape=(3, 3)
            Box vectors for this frame
        int
            current step number, 1 indexed
        float
            current time
        float
            current lambda value

        Raises
        ------
        RuntimeError
            Something must have happened reading the file
        """
        if self.reached_eof:
            raise RuntimeError('Reached last frame in TRR, seek to 0')
        if not self.is_open:
            raise RuntimeError('No file opened')
        if self.mode != 'r':
            raise RuntimeError('File opened in mode: {}. Reading only allow '
                               'in mode "r"'.format('self.mode'))

        cdef int return_code = 1
        cdef int step = 0
        cdef int has_prop = 0
        cdef float time = 0
        cdef float lmbda = 0

        # Use this instead of memviews here to make sure that references are
        # counted correctly
        cdef np.ndarray xyz = np.empty((self.n_atoms, DIMS), dtype=DTYPE)
        cdef np.ndarray velocity = np.empty((self.n_atoms, DIMS), dtype=DTYPE)
        cdef np.ndarray forces = np.empty((self.n_atoms, DIMS), dtype=DTYPE)
        cdef np.ndarray box = np.empty((DIMS, DIMS), dtype=DTYPE)

        return_code = xdrlib.read_trr(self.xfp, self.n_atoms, <int*> &step,
                                      &time, &lmbda, <xdrlib.matrix>box.data,
                                      <xdrlib.rvec*>xyz.data,
                                      <xdrlib.rvec*>velocity.data,
                                      <xdrlib.rvec*>forces.data,
                                      <int*> &has_prop)
        # trr are a bit weird. Reading after the last frame always always
        # results in an integer error while reading. I tried it also with trr
        # produced by different codes.
        if return_code != xdrlib.EOK and return_code != xdrlib.EENDOFFILE \
           and return_code != xdrlib.EINTEGER:
            print self.current_frame
            print step
            print self.reached_eof
            raise RuntimeError('TRR Read Error occured: {}'.format(return_code))

        # In a trr the integer error seems to indicate that the file is ending.
        # There might be corrupted files where this is a legitimate error. But
        # then we just can't read it and stop there which is not to bad.
        if return_code == xdrlib.EENDOFFILE or return_code == xdrlib.EINTEGER:
            self.reached_eof = True
            raise StopIteration

        if return_code == xdrlib.EOK:
            self.current_frame += 1

        has_x = bool(has_prop & HASX)
        has_v = bool(has_prop & HASV)
        has_f = bool(has_prop & HASF)
        return TrrFrame(xyz, velocity, forces, box, step, time, lmbda,
                        has_x, has_v, has_f)

    def write(self, xyz, velocity, forces, box, int step, float time,
              float _lambda, int natoms):
        """write one frame into TRR file.

        Parameters
        ----------
        xyz: ndarray, shape=(n_atoms, 3)
            cartesion coordinates
        box: ndarray, shape=(3, 3)
            Box vectors for this frame
        step: int
            current step number, 1 indexed
        time: float
            current time
        _lambda: float
            current lambda value
        natoms: int
            number of atoms in frame

        Raises
        ------
        RuntimeError
            Couldn't write the file
        ValueError
            The arguments to not match with previous saved frames.
        """
        if self.mode != 'w':
            raise RuntimeError('File opened in mode: {}. Writing only allow '
                               'in mode "w"'.format('self.mode'))

        cdef float* xyz_ptr = NULL
        cdef float* velocity_ptr = NULL
        cdef float* forces_ptr = NULL

        # if any of the xyz/velocity/forces contain any data we use the helper
        # defined here to get a pointer to their first element. This is the only
        # way I know to have a nice pythonic API to the function that can accept
        # array-like inputs or things like None.
        cdef np.ndarray xyz_helper
        cdef np.ndarray velocity_helper
        cdef np.ndarray forces_helper

        if xyz is not None:
            xyz_helper = np.ascontiguousarray(xyz, dtype=DTYPE)
            xyz_ptr = <float*>xyz_helper.data
        if velocity is not None:
            velocity_helper = np.ascontiguousarray(velocity, dtype=DTYPE)
            velocity_ptr = <float*>velocity_helper.data
        if forces is not None:
            forces_helper = np.ascontiguousarray(forces, dtype=DTYPE)
            forces_ptr = <float*>forces_helper.data

        cdef np.ndarray box_helper = np.ascontiguousarray(box, dtype=DTYPE)
        cdef float* box_ptr = <float*>box_helper.data

        if self.current_frame == 0:
            self.box = box
            self.n_atoms = natoms
        else:
            if self.n_atoms != natoms:
                raise ValueError('Previous frames contained {} atoms. You '
                                 'are trying to write {} atoms.'.format(
                                     self.n_atoms, natoms))
            if xyz is not None and self.n_atoms != xyz.shape[0]:
                raise ValueError('Previous frames xyz contained {} atoms. You '
                                 'are trying to write {} atoms.'.format(
                                     self.n_atoms, xyz.shape[0]))
            if velocity is not None and self.n_atoms != velocity.shape[0]:
                raise ValueError('Previous frames velocity contained {} atoms. You '
                                 'are trying to write {} atoms.'.format(
                                     self.n_atoms, velocity.shape[0]))
            if forces is not None and self.n_atoms != forces.shape[0]:
                raise ValueError('Previous frames forces contained {} atoms. You '
                                 'are trying to write {} atoms.'.format(
                                     self.n_atoms, forces.shape[0]))
            if not np.allclose(self.box, box, rtol=0, atol=8):
                raise ValueError('Previous frames contained {} box. You '
                                 'are trying to write {} box.'.format(
                                     self.box, box))

        cdef int return_code
        return_code = xdrlib.write_trr(self.xfp, self.n_atoms, step, time,
                                       _lambda, <xdrlib.matrix> box_ptr,
                                       <xdrlib.rvec*> xyz_ptr,
                                       <xdrlib.rvec*> velocity_ptr,
                                       <xdrlib.rvec*> forces_ptr)
        if return_code != xdrlib.EOK:
            raise RuntimeError('TRR write error: {}'.format(return_code))

        self.current_frame += 1
