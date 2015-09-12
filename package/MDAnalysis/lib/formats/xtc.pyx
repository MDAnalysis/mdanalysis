cimport numpy as np
cimport cython
cimport xdrlib

import cython
import numpy as np
from os.path import exists
from collections import namedtuple

from cython_util cimport ptr_to_ndarray

np.import_array()

ctypedef np.float32_t DTYPE_T
DTYPE = np.float32
cdef int DIMS = 3

XTCFrame = namedtuple('XTCFrame', 'x box step time prec')

cdef class XTCFile:
    """Python wrapper for gromacs-xtc format

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
    >>> import xtc
    >>> with xtc.XTCFile('foo.xtc') as f:
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
        self._has_offsets = False
        self.open(self.fname, mode)

    def __dealloc__(self):
        self.close()

    def close(self):
        """Close the open XTC file

        Raises
        ------
        IOError
            If the XTC file can't be closed for some reason
        """
        cdef int res = 1
        if self.is_open:
            res = xdrlib.xdrfile_close(self.xfp)
            self.is_open = False
            if res != 0:
                raise IOError('Couldn\'t close file: {}'.format(self.fname))
        return True

    def open(self, fname, mode):
        """Open a XTC file

        If another XTC file is currently opened it will be closed

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

            return_code = xdrlib.read_xtc_natoms(fname, &self.n_atoms)

            if return_code != 0:
                raise IOError('XTC read error: {}'.format(return_code))
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

    def tell(self):
        """Get current frame"""
        return self.current_frame

    def calc_offsets(self):
        if not self.is_open:
            return np.array([])
        cdef int n_frames = 0
        cdef int est_nframes = 0
        cdef xdrlib.int64_t* offsets = NULL
        ok = xdrlib.read_xtc_n_frames(self.fname, &n_frames, &est_nframes, &offsets);
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

    def read(self):
        """Read one frame in the XTC file

        Returns
        -------
        xyz: ndarray, shape=(n_atoms, 3)
            cartesion coordinates
        box: ndarray, shape=(3, 3)
            Box vectors for this frame
        int
            current step number, 1 indexed
        float
            current time
        float
            precision used to save file

        Raises
        ------
        RuntimeError
            Something must have happened reading the file
        """
        if self.reached_eof:
            raise RuntimeError('Reached last frame in XTC, seek to 0')
        if not self.is_open:
            raise RuntimeError('No file opened')
        if self.mode != 'r':
            raise RuntimeError('File opened in mode: {}. Reading only allow '
                               'in mode "r"'.format('self.mode'))

        cdef int return_code = 1
        cdef int step
        cdef float time, prec

        cdef np.ndarray xyz = np.empty((self.n_atoms, DIMS), dtype=DTYPE)
        cdef np.ndarray box = np.empty((DIMS, DIMS), dtype=DTYPE)

        return_code = xdrlib.read_xtc(self.xfp, self.n_atoms, <int*> &step,
                                      &time, <xdrlib.matrix>box.data,
                                      <xdrlib.rvec*>xyz.data, <float*> &prec)
        if return_code != xdrlib.EOK and return_code != xdrlib.EENDOFFILE:
            raise RuntimeError('XTC Read Error occured: {}'.format(return_code))

        if return_code == xdrlib.EENDOFFILE:
            self.reached_eof = True
            raise StopIteration

        if return_code == xdrlib.EOK:
            self.current_frame += 1
        return XTCFrame(xyz, box, step, time, prec)

    def write(self, xyz, box, int step, float time, float prec=1000.0):
        """write one frame into XTC file

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
        pref: float (optional)
            precision used to save file

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

        cdef DTYPE_T[:, ::1] xyz_view = np.ascontiguousarray(xyz, dtype=DTYPE)
        cdef DTYPE_T[:, ::1] box_view = np.ascontiguousarray(box, dtype=DTYPE)

        if self.current_frame == 0:
            self.n_atoms = xyz.shape[0]
            self.box = box
        else:
            if self.n_atoms != xyz.shape[0]:
                raise ValueError('Previous frames contained {} atoms. You '
                                 'are trying to write {} atoms.'.format(
                                     self.n_atoms, xyz.shape[1]))
            if not np.allclose(self.box, box, rtol=0, atol=8):
                raise ValueError('Previous frames contained {} box. You '
                                 'are trying to write {} box.'.format(
                                     self.box, box))

        cdef int return_code
        return_code = xdrlib.write_xtc(self.xfp, self.n_atoms, step, time,
                                       <xdrlib.matrix>&box_view[0, 0],
                                       <xdrlib.rvec*>&xyz_view[0, 0], prec)
        if return_code != xdrlib.EOK:
            raise RuntimeError('XTC write error: {}'.format(return_code))

        self.current_frame += 1
