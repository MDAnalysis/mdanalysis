# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
# Beckstein and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
cimport numpy as np
cimport cython
from cython_util cimport ptr_to_ndarray
from libc.stdint cimport int64_t

from libc.stdio cimport SEEK_SET, SEEK_CUR, SEEK_END
_whence_vals = {"SEEK_SET": SEEK_SET, "SEEK_CUR": SEEK_CUR, "SEEK_END": SEEK_END}

cdef extern from 'include/xdrfile.h':
    ctypedef struct XDRFILE:
        pass

    XDRFILE* xdrfile_open (char * path, char * mode)
    int xdrfile_close (XDRFILE * xfp)
    int xdr_seek(XDRFILE *xfp, int64_t pos, int whence)
    int64_t xdr_tell(XDRFILE *xfp)
    ctypedef float matrix[3][3]
    ctypedef float rvec[3]


cdef extern from 'include/xdrfile_xtc.h':
    int read_xtc_natoms(char * fname, int * natoms)
    int read_xtc(XDRFILE * xfp, int natoms, int * step, float * time, matrix box,
                 rvec * x, float * prec)
    int write_xtc(XDRFILE * xfp, int natoms, int step, float time, matrix box,
                  rvec * x, float prec)



cdef extern from 'include/xdrfile_trr.h':
    int read_trr_natoms(char *fname, int *natoms)
    int read_trr(XDRFILE *xfp, int natoms, int *step, float *time, float *_lambda,
                 matrix box, rvec *x, rvec *v, rvec *f, int *has_prop)
    int write_trr(XDRFILE *xfp, int natoms, int step, float time, float _lambda,
                  matrix box, rvec *x, rvec *v, rvec *f)


cdef extern from 'include/xtc_seek.h':
    int read_xtc_n_frames(char *fn, int *n_frames, int *est_nframes, int64_t **offsets)


cdef extern from 'include/trr_seek.h':
    int read_trr_n_frames(char *fn, int *n_frames, int *est_nframes, int64_t **offsets)


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
    ENOTFOUND = 12

error_message = ['OK', 'header', 'string', 'double', 'integer',
                 'float', 'unsigned integer', 'compression',
                 'close', 'magic', 'memory', 'endoffile',
                 'notfound']

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


cdef class _XDRFile:
    """Base python wrapper for gromacs-xdr formats

    This class can be similar to the normal file objects in python. The read()
    function will return a frame and all information in it instead of a single
    line. Additionally the context-manager protocoll is supported as well.

    Parameters
    ----------
    fname : str
        The filename to open.
    mode : ('r', 'w')
        The mode in which to open the file, either 'r' read or 'w' write

    Raises
    ------
    ValueError
        Wrong mode given
    IOError
        Couldn't read the file

    Note
    ----
    This class can't be initialized use one of the subclasses XTCFile, TRRFile
    """
    cdef readonly int n_atoms
    cdef int is_open
    cdef int reached_eof
    cdef XDRFILE *xfp
    cdef readonly fname
    cdef int current_frame
    cdef str mode
    cdef np.ndarray box
    cdef np.ndarray _offsets
    cdef int _has_offsets

    def __cinit__(self, fname, mode='r'):
        self.fname = fname.encode('utf-8')
        self.is_open = False
        self.open(self.fname, mode)

    def __dealloc__(self):
        self.close()

    def close(self):
        """Close the open XDR file

        Raises
        ------
        IOError
            If the TRR file can't be closed for some reason
        """
        res = 1
        if self.is_open:
            res = xdrfile_close(self.xfp)
            self.is_open = False
            if res != EOK:
                raise IOError('Couldn\'t close file: {}  XDR error = {}'.format(
                    self.fname, error_message[res]))
        # forget old offsets in case we open a different file with the same instance.
        self._has_offsets = False
        return True

    def open(self, fname, mode):
        """Open a XDR file

        If another XDR file is currently opened it will be closed

        Parameters
        ----------
        fname : str
            The filename to open.
        mode : ('r', 'w')
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
                raise IOError('File does not exist: {}'.format(self.fname))

            return_code, self.n_atoms = self._calc_natoms(fname)

            if return_code != EOK:
                raise IOError('XDR read error = {}'.format(
                    error_message[return_code]))
            if self.n_atoms <= 0:
                raise IOError('Couldn\'t read number of atoms: {}'.format(
                    fname))

        self.xfp = xdrfile_open(fname, opening_mode)
        if self.xfp is NULL:
            raise IOError('Error opening XTC/TRR file: {}'.format(self.fname))
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
        """Seek to Frame.

        Please note that this function will generate internal file offsets if
        they haven't been set before. For large file this means the first seek
        can be very slow. Later seeks will be very fast.

        Parameters
        ----------
        frame : int
            seek the file to given frame

        Raises
        ------
        IOError
            If you seek for more frames than are available or if the
            seek fails (the low-level system error is reported).
        """
        cdef int64_t offset
        if frame == 0:
            offset = 0
        else:
            if frame >= self.offsets.size:
                raise IOError('Trying to seek over max number of frames')
            offset = self.offsets[frame]
        self.reached_eof = False
        ok = xdr_seek(self.xfp, offset, SEEK_SET)
        if ok != EOK:
            # errno is the direct system error, not from the XDR library.
            # errno=22 can be cryptic, since it means a wrong 'whence'
            #  parameter but it will also be issued when the resulting file
            #  offset is negative, or if the offset overflows the filesystem
            #  limit (ext3 is 16TB, for instance).
            raise IOError("XDR seek failed with system errno={}".format(ok))
        self.current_frame = frame

    def _bytes_seek(self, offset, whence="SEEK_SET"):
        """Low-level access to the stream repositioning xdr_seek call.

        Beware that this function will not update :attr:`current_frame`,
        even on a successful seek, and might render the :class:`Reader`
        unstable.

        Parameters
        ----------
        offset : int
            move the stream *offset* bytes relative to *whence*
            (can be negative).

        whence : str, optional
            must be one of:
             "SEEK_SET": *offset* is applied relative to the beginning
                         of the file.
             "SEEK_CUR": *offset* is applied relative to the current
                         position in the stream.
             "SEEK_END": *offset* is applied relative to the end of
                         the file (a positive *offset* will seek beyond the
                         end of the file, without an error)

        Raises
        ------
        ValueError
            If *whence* is not one of the expected strings.

        IOError
            If the seek fails (the low-level system error is reported).
        """
        cdef int whn
        cdef int64_t offst

        try:
            whn = _whence_vals[whence]
        except KeyError:
            raise ValueError("Parameter 'whence' must be "
                             "one of {}".format(tuple(_whence_vals.keys())))
        offst = offset
        self.reached_eof = False
        ok = xdr_seek(self.xfp, offst, whn)
        if ok != EOK:
            # See the comments to seek() for hints on errno meaning.
            raise IOError("XDR seek failed with system errno={}".format(ok))

    @property
    def offsets(self):
        """get byte offsets from current xdr file

        See Also
        --------
        set_offsets
        """
        if not self._has_offsets:
            self._offsets = self.calc_offsets()
            self._has_offsets = True
        return self._offsets

    def set_offsets(self, offsets):
        """set frame offsets"""
        self._offsets = offsets
        self._has_offsets = True

    def tell(self):
        """Get current frame"""
        return self.current_frame

    def _bytes_tell(self):
        """Low-level call to xdr_tell to get current byte offset."""
        return xdr_tell(self.xfp)


TRRFrame = namedtuple('TRRFrame', 'x v f box step time lmbda hasx hasv hasf')


cdef class TRRFile(_XDRFile):
    """File-like wrapper for gromacs TRR files.

    This class can be similar to the normal file objects in python. The read()
    function will return a frame and all information in it instead of a single
    line. Additionally the context-manager protocoll is supported as well.

    Parameters
    ----------
    fname : str
        The filename to open.
    mode : ('r', 'w')
        The mode in which to open the file, either 'r' read or 'w' write

    Examples
    --------
    >>> from MDAnalysis.lib.formats.libmdaxdr import TRRFile
    >>> with TRRFile('foo.trr') as f:
    >>>     for frame in f:
    >>>         print(frame.x)
    """

    def _calc_natoms(self, fname):
        cdef int n_atoms
        return_code = read_trr_natoms(fname, &n_atoms)
        return return_code, n_atoms


    def calc_offsets(self):
        """read byte offsets from TRR file directly"""
        if not self.is_open:
            return np.array([])
        cdef int n_frames = 0
        cdef int est_nframes = 0
        cdef int64_t* offsets = NULL
        ok = read_trr_n_frames(self.fname, &n_frames, &est_nframes, &offsets)
        if ok != EOK:
            raise RuntimeError("TRR couldn't calculate offsets. "
                               "XDR error = {}".format(error_message[ok]))
        # the read_xtc_n_frames allocates memory for the offsets with an
        # overestimation. This number is saved in est_nframes and we need to
        # tell the new numpy array about the whole allocated memory to avoid
        # memory leaks.
        cdef np.ndarray dims = np.array([est_nframes], dtype=np.int64)
        # this handles freeing the allocated memory correctly.
        cdef np.ndarray nd_offsets = ptr_to_ndarray(<void*> offsets, dims, np.NPY_INT64)
        return nd_offsets[:n_frames]

    def read(self):
        """Read next frame in the TRR file

        Returns
        -------
        frame : libmdaxdr.TRRFrame
            namedtuple with frame information

        See Also
        --------
        TRRFrame
        XTCFile

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

        return_code = 1
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

        return_code = read_trr(self.xfp, self.n_atoms, <int*> &step,
                                      &time, &lmbda, <matrix>box.data,
                                      <rvec*>xyz.data,
                                      <rvec*>velocity.data,
                                      <rvec*>forces.data,
                                      <int*> &has_prop)
        # trr are a bit weird. Reading after the last frame always always
        # results in an integer error while reading. I tried it also with trr
        # produced by different codes (Gromacs, ...).
        if return_code != EOK and return_code != EENDOFFILE \
           and return_code != EINTEGER:
            raise IOError('TRR read error = {}'.format(
                error_message[return_code]))

        # In a trr the integer error seems to indicate that the file is ending.
        # There might be corrupted files where this is a legitimate error. But
        # then we just can't read it and stop there which is not too bad.
        if return_code == EENDOFFILE or return_code == EINTEGER:
            self.reached_eof = True
            raise StopIteration

        if return_code == EOK:
            self.current_frame += 1

        has_x = bool(has_prop & HASX)
        has_v = bool(has_prop & HASV)
        has_f = bool(has_prop & HASF)
        return TRRFrame(xyz, velocity, forces, box, step, time, lmbda,
                        has_x, has_v, has_f)

    def write(self, xyz, velocity, forces, box, int step, float time,
              float _lambda, int natoms):
        """write one frame into TRR file.

        Parameters
        ----------
        xyz : ndarray, shape=(n_atoms, 3)
            cartesion coordinates
        box : ndarray, shape=(3, 3)
            Box vectors for this frame
        step : int
            current step number, 1 indexed
        time : float
            current time
        _lambda : float
            current lambda value
        natoms : int
            number of atoms in frame

        Raises
        ------
        RuntimeError
            Couldn't write the file
        ValueError
            The arguments do not match with previous saved frames.
        """
        if self.mode != 'w' :
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

        return_code = write_trr(self.xfp, self.n_atoms, step, time,
                                       _lambda, <matrix> box_ptr,
                                       <rvec*> xyz_ptr,
                                       <rvec*> velocity_ptr,
                                       <rvec*> forces_ptr)
        if return_code != EOK:
            raise IOError('TRR write error = {}'.format(
                error_message[return_code]))

        self.current_frame += 1


XTCFrame = namedtuple('XTCFrame', 'x box step time prec')


cdef class XTCFile(_XDRFile):
    """File-like wrapper for gromacs XTC files.

    This class can be similar to the normal file objects in python. The read()
    function will return a frame and all information in it instead of a single
    line. Additionally the context-manager protocoll is supported as well.

    Parameters
    ----------
    fname : str
        The filename to open.
    mode : ('r', 'w')
        The mode in which to open the file, either 'r' read or 'w' write

    Examples
    --------
    >>> from MDAnalysis.lib.formats.libmdaxdr import XTCFile
    >>> with XTCFile('foo.trr') as f:
    >>>     for frame in f:
    >>>         print(frame.x)
    """
    cdef float precision

    def _calc_natoms(self, fname):
        cdef int n_atoms
        return_code = read_xtc_natoms(fname, &n_atoms)
        return return_code, n_atoms


    def calc_offsets(self):
        """Calculate offsets from XTC file directly"""
        if not self.is_open:
            return np.array([])
        cdef int n_frames = 0
        cdef int est_nframes = 0
        cdef int64_t* offsets = NULL
        ok = read_xtc_n_frames(self.fname, &n_frames, &est_nframes, &offsets)
        if ok != EOK:
            raise RuntimeError("XTC couldn't calculate offsets. "
                               "XDR error = {}".format(error_message[ok]))
        # the read_xtc_n_frames allocates memory for the offsets with an
        # overestimation. This number is saved in est_nframes and we need to
        # tell the new numpy array about the whole allocated memory to avoid
        # memory leaks.
        cdef np.ndarray dims = np.array([est_nframes], dtype=np.int64)
        # this handles freeing the allocated memory correctly.
        cdef np.ndarray nd_offsets = ptr_to_ndarray(<void*> offsets, dims, np.NPY_INT64)
        return nd_offsets[:n_frames]

    def read(self):
        """Read next frame in the XTC file

        Returns
        -------
        frame : libmdaxdr.XTCFrame
            namedtuple with frame information

        See Also
        --------
        XTCFrame
        TRRFile

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

        return_code = 1
        cdef int step
        cdef float time, prec

        cdef np.ndarray xyz = np.empty((self.n_atoms, DIMS), dtype=DTYPE)
        cdef np.ndarray box = np.empty((DIMS, DIMS), dtype=DTYPE)

        return_code = read_xtc(self.xfp, self.n_atoms, <int*> &step,
                                      &time, <matrix>box.data,
                                      <rvec*>xyz.data, <float*> &prec)
        if return_code != EOK and return_code != EENDOFFILE:
            raise IOError('XTC read error = {}'.format(
                error_message[return_code]))

        if return_code == EENDOFFILE:
            self.reached_eof = True
            raise StopIteration

        if return_code == EOK:
            self.current_frame += 1
        return XTCFrame(xyz, box, step, time, prec)

    def write(self, xyz, box, int step, float time, float precision=1000):
        """write one frame to the XTC file

        Parameters
        ----------
        xyz : ndarray, shape=(n_atoms, 3)
            cartesion coordinates
        box : ndarray, shape=(3, 3)
            Box vectors for this frame
        step : int
            current step number, 1 indexed
        time : float
            current time
        precision : float (optional)
            precision of saved trajectory, see Notes


        Notes
        -----
        Gromacs specifies the precision a little bit strange. The coordinates
        will be multiplied by precision and then converted to an integer. This
        means a precision of 1000 yields an accuracy of 3 decimal places.

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
            self.precision = precision
        else:
            if self.n_atoms != xyz.shape[0]:
                raise ValueError('Previous frames contained {} atoms. You '
                                 'are trying to write {} atoms.'.format(
                                     self.n_atoms, xyz.shape[1]))
            if self.precision != precision:
                raise ValueError('Previous frames used precision of {}. You '
                                 'are trying to use {}'.format(
                                     self.precision, precision))

        return_code = write_xtc(self.xfp, self.n_atoms, step, time,
                                       <matrix>&box_view[0, 0],
                                       <rvec*>&xyz_view[0, 0], precision)
        if return_code != EOK:
            raise IOError('XTC write error = {}'.format(
                error_message[return_code]))

        self.current_frame += 1
