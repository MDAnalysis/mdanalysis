# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
"""\
Low-level Gromacs XDR trajectory reading — :mod:`MDAnalysis.lib.formats.libmdaxdr`
----------------------------------------------------------------------------------

:mod:`libmdaxdr` contains the classes :class:`XTCFile` and
:class:`TRRFile`. Both can be used to read and write frames from and to
Gromacs_ XTC and TRR files. These classes are used internally by MDAnalysis in
:mod:`MDAnalysis.coordinates.XTC` and :mod:`MDAnalysis.coordinates.TRR`. They
behave similar to normal file objects.

For example, one can use a :class:`XTCFile` to directly calculate mean
coordinates (where the coordinates are stored in `x` attribute of the
:class:`namedtuple` `frame`):

.. code-block:: python
   :emphasize-lines: 1,2,5

   with XTCFile("trajectory.xtc") as xtc:
      n_atoms = xtc.n_atoms
      mean = np.zeros((n_atoms, 3))
      # iterate over trajectory
      for frame in xtc:
          mean += frame.x

The :class:`XTCFile` class can be useful as a compressed storage format.

Besides iteration, :class:`XTCFile` and :class:`TRRFile` one can also seek to
arbitrary frames using the :meth:`~XTCFile.seek` method. This is provided by
lazily generating a offset list for stored frames. The offset list is generated
the first time :func:`len` or :`~XTCFile.seek` is called.

(For more details on how to use :class:`XTCFile` and :class:`TRRFile` on their
own please see the source code in `lib/formats/libmdaxdr.pyx`_ for the time being.)


.. _Gromacs: http://www.gromacs.org
.. _`lib/formats/libmdaxdr.pyx`:
   https://github.com/MDAnalysis/mdanalysis/blob/develop/package/MDAnalysis/lib/formats/libmdaxdr.pyx
"""

cimport numpy as np
cimport cython
from MDAnalysis.lib.formats.cython_util cimport ptr_to_ndarray
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
    IOError

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
    cdef readonly int _has_offsets

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
        IOError
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
            raise IOError('mode must be one of "r" or "w", you '
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
        """Return next frame"""
        if self.reached_eof:
            raise StopIteration
        return self.read()

    def __len__(self):
        if not self.is_open:
            raise IOError('No file currently opened')
        return self.offsets.size

    def __reduce__(self):
        """low level pickle function. It gives instructions for unpickling which class
        should be created with arguments to `__cinit__` and the current state

        """
        return (self.__class__, (self.fname.decode(), self.mode),
                self.__getstate__())

    def __getstate__(self):
        """
        current state
        """
        return (self.is_open, self.current_frame, self._offsets,
                self._has_offsets)

    def __setstate__(self, state):
        """
        reset state offsets and current frame
        """
        is_open = state[0]
        # by default the file is open
        if not is_open:
            self.close()
            return

        # set them now to avoid recalculation
        self._offsets = state[2]
        self._has_offsets = state[3]

        # where was I
        current_frame = state[1]
        if current_frame < self.offsets.size:
            self.seek(current_frame)
        elif current_frame == self.offsets.size:
            #  cannot seek to self.offsets.size (a.k.a len(_XDRFile));
            #  instead, we seek to the previous frame and read next. 
            #  which is the state of the file when we need to serialize
            #  at the end of the trajectory.
            self.seek(current_frame - 1)
            _ = self.read()
        else:       # pragma: no cover
            raise RuntimeError("Invalid frame number {} > {} -- this should"
                               "not happen.".format(current_frame,
                                                    self.offsets.size)
                              )

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
        """
        cdef int64_t offset
        if frame == 0:
            offset = 0
        elif frame < 0:
            raise IOError("Can't seek to negative frame")
        else:
            if frame >= self.offsets.size:
                raise EOFError('Trying to seek over max number of frames')
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
        IOError
        """
        cdef int whn
        cdef int64_t offst

        try:
            whn = _whence_vals[whence]
        except KeyError:
            raise IOError("Parameter 'whence' must be "
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

    Notes
    -----
    This class can be pickled. The pickle will store filename, mode, current
    frame and offsets
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
            raise IOError("TRR couldn't calculate offsets. "
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
        IOError
        """
        if self.reached_eof:
            raise EOFError('Reached last frame in TRR, seek to 0')
        if not self.is_open:
            raise IOError('No file opened')
        if self.mode != 'r':
            raise IOError('File opened in mode: {}. Reading only allow '
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
        xyz : array_like, shape=(`n_atoms`, 3), optional
            cartesion coordinates. Only written if not ``None``.
        velocity : array_like, shape=(`n_atoms`, 3), optional
            cartesion velocities. Only written if not ``None``.
        forces : array_like, shape=(`n_atoms`, 3), optional
            cartesion forces. Only written if not ``None``.
        box : array_like, shape=(3, 3)
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
        IOError
        """
        if self.mode != 'w' :
            raise IOError('File opened in mode: {}. Writing only allow '
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
            xyz = np.asarray(xyz)
            xyz_helper = np.ascontiguousarray(xyz, dtype=DTYPE)
            xyz_ptr = <float*>xyz_helper.data
        if velocity is not None:
            velocity = np.asarray(velocity)
            velocity_helper = np.ascontiguousarray(velocity, dtype=DTYPE)
            velocity_ptr = <float*>velocity_helper.data
        if forces is not None:
            forces = np.asarray(forces)
            forces_helper = np.ascontiguousarray(forces, dtype=DTYPE)
            forces_ptr = <float*>forces_helper.data

        box = np.asarray(box)
        cdef np.ndarray box_helper = np.ascontiguousarray(box, dtype=DTYPE)
        cdef float* box_ptr = <float*>box_helper.data

        if self.current_frame == 0:
            self.box = box
            self.n_atoms = natoms
        else:
            if self.n_atoms != natoms:
                raise IOError('Previous frames contained {} atoms. You '
                              'are trying to write {} atoms.'.format(
                                  self.n_atoms, natoms))
            if xyz is not None and self.n_atoms != xyz.shape[0]:
                raise IOError('Previous frames xyz contained {} atoms. You '
                              'are trying to write {} atoms.'.format(
                                  self.n_atoms, xyz.shape[0]))
            if velocity is not None and self.n_atoms != velocity.shape[0]:
                raise IOError('Previous frames velocity contained {} atoms. You '
                              'are trying to write {} atoms.'.format(
                                  self.n_atoms, velocity.shape[0]))
            if forces is not None and self.n_atoms != forces.shape[0]:
                raise IOError('Previous frames forces contained {} atoms. You '
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

    Notes
    -----
    This class can be pickled. The pickle will store filename, mode, current
    frame and offsets
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
            raise IOError("XTC couldn't calculate offsets. "
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
        IOError
        """
        if self.reached_eof:
            raise EOFError('Reached last frame in XTC, seek to 0')
        if not self.is_open:
            raise IOError('No file opened')
        if self.mode != 'r':
            raise IOError('File opened in mode: {}. Reading only allow '
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
        xyz : array_like, shape=(`n_atoms`, 3)
            cartesion coordinates
        box : array_like, shape=(3, 3)
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
        IOError

        """
        if self.mode != 'w':
            raise IOError('File opened in mode: {}. Writing only allow '
                          'in mode "w"'.format('self.mode'))

        xyz = np.asarray(xyz)
        box = np.asarray(box)

        cdef DTYPE_T[:, ::1] xyz_view = np.ascontiguousarray(xyz, dtype=DTYPE)
        cdef DTYPE_T[:, ::1] box_view = np.ascontiguousarray(box, dtype=DTYPE)

        if self.current_frame == 0:
            self.n_atoms = xyz.shape[0]
            self.box = box
            self.precision = precision
        else:
            if self.n_atoms != xyz.shape[0]:
                raise IOError('Previous frames contained {} atoms. You '
                              'are trying to write {} atoms.'.format(
                                  self.n_atoms, xyz.shape[1]))
            if self.precision != precision:
                raise IOError('Previous frames used precision of {}. You '
                              'are trying to use {}'.format(
                                  self.precision, precision))

        return_code = write_xtc(self.xfp, self.n_atoms, step, time,
                                       <matrix>&box_view[0, 0],
                                       <rvec*>&xyz_view[0, 0], precision)
        if return_code != EOK:
            raise IOError('XTC write error = {}'.format(
                error_message[return_code]))

        self.current_frame += 1
