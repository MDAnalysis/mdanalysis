/* -*- Mode: C; indent-tabs-mode:nil -*-  */
/* SWIG interface for Gromacs libxdrfile2 with the XTC and TRR code 
   Copyright (c) 2010 Oliver Beckstein <orbeckst@gmail.com>
   Copyright (c) 2013,2014 Manuel Melo <manuel.nuno.melo@gmail.com>
   Published under the GNU GENERAL PUBLIC LICENSE Version 2 (or higher)

   swig -python -outdir MDAnalysis/coordinates/xdrfile src/xdrfile/libxdrfile2.i
*/
%define DOCSTRING
"
:Author:  Oliver Beckstein <orbeckst@gmail.com>
:Author:  Manuel Melo <manuel.nuno.melo@gmail.com>
:Year:    2014
:Licence: GNU GENERAL PUBLIC LICENSE Version 2 (or higher)


The Gromacs XTC/TRR library :mod:`libxdrfile2`
==============================================

:mod:`libxdrfile2`, a derivative of the Gromacs_ `libxdrfile library`_,
provides an interface to some high-level functions for XTC/TRR trajectory
handling.  Only functions required for reading and processing whole
trajectories are exposed at the moment; low-level routines to read individual
numbers are not provided. In addition, :mod:`libxdrfile2` exposes functions to
allow fast frame indexing and XDR file seeking.

The functions querying the numbers of atoms in a trajectory frame
(:func:`read_xtc_natoms` and :func:`read_trr_natoms`) open a file themselves
and only require the file name.

All other functions operate on a *XDRFILE* object, which is a special file
handle for xdr files.  Any xdr-based trajectory file (XTC or TRR format) always
has to be opened with :func:`xdrfile_open`. When done, close the trajectory
with :func:`xdrfile_close`.

The functions fill or read existing arrays of coordinates; they never allocate
these arrays themselves. Hence they need to be setup outside libxdrfile2 as
numpy arrays. The exception to these are the indexing ones functions that take
care of array allocation and transference to a garbage-collectable memory
object.


.. _Gromacs: http://www.gromacs.org
.. _libxdrfile library: 
   http://www.gromacs.org/Developer_Zone/Programming_Guide/XTC_Library

.. versionchanged:: 0.8.0
   :mod:`libxdrfile2` is now used instead of
   :mod:`libxdrfile`. :mod:`libxdrfile2` is based on :mod:`libxdrfile` but has
   xdr seeking and indexing capabilities.  Unlike :mod:`libxdrfile` before it,
   :mod:`libxdrfile2` is distributed under the GNU GENERAL PUBLIC LICENSE,
   version 2 (or higher).


Example: Reading from a XTC
---------------------------

In the example we read coordinate frames from an existing XTC trajectory::

  import numpy as np
  from libxdrfile2 import xdrfile_open, xdrfile_close, read_xtc_natoms, read_xtc, DIM, exdrOK
  xtc = 'md.xtc'
  
  # get number of atoms
  natoms = read_xtc_natoms(xtc)

  # allocate coordinate array of the right size and type
  # (the type float32 is crucial to match the underlying C-code!!)
  x = np.zeros((natoms, DIM), dtype=np.float32)
  # allocate unit cell box
  box = np.zeros((DIM, DIM), dtype=np.float32)

  # open file
  XTC = xdrfile_open(xtc, 'r')

  # loop through file until return status signifies end or a problem
  # (it should become exdrENDOFFILE on the last iteration)
  status = exdrOK
  while status == exdrOK:
     status,step,time,prec = read_xtc(XTC, box, x)
     # do something with x
     centre = x.mean(axis=0)
     print 'Centre of geometry at %(time)g ps: %(centre)r' % vars()

  # finally close file
  xdrfile_close(XTC)

Note that only the *contents* of the coordinate and unitcell arrays *x* and
*box* change.


Functions and constants
-----------------------

The module defines a number of constants such as :data:`DIM` or the
`Status symbols`_.

.. data:: DIM

          The number of cartesian dimensions for which the underlying C-code
          was compiled; this is most certainly 3.


Status symbols
~~~~~~~~~~~~~~

A number of symbols are exported; they all start with the letters
``exdr``. Important ones are:

.. data:: exdrOK

          Success of xdr file read/write operation.

.. data:: exdrCLOSE
 
          xdr file is closed

.. data:: exdrENDOFFILE

          end of file was reached (response of :func:`read_xtc` and
          :func:`read_trr` after the last read frame)

.. data:: exdrFILENOTFOUND

          :func:`xdrfile_open` cannot find the requested file

Other symbols that are used internally are:

.. data:: exdrHEADER

          header

.. data:: exdrSTRING

          string

.. data:: exdrDOUBLE

          double precision floating point number

.. data:: exdrINT

          integer

.. data:: exdrFLOAT

          floating point number

.. data:: exdrUINT

          unsigned integer

.. data:: exdr3DX

          compressed 3D coordinates

.. data:: exdrMAGIC

          magic number

.. data:: exdrNOMEM

          not enough memory to allocate space for a XDR data structure.      

Opening and closing of XDR files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Two low-level functions are used to obtain a *XDRFILE* object (a file handle)
to access xdr files such as XTC or TRR trajectories.

.. function:: xdrfile_open(path, mode) -> XDRFILE

              Open *path* and returns a *XDRFILE* handle that is required by other
              functions.

              :Arguments:
                  *path*
                     file name
                  *mode*
                     'r' for reading and 'w' for writing
              :Returns: *XDRFILE* handle

.. function:: xdrfile_close(XDRFILE) -> status

              Close the xdrfile pointed to by *XDRFILE*. 

              .. Warning:: Closing an already closed file will lead to a 
                           crash with a double-free pointer error.

XTC functions
~~~~~~~~~~~~~

The XTC trajectory format is a lossy compression format that only stores
coordinates. Compression level is determined by the *precision* argument to the
:func:`write_xtc` function. Coordinates (Gromacs_ uses nm natively) are
multiplied by *precision* and truncated to the integer part. A typical value is
1000.0, which gives an accuracy of 1/100 of an Angstroem.

The advantage of XTC over TRR is its significantly reduced size.


.. function:: read_xtc_natoms(fn) -> natoms

              Read the number of atoms *natoms* from a xtc file *fn*.

              :Arguments:
                *fn*
                   file name of an xtc file

              :Raises: :exc:`IOError` if the supplied filed is not a XTC 
                       or if it is not readable.

.. function:: read_xtc_n_frames(fn) -> (n_frames, offsets)

              Read through the whole trajectory headers to obtain the 
              total number of frames. The process is speeded up by reading 
              frame headers for the amount of data in the frame, 
              and then skipping directly to the next header. An array of 
              frame offsets is also returned, which can later be used to 
              seek direcly to arbitrary frames in the trajectory. 

              :Arguments:
                *fn*
                   file name of an xtc file

              :Returns:
                a tuple containing:
                  *n_frames*
                     an int with the total frame count in the trajectory
                  *offsets*
                     a numpy array of int64 recording the starting byte offset of each frame

              :Raises: :exc:`IOError` if the supplied filed is not a XTC 
                       or if it is not readable.

.. function:: read_xtc(XDRFILE, box, x) -> (status, step, time, precision)

              Read the next frame from the opened xtc trajectory into *x*.

              :Arguments:
                *XDRFILE*
                   open *XDRFILE* object
                *box*
                   pre-allocated numpy ``array((DIM,DIM),dtype=numpy.float32)`` which
                   is filled with the unit cell box vectors
                *x*
                   pre-allocated numpy ``array((natoms, DIM),dtype=numpy.float32)``
                   which is updated with the coordinates from the frame

              :Returns:
                a tuple containing:
                  *status*
                     integer status (0 = exdrOK), see `Status symbols`_ for other
                     values)
                  *step*
                     simulation step
                  *time*
                     simulation time in ps
                  *precision*
                     precision of the lossy xtc format (typically 1000.0)

.. function:: write_xtc(XDRFILE, step, time, box, x, prec) -> status

              Write the next frame *x* to the opened xtc trajectory.

              :Arguments:
                *XDRFILE*
                   open *XDRFILE* object (writable)
                *step*
                   simulation step
                *time*
                   time step in ps
                *box*
                   numpy ``array((DIM,DIM),dtype=numpy.float32)`` which contains 
                   the unit cell box vectors
                *x*
                   numpy ``array((natoms, DIM),dtype=nump.float32)``
                   which contains the coordinates from the frame
                *precision*
                   precision of the lossy xtc format (typically 1000.0)

              :Returns: *status*, integer status (0 = OK), see the ``libxdrfile2.exdr*`` 
                        constants under `Status symbols`_ for other values)

TRR functions
~~~~~~~~~~~~~

TRR is the Gromacs_ native full-feature trajectory storage format. It can contain position 
coordinates, velocities and forces, and the lambda value for free energy perturbation 
calculations. Velocities and forces are optional in the sense that they can be all zero.

.. function:: read_trr_natoms(fn) -> natoms

              Read the number of atoms *natoms* from a trr file *fn*.

              :Arguments:
                *fn*
                   file name of a trr file

              :Raises: :exc:`IOError` if the supplied filed is not a TRR
                       or if it is not readable.

.. function:: read_trr_n_frames(fn) -> (n_frames, offsets)

              Read through the whole trajectory headers to obtain the total number of frames. 
              The process is speeded up by reading frame headers for the amount of data in the frame,
              and then skipping directly to the next header. An array of frame offsets is also
              returned, which can later be used to seek direcly to arbitrary frames in the trajectory. 

              :Arguments:
                *fn*
                   file name of an xtc file

              :Returns:
                a tuple containing:
                  *n_frames*
                     an int with the total frame count in the trajectory
                  *offsets*
                     a numpy array of int64 recording the starting byte offset of each frame

              :Raises: :exc:`IOError` if the supplied filed is not a TRR or if it is not readable.

.. function:: read_trr(XDRFILE, box, x, v, f) -> (status, step, time, lambda)

              Read the next frame from the opened trr trajectory into *x*, *v*, and *f*.

              :Arguments:
                *XDRFILE*
                   open *XDRFILE* object
                *box*
                   pre-allocated numpy ``array((DIM,DIM),dtype=numpy.float32)`` which
                   is filled with the unit cell box vectors
                *x*
                   pre-allocated numpy ``array((natoms, DIM),dtype=nump.float32)``
                   which is updated with the **coordinates** from the frame
                *v*
                   pre-allocated numpy ``array((natoms, DIM),dtype=nump.float32)``
                   which is updated with the **velocities** from the frame
                *f*
                   pre-allocated numpy ``array((natoms, DIM),dtype=nump.float32)``
                   which is updated with the **forces** from the frame

              :Returns:
                a tuple containing:
                  *status*
                     integer status (0 = exdrOK), see the ``libxdrfile2.exdr*`` constants 
                     under `Status symbols`_ for other values)
                  *step*
                     simulation step
                  *time*
                     simulation time in ps
                  *lambda*
                     current lambda value (only interesting for free energy perturbation)
                  *has_x*
                     boolean indicating whether coordinates were read from the TRR
                  *has_v*
                     boolean indicating whether velocities were read from the TRR
                  *has_f*
                     boolean indicating whether forces were read from the TRR

.. function:: write_trr(XDRFILE, step, time, lambda, box, x, v, f) -> status

              Write the next frame to the opened trr trajectory.

              :Arguments:
                *XDRFILE*
                   open *XDRFILE* object (writable)
                *step*
                   simulation step
                *time*
                   time step in ps
                *lambda*
                   free energy lambda value (typically 0.0)
                *box*
                   numpy ``array((DIM,DIM),dtype=numpy.float32)`` which contains 
                   the unit cell box vectors
                *x*
                   numpy ``array((natoms, DIM),dtype=nump.float32)``
                   which contains the **coordinates** from the frame
                *v*
                   numpy ``array((natoms, DIM),dtype=nump.float32)``
                   which contains the **velocities** from the frame
                *f*
                   numpy ``array((natoms, DIM),dtype=nump.float32)``
                   which contains the **forces** from the frame

              .. versionchanged:: 0.8.0
                   either one of *x*, *v*, or *f* can now be set as a natom,0-DIM
                   numpy ``array((natom, 0),dtype=nump.float32)``. This will cause the
                   corresponding property to be skipped when writing to file.
 
              :Returns: *status*, integer status (0 = OK), see the ``libxdrfile2.exdr*`` 
                        constants under `Status symbols`_ for other values)

"
%enddef

%module(docstring=DOCSTRING) libxdrfile2


%{
/* Python SWIG interface to some functions in Gromacs libxdr v 2.0
   Copyright (c) 2010 Oliver Beckstein <orbeckst@gmail.com>
   Published under the GNU LESSER GENERAL PUBLIC LICENSE Version 3 (or higher)
   See http://www.mdanalysis.org for details.
 */
#define SWIG_FILE_WITH_INIT
#include <stdio.h>
#include "xdrfile.h"
#include "xdrfile_trr.h"
#include "xdrfile_xtc.h"
%}

%include "numpy.i"

%init %{
import_array();
%}


/* From Gromacs xdrfile.c 

   I am only wrapping 'high level' functions and modify call
   signatures so that one does not need anything like pointers from
   python.
*/


/* status codes */
enum { exdrOK, exdrHEADER, exdrSTRING, exdrDOUBLE, 
       exdrINT, exdrFLOAT, exdrUINT, exdr3DX, exdrCLOSE, exdrMAGIC,
       exdrNOMEM, exdrENDOFFILE, exdrFILENOTFOUND, exdrNR };

/* These com from stdio.h, for file seeking. Gives all the flexibility to _fseek(). */
enum { SEEK_SET, SEEK_CUR, SEEK_END };

/* open/close xdr files */
%feature("autodoc", "0") xdrfile_open;
extern XDRFILE* xdrfile_open(const char *path, const char *mode);

%feature("autodoc", "0") xdrfile_close;
extern int xdrfile_close(XDRFILE *fp);


/* from xdrfile_xtc.c */
/* This function returns the number of atoms in the xtc file in *natoms
     extern int read_xtc_natoms(char *fn,int *natoms);
   ... but the wrapped function returns natoms as the python return value
*/
%feature("autodoc", "0") my_read_xtc_natoms;
%rename (read_xtc_natoms) my_read_xtc_natoms;
%exception my_read_xtc_natoms {
  $action
  if (PyErr_Occurred()) SWIG_fail;
}
%inline %{
  int my_read_xtc_natoms(char *fn) {
    int natoms;
    int status;
    status = read_xtc_natoms(fn, &natoms);
    if (status != exdrOK) {
      PyErr_Format(PyExc_IOError, "[%d] Error reading natoms from xtc '%s'", status, fn);
      return 0;
    }
    return natoms;
  }
%}

%feature("autodoc", "0") my_read_xtc_n_frames;
%rename (read_xtc_n_frames) my_read_xtc_n_frames;
%exception my_read_xtc_n_frames {
  $action
  if (PyErr_Occurred()) SWIG_fail;
}
%inline %{
PyObject * my_read_xtc_n_frames(char *fn) {
    int n_frames, status;
    int64_t *offsets[1];
    PyObject *npoffsets = NULL;
    status = read_xtc_n_frames(fn, &n_frames, offsets);
    if (status != exdrOK) {
      PyErr_Format(PyExc_IOError, "[%d] Error reading n_frames by seeking through xtc '%s'", status, fn);
      return 0;
    }
    npy_intp nfrms[1] = { n_frames };
    npoffsets = PyArray_SimpleNewFromData(1, nfrms, NPY_INT64, *offsets);
    if (npoffsets==NULL)
    {
      free(*offsets);
      Py_XDECREF(npoffsets);
      PyErr_Format(PyExc_IOError, "Error copying frame index into Python.");
      return 0;
    }
    /* From http://web.archive.org/web/20130304224839/http://blog.enthought.com/python/numpy/simplified-creation-of-numpy-arrays-from-pre-allocated-memory/ */
    PyArray_BASE(npoffsets) = PyCObject_FromVoidPtr(*offsets, free);
    PyObject *tuple = PyTuple_New(2); 
    PyTuple_SET_ITEM(tuple, 0, PyInt_FromLong((long)n_frames));
    PyTuple_SET_ITEM(tuple, 1, npoffsets);
    return tuple;
  }
%}


/* This function returns the number of atoms in the trr file in *natoms 
     extern int read_trr_natoms(char *fn,int *natoms);
 ... but the wrapped function returns natoms as the python return value 
*/
%feature("autodoc", "0") my_read_trr_natoms;
%rename (read_trr_natoms) my_read_trr_natoms;
%exception my_read_trr_natoms {
  $action
  if (PyErr_Occurred()) SWIG_fail;
}
%inline %{
  int my_read_trr_natoms(char *fn) {
    int natoms;
    int status;
    status = read_trr_natoms(fn, &natoms);
    if (status != exdrOK) {
      PyErr_Format(PyExc_IOError, "[%d] Error reading natoms from trr '%s'", status, fn);
      return 0;
    }
    return natoms;
  }
%}


%feature("autodoc", "0") my_read_trr_n_frames;
%rename (read_trr_n_frames) my_read_trr_n_frames;
%exception my_read_trr_n_frames {
  $action
  if (PyErr_Occurred()) SWIG_fail;
}
%inline %{
PyObject * my_read_trr_n_frames(char *fn) {
    int n_frames, status;
    int64_t *offsets[1];
    PyObject *npoffsets = NULL;
    status = read_trr_n_frames(fn, &n_frames, offsets);
    if (status != exdrOK) {
      PyErr_Format(PyExc_IOError, "[%d] Error reading n_frames by seeking through trr '%s'", status, fn);
      return 0;
    }
    npy_intp nfrms[1] = { n_frames };
    npoffsets = PyArray_SimpleNewFromData(1, nfrms, NPY_INT64, *offsets);
    if (npoffsets==NULL)
    {
      free(*offsets);
      Py_XDECREF(npoffsets);
      PyErr_Format(PyExc_IOError, "Error copying frame index into Python.");
      return 0;
    }
    /* From http://web.archive.org/web/20130304224839/http://blog.enthought.com/python/numpy/simplified-creation-of-numpy-arrays-from-pre-allocated-memory/ */
    PyArray_BASE(npoffsets) = PyCObject_FromVoidPtr(*offsets, free);
    PyObject *tuple = PyTuple_New(2); 
    PyTuple_SET_ITEM(tuple, 0, PyInt_FromLong((long)n_frames));
    PyTuple_SET_ITEM(tuple, 1, npoffsets);
    return tuple;
  }
%}

  

#define DIM 3
typedef float matrix[DIM][DIM];
typedef float rvec[DIM];


/* Reading from xdr files */

%apply (float INPLACE_ARRAY2[ANY][ANY]) {(matrix box)}
%apply (int DIM1, int DIM2, float* INPLACE_ARRAY2) {(int natoms,  int _DIM,  float *x),
                                                    (int vnatoms, int v_DIM, float *v),
                                                    (int fnatoms, int f_DIM, float *f)}

/* Read one frame of an open xtc file */
/*
extern int read_xtc(XDRFILE *xd,int natoms,int *step,float *time,
                    matrix box,rvec *x,float *prec); 
*/
%feature("autodoc", "read_xtc(XDRFILE, box, x) -> (status, step, time, precision)") my_read_xtc;
%rename (read_xtc) my_read_xtc;
%inline %{
PyObject * my_read_xtc(XDRFILE *xd, matrix box, int natoms, int _DIM, float *x) {
  /* _DIM = 3 always, need to reorder for numpy.i SWIG */
  int status, step;
  float time, prec;
  PyObject *tuple = PyTuple_New(4); 
  status = read_xtc(xd, natoms, &step, &time, box, (rvec *)x, &prec);
  PyTuple_SET_ITEM(tuple, 0, PyInt_FromLong((long)status));
  PyTuple_SET_ITEM(tuple, 1, PyInt_FromLong((long)step));
  PyTuple_SET_ITEM(tuple, 2, PyFloat_FromDouble((double)time));
  PyTuple_SET_ITEM(tuple, 3, PyFloat_FromDouble((double)prec));
  return tuple; // return  (status, step, time, prec)
}
%}

%feature("autodoc", "read_trr(XDRFILE, box, x, v, f) -> (status, step, time, lambda)") my_read_trr;
%rename (read_trr) my_read_trr;
%inline %{
PyObject * my_read_trr(XDRFILE *xd, matrix box, 
                int natoms,  int _DIM,  float *x,
                int vnatoms, int v_DIM, float *v,
                int fnatoms, int f_DIM, float *f) {
  /* _DIM = 3 always, need to reorder for numpy.i SWIG */
  int status, step, has_prop=0;
  float time, lmbda;
  PyObject *tuple = PyTuple_New(7); 
  status = read_trr(xd, natoms, &step, &time, &lmbda, box, (rvec *)x, (rvec *)v, (rvec *)f, &has_prop);
  PyTuple_SET_ITEM(tuple, 0, PyInt_FromLong((long)status));
  PyTuple_SET_ITEM(tuple, 1, PyInt_FromLong((long)step));
  PyTuple_SET_ITEM(tuple, 2, PyFloat_FromDouble((double)time));
  PyTuple_SET_ITEM(tuple, 3, PyFloat_FromDouble((double)lmbda));
  PyTuple_SET_ITEM(tuple, 4, PyBool_FromLong((long)(has_prop & HASX)));
  PyTuple_SET_ITEM(tuple, 5, PyBool_FromLong((long)(has_prop & HASV)));
  PyTuple_SET_ITEM(tuple, 6, PyBool_FromLong((long)(has_prop & HASF)));
  return tuple; // return  (status, step, time, lmbda, has_x, has_v, has_f)
}
%}

%clear (matrix box);
%clear (int natoms,  int _DIM,  float *x);
%clear (int vnatoms, int v_DIM, float *v);
%clear (int fnatoms, int f_DIM, float *f);


/* Writing of xdr files */

%apply (float IN_ARRAY2[ANY][ANY]) {(matrix box)}
%apply (int DIM1, int DIM2, float* IN_ARRAY2) {(int natoms,  int _DIM,  float *x),
                                               (int vnatoms, int v_DIM, float *v),
                                               (int fnatoms, int f_DIM, float *f)}
  
/* Write a frame to xtc file */
/*
extern int write_xtc(XDRFILE *xd, int natoms,int step,float time,
                     matrix box,rvec *x,float prec);
*/
%feature("autodoc", "write_xtc(XDRFILE, step, time, box, x, prec) -> status") my_write_xtc;
%rename (write_xtc) my_write_xtc;
%inline %{
int my_write_xtc(XDRFILE *xd, int step, float time,
                 matrix box, int natoms, int _DIM, float *x, float prec) {
  /* _DIM = 3 always, need to reorder for numpy.i SWIG */
  return write_xtc(xd, natoms, step, time, box, (rvec *)x, prec);
}
%}

%feature("autodoc", "write_trr(XDRFILE, step, time, lambda, box, x, v, f) -> status") my_write_trr;
%rename (write_trr) my_write_trr;
%inline %{
int my_write_trr(XDRFILE *xd, int step, float time, float lmbda, matrix box, 
                 int natoms,  int _DIM,  float *x, 
                 int vnatoms, int v_DIM, float *v, 
                 int fnatoms, int f_DIM, float *f) { 
  /* Preparing for the case of empty arrays - NULL pointers tell the library to skip this property. */
  if (_DIM == 0) x = NULL;
  if (v_DIM == 0) v = NULL;
  if (f_DIM == 0) f = NULL;
  return write_trr(xd, natoms, step, time, lmbda, box, (rvec *)x, (rvec *)v, (rvec *)f);
}
%}

%feature("autodoc", "0") xdr_seek;
extern int xdr_seek(XDRFILE *xd, long long pos, int whence);

%feature("autodoc", "0") xdr_tell;
extern long long xdr_tell(XDRFILE *xd);

%clear (matrix box);
%clear (int natoms,  int _DIM,  float *x);
%clear (int vnatoms, int v_DIM, float *v);
%clear (int fnatoms, int f_DIM, float *f);


