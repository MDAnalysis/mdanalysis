/* -*- C -*-  (not really, but good for syntax highlighting) */
/* SWIG interface for Gromacs libxdrfile 1.1 with the xtc and trr code 
   Copyright (c) 2010 Oliver Beckstein <orbeckst@gmail.com>
   Published under the GNU LESSER GENERAL PUBLIC LICENSE Version 3 (or higher)

   swig -python -outdir MDAnalysis/coordinates/xdrfile src/xdrfile/libxdrfile.i
*/
%define DOCSTRING
"
:Author:  Oliver Beckstein <orbeckst@gmail.com>
:Year:    2010
:Licence: GNU LESSER GENERAL PUBLIC LICENSE Version 3 (or higher)


The Gromacs xtc/trr library :mod:`libxdrfile`
=============================================

:mod:`libxdrfile` provides an interface to some high-level functions in the
Gromacs_ `XTC Library`_ version 1.1.1. Only functions required for reading and
processing whole trajectories are exposed at the moment; low-level routines to
read individual numbers are not provided.

The functions querying the numbers of atoms in a trajectory frame
(:func:`read_xtc_natoms` and :func:`read_trr_natoms`) open a file themselves and
only require the file name.

All other functions operate on a *XDRFILE* object, which is a special file
handle for xdr files.  Any xdr-based trajectory file (xtc or trr format) always
has to be opened with :func:`xdrfile_open`. When done, close the trajectory
with :func:`xdrfile_close`.

The functions fill or read existing arrays of coordinates; they never allocate
these arrays themselves. Hence they need to be setup outside libxdrfile as
numpy arrays.


.. _Gromacs: http://www.gromacs.org
.. _XTC Library: http://www.gromacs.org/Developer_Zone/Programming_Guide/XTC_Library


Example: Reading from a xtc
---------------------------

In the example we read coordinate frames from an existing xtc trajectory::

  import numpy as np
  from libxdrfile import xdrfile_open, xdrfile_close, read_xtc_natoms, read_xtc, DIM, exdrOK
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
``exdr``. Important ones are listed here:

.. data:: exdrOK

          Success of xdr file read/write operation.

.. data:: exdrCLOSE
 
          xdr file is closed

.. data:: exdrENDOFFILE

          end of file was reached (response of :func:`read_xtc` and
          :func:`read_trr` after the last read frame)

.. data:: exdrFILENOTFOUND

          :func:`xdrfile_open` cannot find the requested file


Opening and closing of XDR files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Two low-level functions are used to obtain a *XDRFILE* object (a file handle)
to access xdr files such as xtc or trr trajectories.

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

.. function:: read_xtc_numframes(fn) -> numframes

              Read through the whole trajectory (!) to obtaine the total number of frames. 
              This can take a long time but it might still be advantageous to obtain 
              *numframes* in this way before setting up a complicated analysis. Unlike the DCD
              format, there is no way to obtain the total number of frames in the trajectory 
              except iterating through the whole trajectory.

              :Arguments:
                *fn*
                   file name of an xtc file

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

              :Returns: The function returns a tuple containing
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

              :Returns: *status*, integer status (0 = OK), see the ``libxdrfile.exdr*`` 
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

.. function:: read_trr_numframes(fn) -> numframes

              Read through the whole trajectory (!) to obtaine the total number of frames. 
              This can take a long time but it might still be advantageous to obtain 
              *numframes* in this way before setting up a complicated analysis. (This is a 
              poor implementation that loops through the *whole* trajectory and counts the 
              frames---please supply a better one.)

              :Arguments:
                *fn*
                   file name of an xtc file

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

              :Returns: The function returns a tuple containing
                *status*
                   integer status (0 = exdrOK), see the ``libxdrfile.exdr*`` constants 
                   under `Status symbols`_ for other values)
                *step*
                   simulation step
                *time*
                   simulation time in ps
                *lambda*
                   current lambda value (only interesting for free energy perturbation)

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
 
              :Returns: *status*, integer status (0 = OK), see the ``libxdrfile.exdr*`` 
                        constants under `Status symbols`_ for other values)

"
%enddef

%module(docstring=DOCSTRING) libxdrfile


%{
/* Python SWIG interface to some functions in Gromacs libxdr v 1.1
   Copyright (c) 2010 Oliver Beckstein <orbeckst@gmail.com>
   Published under the GNU LESSER GENERAL PUBLIC LICENSE Version 3 (or higher)
   See http://mdanalysis.googlecode.com for details.
 */
#define SWIG_FILE_WITH_INIT
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

%feature("autodoc", "0") my_read_xtc_numframes;
%rename (read_xtc_numframes) my_read_xtc_numframes;
%exception my_read_xtc_numframes {
  $action
  if (PyErr_Occurred()) SWIG_fail;
}
%inline %{
  int my_read_xtc_numframes(char *fn) {
    int numframes;
    int status;
    status = read_xtc_numframes(fn, &numframes);
    if (status != exdrOK) {
      PyErr_Format(PyExc_IOError, "[%d] Error reading numframes by iterating through xtc '%s'", status, fn);
      return 0;
    }
    return numframes;
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

%feature("autodoc", "0") my_read_trr_numframes;
%rename (read_trr_numframes) my_read_trr_numframes;
%exception my_read_trr_numframes {
  $action
  if (PyErr_Occurred()) SWIG_fail;
}
%inline %{
  int my_read_trr_numframes(char *fn) {
    int numframes;
    int status;
    status = read_trr_numframes(fn, &numframes);
    if (status != exdrOK) {
      PyErr_Format(PyExc_IOError, "[%d] Error reading numframes from trr '%s'", status, fn);
      return 0;
    }
    return numframes;
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
  int status, step;
  float time, lmbda;
  PyObject *tuple = PyTuple_New(4); 
  status = read_trr(xd, natoms, &step, &time, &lmbda, box, (rvec *)x, (rvec *)v, (rvec *)f);
  PyTuple_SET_ITEM(tuple, 0, PyInt_FromLong((long)status));
  PyTuple_SET_ITEM(tuple, 1, PyInt_FromLong((long)step));
  PyTuple_SET_ITEM(tuple, 2, PyFloat_FromDouble((double)time));
  PyTuple_SET_ITEM(tuple, 3, PyFloat_FromDouble((double)lmbda));
  return tuple; // return  (status, step, time, lmbda)
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

%feature("autodoc", "write_xtc(XDRFILE, step, time, lambda, box, x, v, f) -> status") my_write_trr;
%rename (write_trr) my_write_trr;
%inline %{
int my_write_trr(XDRFILE *xd, int step, float time, float lmbda, matrix box, 
		 int natoms,  int _DIM,  float *x, 
		 int vnatoms, int v_DIM, float *v, 
		 int fnatoms, int f_DIM, float *f) { 
  return write_trr(xd, natoms, step, time, lmbda, box, (rvec *)x, (rvec *)v, (rvec *)f);
}
%}

%clear (matrix box);
%clear (int natoms,  int _DIM,  float *x);
%clear (int vnatoms, int v_DIM, float *v);
%clear (int fnatoms, int f_DIM, float *f);


