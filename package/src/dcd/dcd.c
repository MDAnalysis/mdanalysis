/* -*- Mode: C; tab-width: 4; indent-tabs-mode:nil; -*- */
/* vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 */
/*
  MDAnalysis --- http://mdanalysis.googlecode.com
  Copyright (c) 2006-2011 Naveen Michaud-Agrawal,
                Elizabeth J. Denning, Oliver Beckstein,
                and contributors (see website for details)
  Released under the GNU Public Licence, v2 or any higher version

  Please cite your use of MDAnalysis in published work:

      N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
      O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
      Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
      in press.
*/


/* I have to fix error handling in these functions */

#include <Python.h>
#include <numpy/arrayobject.h>

#include "readdcd.h"

typedef struct {
  fio_fd fd;
  fio_size_t header_size;
  int natoms;
  int nsets;
  int setsread;
  int istart;
  int nsavc;
  double delta;
  int nfixed;
  int *freeind;
  float *fixedcoords;
  int reverse;
  int charmm;
  int first;
  int with_unitcell;
} dcdhandle;


int jump_to_frame(dcdhandle *dcd, int frame);

static PyObject *
__write_dcd_header(PyObject *self, PyObject *args)
{
  /* At this point we assume the file has been opened for writing */
  PyObject *temp = NULL;
  dcdhandle *dcd = NULL;
  fio_fd fd;
  int rc = 0;
  int natoms = 0;
  const char *remarks = "DCD";
  int istart = 0;             /* starting timestep of DCD file                  */
  int nsavc = 1;              /* number of timesteps between written DCD frames */
  double delta = 1.0;         /* length of a timestep                           */
  int with_unitcell = 1;      /* contains unit cell information (charmm format) */
  int charmm = DCD_IS_CHARMM; /* charmm-formatted DCD file                      */
  if (with_unitcell)
    charmm |= DCD_HAS_EXTRA_BLOCK;

  if (! self) {
    /* we were in fact called as a module function, try to retrieve 
       a matching object from args */
    if( !PyArg_ParseTuple(args, "Oi|iids", &self, &natoms, &istart, &nsavc, &delta, &remarks) )
      return NULL;
  } else {
    /* we were obviously called as an object method so args should 
       only have the int value. */
    if( !PyArg_ParseTuple(args, "i|iids", &natoms, &istart, &nsavc, &delta, &remarks) )
      return NULL;
  }

  // Get the file object from the class
  if (!PyObject_HasAttrString(self, "dcdfile")) {
    // Raise exception
    PyErr_SetString(PyExc_AttributeError, "dcdfile is not an attribute");
    return NULL;
  }
  
  if ((temp = PyObject_GetAttrString(self, "dcdfile")) == NULL) { // This gives me a New Reference
    // Raise exception
    PyErr_SetString(PyExc_AttributeError, "dcdfile is not an attribute");
    return NULL;
  }
  
  if (!PyFile_CheckExact(temp)) {
    // Raise exception
    PyErr_SetString(PyExc_TypeError, "dcdfile does not refer to a file object");
    Py_DECREF(temp);
    return NULL;
  }
  fd = fileno(PyFile_AsFile(temp)); 
  // No longer need the reference to temp
  Py_DECREF(temp);
  
  dcd = (dcdhandle *)malloc(sizeof(dcdhandle));
  memset(dcd, 0, sizeof(dcdhandle));
  dcd->fd = fd;

  if ((rc = write_dcdheader(dcd->fd, remarks, natoms, istart, nsavc, delta, with_unitcell, charmm)) < 0) {
    PyErr_SetString(PyExc_IOError, "Cannot write header of DCD file");
    free(dcd);
    return NULL;
  }

  dcd->natoms = natoms;
  dcd->nsets = 0;
  dcd->istart = istart;
  dcd->nsavc = nsavc;
  dcd->delta = delta;
  dcd->with_unitcell = with_unitcell;
  dcd->charmm = charmm;
  temp = PyCObject_FromVoidPtr(dcd, free); // Creates a New Reference
  if (PyObject_SetAttrString(self, "_dcd_C_ptr", temp) == -1) {
    // Raise exception - who knows what exception to raise??
    PyErr_SetString(PyExc_AttributeError, "Could not create attribute _dcd_C_ptr");
    Py_DECREF(temp);
    return NULL;
  }
  Py_DECREF(temp);
  
  // For debugging purposes
  temp = PyBuffer_FromMemory(dcd, sizeof(dcdhandle)); // Creates a New Reference
  if (PyObject_SetAttrString(self, "_dcd_C_str", temp) == -1) {
    // Raise exception - who knows what exception to raise??
    PyErr_SetString(PyExc_AttributeError, "Could not create attribute _dcd_C_str");
    Py_DECREF(temp);
    return NULL;
  }
  Py_DECREF(temp);
  
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *
__write_next_frame(PyObject *self, PyObject *args)
{
  dcdhandle *dcd;
  PyObject *temp;
  PyArrayObject *x, *y, *z, *uc;
  int rc, curstep;
  float* uc_array;
  double unitcell[6];
	
  if (!self) {
    /* we were in fact called as a module function, try to retrieve 
       a matching object from args */
    if( !PyArg_ParseTuple(args, "OO!O!O!O!", &self, &PyArray_Type, &x, &PyArray_Type, &y, &PyArray_Type, &z, &PyArray_Type, &uc) )
      return NULL;
  } else {
    /* we were obviously called as an object method so args should 
       only have the int value. */
    if( !PyArg_ParseTuple(args, "O!O!O!O!", &PyArray_Type, &x, &PyArray_Type, &y, &PyArray_Type, &z, &PyArray_Type, &uc) )
      return NULL;
  }

  if ((temp = PyObject_GetAttrString(self, "_dcd_C_ptr")) == NULL) { // This gives me a New Reference
    // Raise exception
    PyErr_SetString(PyExc_AttributeError, "_dcd_C_ptr is not an attribute");
    return NULL;
  }

  dcd = (dcdhandle*)PyCObject_AsVoidPtr(temp);
  Py_DECREF(temp);
  dcd->nsets++;
  curstep = dcd->istart + dcd->nsets * dcd->nsavc;

  uc_array = (float*) uc->data;
  unitcell[0] = uc_array[0];  /* A */
  unitcell[2] = uc_array[2];  /* B */
  unitcell[5] = uc_array[5];  /* C */
  /* write angle cosines with NAMD ordering [orbeckst] */
  /* (changed in MDAnalysis 0.9.0) */
  unitcell[4] = sin((M_PI_2 / 90.0 ) * (90.0 - uc_array[4]));  /* cos(alpha) */
  unitcell[3] = sin((M_PI_2 / 90.0 ) * (90.0 - uc_array[3]));  /* cos(beta) */
  unitcell[1] = sin((M_PI_2 / 90.0 ) * (90.0 - uc_array[1]));  /* cos(gamma) */

  if ((rc = write_dcdstep(dcd->fd, dcd->nsets, curstep, dcd->natoms, (float*)x->data, (float*)y->data, (float*)z->data,
			  dcd->with_unitcell ? unitcell: NULL,
			  dcd->charmm)) < 0)
    {
      PyErr_SetString(PyExc_IOError, "Could not write timestep to dcd file");
      return NULL;
    }

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *
__finish_dcd_write(PyObject *self, PyObject *args)
{
  //PyObject* temp;
  //dcdhandle *dcd;

  if (! self) { 
    /* we were in fact called as a module function, try to retrieve 
       a matching object from args */
    if( !PyArg_ParseTuple(args, "O", &self) )
      return NULL;
  } else { 
    /* we were obviously called as an object method so args should 
       only have the int value. */
    if( !PyArg_ParseTuple(args, "") )
      return NULL; 
  } 

  /*if ( !PyObject_HasAttrString(self, "_dcd_C_ptr") ) {
  // Raise exception
  PyErr_SetString(PyExc_AttributeError, "_dcd_C_ptr is not an attribute");
  return NULL;
  }
  
  if ((temp = PyObject_GetAttrString(self, "_dcd_C_ptr")) == NULL) { // This gives me a New Reference
  // Raise exception
  PyErr_SetString(PyExc_AttributeError, "_dcd_C_ptr is not an attribute");
  return NULL;
  } 
  dcd = (dcdhandle*)PyCObject_AsVoidPtr(temp);
  free(dcd);
  Py_DECREF(temp);*/

  if ( PyObject_DelAttrString(self, "_dcd_C_ptr") == -1) {
    // Raise exception
    PyErr_SetString(PyExc_AttributeError, "_dcd_C_ptr is not an attribute");
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *
__read_dcd_header(PyObject *self, PyObject *args)
{

  /* At this point we assume the file has been checked for existence and opened, and the DCD class
   * contains a reference to the file (which can be used to retrieve the file descriptor)
   */
  PyObject *temp;
  fio_fd fd;
  int rc;
  char *remarks = NULL;
  int len_remarks = 0;
  dcdhandle *dcd = NULL;

  if (! self) {
    /* we were in fact called as a module function, try to retrieve 
       a matching object from args */ 
    if( !PyArg_ParseTuple(args, "O", &self) ) 
      return NULL; 
  } else { 
    /* we were obviously called as an object method so args should 
       only have the int value. */ 
    if( !PyArg_ParseTuple(args, "") ) 
      return NULL; 
  } 

  // Get the file object from the class
  if (!PyObject_HasAttrString(self, "dcdfile")) {
    // Raise exception
    PyErr_SetString(PyExc_AttributeError, "dcdfile is not an attribute");
    return NULL;
  }

  if ((temp = PyObject_GetAttrString(self, "dcdfile")) == NULL) { // This gives me a New Reference
    // Raise exception
    PyErr_SetString(PyExc_AttributeError, "dcdfile is not an attribute");
    return NULL;
  }
	
  if (!PyFile_CheckExact(temp)) {
    // Raise exception
    PyErr_SetString(PyExc_TypeError, "dcdfile does not refer to a file object");
    goto error;
  }
  fd = fileno(PyFile_AsFile(temp));
  // No longer need the reference to temp
  Py_DECREF(temp);
	
  dcd = (dcdhandle *)malloc(sizeof(dcdhandle));
  memset(dcd, 0, sizeof(dcdhandle));
  dcd->fd = fd;

  if ((rc = read_dcdheader(dcd->fd, &dcd->natoms, &dcd->nsets, &dcd->istart,
			   &dcd->nsavc, &dcd->delta, &dcd->nfixed, &dcd->freeind,
			   &dcd->fixedcoords, &dcd->reverse, &dcd->charmm, &remarks, &len_remarks))) {
    // Raise exception
    PyErr_SetString(PyExc_IOError, "Cannot read DCD header");
    goto error;
  }

  /* 
   * Check that the file is big enough to really hold all the frames it claims to have
   */
  {
    off_t ndims, firstframesize, framesize, extrablocksize;
    off_t filesize;
    struct stat stbuf;

    extrablocksize = dcd->charmm & DCD_HAS_EXTRA_BLOCK ? 48 + 8 : 0;
    ndims = dcd->charmm & DCD_HAS_4DIMS ? 4 : 3;
    firstframesize = (dcd->natoms+2) * ndims * sizeof(float) + extrablocksize;
    framesize = (dcd->natoms-dcd->nfixed+2) * ndims * sizeof(float)
      + extrablocksize;

    // Get filesize from file descriptor
    memset(&stbuf, 0, sizeof(struct stat));
    if (fstat(dcd->fd, &stbuf)) {
      PyErr_SetString(PyExc_IOError, "Could not stat file");
      goto error;
    }

    // Size of header
    dcd->header_size = fio_ftell(dcd->fd);

    /* 
     * It's safe to use ftell, even though ftell returns a long, because the 
     * header size is < 4GB.
     */
    filesize = stbuf.st_size - fio_ftell(dcd->fd) - firstframesize;
    if (filesize < 0) {
      PyErr_SetString(PyExc_IOError, "DCD file appears to contain no timesteps");
      goto error;
    }
    dcd->nsets = filesize / framesize + 1;
    dcd->setsread = 0;
  }

  // We are at the first frame
  dcd->first = 1;

  temp = Py_BuildValue("s#", remarks, len_remarks);
  if (temp == NULL) goto error;
  if (PyObject_SetAttrString(self, "remarks", temp) == -1) {
    PyErr_SetString(PyExc_AttributeError, "Could not create attribute remarks");
    goto error;
  }
  Py_DECREF(temp);
  temp = Py_BuildValue("i", dcd->natoms);
  if (temp == NULL) goto error;
  if (PyObject_SetAttrString(self, "n_atoms", temp) == -1) {
    PyErr_SetString(PyExc_AttributeError, "Could not create attribute n_atoms");
    goto error;
  }
  Py_DECREF(temp);
  temp = Py_BuildValue("i", dcd->nsets);
  if (temp == NULL) goto error;
  if (PyObject_SetAttrString(self, "n_frames", temp) == -1) {
    PyErr_SetString(PyExc_AttributeError, "Could not create attribute n_frames");
    goto error;
  }
  Py_DECREF(temp);
  temp = Py_BuildValue("i", dcd->nfixed);
  if (temp == NULL) goto error;
  if (PyObject_SetAttrString(self, "fixed", temp) == -1) {
    PyErr_SetString(PyExc_AttributeError, "Could not create attribute fixed");
    goto error;
  }
  Py_DECREF(temp);
  temp = Py_BuildValue("i", dcd->istart);
  if (temp == NULL) goto error;
  if (PyObject_SetAttrString(self, "start_timestep", temp) == -1) {
    PyErr_SetString(PyExc_AttributeError, "Could not create attribute fixed");
    goto error;
  }
  Py_DECREF(temp);
  temp = Py_BuildValue("i", dcd->nsavc);
  if (temp == NULL) goto error;
  if (PyObject_SetAttrString(self, "skip_timestep", temp) == -1) {
    PyErr_SetString(PyExc_AttributeError, "Could not create attribute fixed");
    goto error;
  }
  Py_DECREF(temp);
  temp = Py_BuildValue("d", dcd->delta);
  if (temp == NULL) goto error;
  if (PyObject_SetAttrString(self, "delta", temp) == -1) {
    PyErr_SetString(PyExc_AttributeError, "Could not create attribute fixed");
    goto error;
  }
  Py_DECREF(temp);
	
  temp = PyCObject_FromVoidPtr(dcd, NULL);
  if (temp == NULL) goto error;
  if (PyObject_SetAttrString(self, "_dcd_C_ptr", temp) == -1) {
    PyErr_SetString(PyExc_AttributeError, "Could not create attribute _dcd_C_ptr");
    goto error;
  }

  if ((dcd->charmm & DCD_IS_CHARMM) && (dcd->charmm & DCD_HAS_EXTRA_BLOCK)) {
    dcd->with_unitcell = 1;
    temp = Py_True;
  } else {
    temp = Py_False;
  }
  Py_INCREF(temp);
  if (PyObject_SetAttrString(self, "periodic", temp) == -1) {
    PyErr_SetString(PyExc_AttributeError, "Could not create attribute periodic");
    goto error;
  }
  Py_DECREF(temp);
	
  // For debugging purposes
  temp = PyBuffer_FromMemory(dcd, sizeof(dcdhandle));
  if (temp == NULL) goto error;
  if (PyObject_SetAttrString(self, "_dcd_C_str", temp) == -1) {
    PyErr_SetString(PyExc_AttributeError, "Could not create attribute _dcd_C_str");
    goto error;
  }
  Py_DECREF(temp);

  Py_INCREF(Py_None);
  return Py_None;

 error:
  Py_XDECREF(temp);
  if (dcd != NULL) free(dcd);
  if (remarks != NULL) free(remarks);
  return NULL;
}

static PyObject *
__read_timeseries(PyObject *self, PyObject *args)
{
  PyObject *temp = NULL;
  PyArrayObject *coord = NULL;
  PyListObject *atoms = NULL;
  int lowerb = 0, upperb = 0, range=0;
  float *tempX = NULL, *tempY = NULL, *tempZ = NULL;
  int rc;
  int i, j, index;
  int n_atoms = 0, n_frames = 0;
  int start = 0, stop = -1, skip = 1, numskip = 0;
  dcdhandle *dcd = NULL;
  int *atomlist = NULL;
  npy_intp dimensions[3];
  float unitcell[6];
  const char* format = "afc";
	
  if (!self) {
    /* we were in fact called as a module function, try to retrieve 
       a matching object from args */
    if( !PyArg_ParseTuple(args, "OO!|iiis", &self, &PyList_Type, &atoms, &start, &stop, &skip, &format) )
      return NULL; 
  } else {
    /* we were obviously called as an object method so args should 
       only have the int value. */
    if( !PyArg_ParseTuple(args, "O!|iiis", &PyList_Type, &atoms, &start, &stop, &skip, &format) )
      return NULL;
  }

  if ((temp = PyObject_GetAttrString(self, "_dcd_C_ptr")) == NULL) { // This gives me a New Reference
    // Raise exception
    PyErr_SetString(PyExc_AttributeError, "_dcd_C_ptr is not an attribute");
    return NULL;
  }

  dcd = (dcdhandle*)PyCObject_AsVoidPtr(temp);
  Py_DECREF(temp);

  // Assume that start and stop are valid
  if (stop == -1) { stop = dcd->nsets; }
  n_frames = (stop-start+1) / skip;
  //n_frames = dcd->nsets / skip;
  n_atoms = PyList_Size((PyObject*)atoms);
  if (n_atoms == 0) {
    PyErr_SetString(PyExc_Exception, "No atoms passed into _read_timeseries function");
    return NULL;
  }
  atomlist = (int*)malloc(sizeof(int)*n_atoms);
  memset(atomlist, 0, sizeof(int)*n_atoms);

  // Get the atom indexes
  for (i=0;i<n_atoms;i++) {
    temp = PyList_GetItem((PyObject*)atoms, i); // Borrowed Reference
    if (temp==NULL) goto error;
    // Make sure temp is an integer
    /* TODO: this is not the proper check but [OB] cannot figure out how
       to check if this is a numpy.int64 or similar; PyInt_Check would fail
       on those (Issue 18)
    */
    if (!PyArray_IsAnyScalar((PyObject*)temp)) {
      PyErr_SetString(PyExc_ValueError, "Atom number is not an integer");
      goto error;
    }
    atomlist[i] = PyInt_AsLong(temp);
  }

  lowerb = atomlist[0]; upperb = atomlist[n_atoms-1];
  range = upperb-lowerb+1;
	
  // Figure out the format string
  if (strncasecmp(format, "afc", 3) == 0) {
    dimensions[0] = n_atoms; dimensions[1] = n_frames; dimensions[2] = 3;
  } else
    if (strncasecmp(format, "acf", 3) == 0) {
      dimensions[0] = n_atoms; dimensions[1] = 3; dimensions[2] = n_frames;
    } else
      if (strncasecmp(format, "fac", 3) == 0) {
	dimensions[0] = n_frames; dimensions[1] = n_atoms; dimensions[2] = 3;
      } else
	if (strncasecmp(format, "fca", 3) == 0) {
	  dimensions[0] = n_frames; dimensions[1] = 3; dimensions[2] = n_atoms;
	} else
	  if (strncasecmp(format, "caf", 3) == 0) {
	    dimensions[0] = 3; dimensions[1] = n_atoms; dimensions[2] = n_frames;
	  } else
	    if (strncasecmp(format, "cfa", 3) == 0) {
	      dimensions[0] = 3; dimensions[1] = n_frames; dimensions[2] = n_atoms;
	    }
	
  coord = (PyArrayObject*) PyArray_SimpleNew(3, dimensions, NPY_DOUBLE);
  if (coord == NULL) goto error;

  // Reset trajectory
  rc = fio_fseek(dcd->fd, dcd->header_size, FIO_SEEK_SET);
  dcd->setsread = 0;
  dcd->first = 1;

  // Jump to starting frame
  jump_to_frame(dcd, start);
	
  // Now read through trajectory and get the atom pos
  tempX = (float*)malloc(sizeof(float)*range);
  tempY = (float*)malloc(sizeof(float)*range);
  tempZ = (float*)malloc(sizeof(float)*range);
  if ((tempX == NULL) | (tempY == NULL) || (tempZ == NULL)) {
    PyErr_SetString(PyExc_MemoryError, "Can't allocate temporary space for coordinate arrays");
    goto error;
  }
	
  for (i=0;i<n_frames;i++)
    {
      if (skip > 1) {
	/*if (dcd->first && dcd->nfixed) {
	  rc = read_dcdsubset(dcd->fd, dcd->natoms, lowerb, upperb, tempX, tempY, tempZ,
	  unitcell, dcd->nfixed, dcd->first, dcd->freeind, dcd->fixedcoords,
	  dcd->reverse, dcd->charmm);
	  dcd->first=0;
	  if (rc < 0) {
          // return an exception
          PyErr_SetString(PyExc_IOError, "Error reading first frame from DCD file");
          //fprintf(stderr, "read_dcdstep returned %d\n", rc);
          return NULL;
          //return MOLFILE_ERROR;
	  }
	  dcd->setsread++;
	  }*/
	// Figure out how many steps to skip
	numskip = skip - (dcd->setsread % skip) - 1;
	rc = skip_dcdstep(dcd->fd, dcd->natoms, dcd->nfixed, dcd->charmm, numskip);
    	if (rc < 0) {
	  // return an exception
	  PyErr_SetString(PyExc_IOError, "Error skipping frame from DCD file");
	  goto error;
    	}
    	dcd->setsread+=numskip;                                                                              
      }
      rc = read_dcdsubset(dcd->fd, dcd->natoms, lowerb, upperb, tempX, tempY, tempZ,
			  unitcell, dcd->nfixed, dcd->first, dcd->freeind, dcd->fixedcoords,
			  dcd->reverse, dcd->charmm);
      dcd->first = 0;
      dcd->setsread++;
      if (rc < 0) {
    	// return an exception
    	PyErr_SetString(PyExc_IOError, "Error reading frame from DCD file");
    	goto error;
      }
      // Copy into Numeric array only those atoms we are interested in
      for (j=0;j<n_atoms;j++) {
	index = atomlist[j]-lowerb;
	/* 
	 * coord[a][b][c] = *(float*)(coord->data + a*coord->strides[0] + b*coord->strides[1] + c*coord->strides[2])
	 */
	if (strncasecmp(format, "afc", 3) == 0) {
	  *(double*)(coord->data + j*coord->strides[0] + i*coord->strides[1] + 0*coord->strides[2]) = tempX[index];
	  *(double*)(coord->data + j*coord->strides[0] + i*coord->strides[1] + 1*coord->strides[2]) = tempY[index];
	  *(double*)(coord->data + j*coord->strides[0] + i*coord->strides[1] + 2*coord->strides[2]) = tempZ[index];
	} else
	  if (strncasecmp(format, "acf", 3) == 0) {
	    *(double*)(coord->data + j*coord->strides[0] + 0*coord->strides[1] + i*coord->strides[2]) = tempX[index];
	    *(double*)(coord->data + j*coord->strides[0] + 1*coord->strides[1] + i*coord->strides[2]) = tempY[index];
	    *(double*)(coord->data + j*coord->strides[0] + 2*coord->strides[1] + i*coord->strides[2]) = tempZ[index];
	  } else
	    if (strncasecmp(format, "fac", 3) == 0) {
	      *(double*)(coord->data + i*coord->strides[0] + j*coord->strides[1] + 0*coord->strides[2]) = tempX[index];
	      *(double*)(coord->data + i*coord->strides[0] + j*coord->strides[1] + 1*coord->strides[2]) = tempY[index];
	      *(double*)(coord->data + i*coord->strides[0] + j*coord->strides[1] + 2*coord->strides[2]) = tempZ[index];
	    } else
	      if (strncasecmp(format, "fca", 3) == 0) {
		*(double*)(coord->data + i*coord->strides[0] + 0*coord->strides[1] + j*coord->strides[2]) = tempX[index];
		*(double*)(coord->data + i*coord->strides[0] + 1*coord->strides[1] + j*coord->strides[2]) = tempY[index];
		*(double*)(coord->data + i*coord->strides[0] + 2*coord->strides[1] + j*coord->strides[2]) = tempZ[index];
	      } else
  		if (strncasecmp(format, "caf", 3) == 0) {
		  *(double*)(coord->data + 0*coord->strides[0] + j*coord->strides[1] + i*coord->strides[2]) = tempX[index];
		  *(double*)(coord->data + 1*coord->strides[0] + j*coord->strides[1] + i*coord->strides[2]) = tempY[index];
		  *(double*)(coord->data + 2*coord->strides[0] + j*coord->strides[1] + i*coord->strides[2]) = tempZ[index];
  		} else
		  if (strncasecmp(format, "cfa", 3) == 0) {
		    *(double*)(coord->data + 0*coord->strides[0] + i*coord->strides[1] + j*coord->strides[2]) = tempX[index];
		    *(double*)(coord->data + 1*coord->strides[0] + i*coord->strides[1] + j*coord->strides[2]) = tempY[index];
		    *(double*)(coord->data + 2*coord->strides[0] + i*coord->strides[1] + j*coord->strides[2]) = tempZ[index];
		  }
      }
      // Check if we've been interupted by the user
      if (PyErr_CheckSignals() == 1) goto error;
    }

  // Reset trajectory
  rc = fio_fseek(dcd->fd, dcd->header_size, FIO_SEEK_SET);
  dcd->setsread = 0;
  dcd->first = 1;	
  free(atomlist);
  free(tempX);
  free(tempY);
  free(tempZ);
  return PyArray_Return(coord);

 error:
  // Reset trajectory
  rc = fio_fseek(dcd->fd, dcd->header_size, FIO_SEEK_SET);
  dcd->setsread = 0;
  dcd->first = 1;	
  Py_XDECREF(coord);
  if (atomlist != NULL) free(atomlist);
  if (tempX != NULL) free(tempX);
  if (tempY != NULL) free(tempY);
  if (tempZ != NULL) free(tempZ);
  return NULL;
}

static PyObject *
__read_next_frame(PyObject *self, PyObject *args)
{
  dcdhandle *dcd;
  PyObject *temp;
  PyArrayObject *x, *y, *z, *uc;
  int skip=1;
  int rc,numskip;
  float* unitcell;
  float alpha, beta, gamma;
	
  if (!self) {
    /* we were in fact called as a module function, try to retrieve 
       a matching object from args */
    if( !PyArg_ParseTuple(args, "OO!O!O!O!|i", &self, &PyArray_Type, &x, &PyArray_Type, &y, &PyArray_Type, &z, &PyArray_Type, &uc, &skip) )
      return NULL; 
  } else {
    /* we were obviously called as an object method so args should 
       only have the int value. */
    if( !PyArg_ParseTuple(args, "O!O!O!O!|i", &PyArray_Type, &x, &PyArray_Type, &y, &PyArray_Type, &z, &PyArray_Type, &uc, &skip) ) 
      return NULL;
  }

  if ((temp = PyObject_GetAttrString(self, "_dcd_C_ptr")) == NULL) { // This gives me a New Reference
    // Raise exception
    PyErr_SetString(PyExc_AttributeError, "_dcd_C_ptr is not an attribute");
    return NULL;
  }

  dcd = (dcdhandle*)PyCObject_AsVoidPtr(temp);
  Py_DECREF(temp);

  unitcell = (float*) uc->data;
  unitcell[0] = unitcell[2] = unitcell[5] = 0.0f;
  unitcell[1] = unitcell[3] = unitcell[4] = 90.0f;                                                                     

  /* Check for EOF here; that way all EOF's encountered later must be errors */
  if (dcd->setsread == dcd->nsets) {
    // Is this an exception?
    PyErr_SetString(PyExc_IOError, "End of file reached for dcd file");
    return NULL;
    //return MOLFILE_EOF;
  }
  if (skip > 1) {
    if (dcd->first && dcd->nfixed) {
      /* We can't just skip it because we need the fixed atom coordinates */
      rc = read_dcdstep(dcd->fd, dcd->natoms, (float*)x->data, (float*)y->data, (float*)z->data,
			unitcell, dcd->nfixed, dcd->first, dcd->freeind, dcd->fixedcoords,
			dcd->reverse, dcd->charmm);
      dcd->first = 0;
      if (rc < 0) {
        // return an exception
        PyErr_SetString(PyExc_IOError, "Error reading first frame from DCD file");
        //fprintf(stderr, "read_dcdstep returned %d\n", rc);
        return NULL;
        //return MOLFILE_ERROR;
      }
      dcd->setsread++;
      temp = Py_BuildValue("i", dcd->setsread);
      return temp;
      //return rc; /* XXX this needs to be updated */
    }
    dcd->first = 0;
    // Figure out how many steps to skip
    numskip = skip - (dcd->setsread % skip) - 1;
    /* XXX this needs to be changed */
    rc = skip_dcdstep(dcd->fd, dcd->natoms, dcd->nfixed, dcd->charmm, numskip);
    if (rc < 0) {
      // return an exception
      PyErr_SetString(PyExc_IOError, "Error skipping frame from DCD file");
      //fprintf(stderr, "read_dcdstep returned %d\n", rc);
      return NULL;
      //return MOLFILE_ERROR;
    }
    dcd->setsread+=numskip;
  }
  rc = read_dcdstep(dcd->fd, dcd->natoms, (float*)x->data, (float*)y->data, (float*)z->data, unitcell,
		    dcd->nfixed, dcd->first, dcd->freeind, dcd->fixedcoords,
		    dcd->reverse, dcd->charmm);
  dcd->first = 0;
  dcd->setsread++;
  if (rc < 0) {
    // return an exception
    PyErr_SetString(PyExc_IOError, "Error reading frame from DCD file");
    //fprintf(stderr, "read_dcdstep returned %d\n", rc);
    return NULL;
    //return MOLFILE_ERROR;
  }                                                                                            

  if (unitcell[1] >= -1.0 && unitcell[1] <= 1.0 &&
      unitcell[3] >= -1.0 && unitcell[3] <= 1.0 &&
      unitcell[4] >= -1.0 && unitcell[4] <= 1.0) {
    /* This file was generated by Charmm, or by NAMD > 2.5, with the angle */
    /* cosines of the periodic cell angles written to the DCD file.        */
    /* This formulation improves rounding behavior for orthogonal cells    */
    /* so that the angles end up at precisely 90 degrees, unlike acos().   */
    /* (changed in MDAnalysis 0.9.0 to have NAMD ordering of the angles;   */
    /* see Issue 187) */
    alpha = 90.0 - asin(unitcell[4]) * 90.0 / M_PI_2;
    beta  = 90.0 - asin(unitcell[3]) * 90.0 / M_PI_2;
    gamma = 90.0 - asin(unitcell[1]) * 90.0 / M_PI_2;
  } else {
    /* This file was likely generated by NAMD 2.5 and the periodic cell    */
    /* angles are specified in degrees rather than angle cosines.          */
    alpha = unitcell[4];
    beta  = unitcell[3];
    gamma = unitcell[1];
  }
  unitcell[4] = alpha;
  unitcell[3] = beta;
  unitcell[1] = gamma;
	
  // Return the frame read
  temp = Py_BuildValue("i", dcd->setsread);
  return temp;
}

static PyObject *
__finish_dcd_read(PyObject *self, PyObject *args)
{
  PyObject* temp;
  dcdhandle *dcd;

  if (! self) {
    /* we were in fact called as a module function, try to retrieve 
       a matching object from args */ 
    if( !PyArg_ParseTuple(args, "O", &self) )
      return NULL; 
  } else { 
    /* we were obviously called as an object method so args should 
       only have the int value. */ 
    if( !PyArg_ParseTuple(args, "") )
      return NULL; 
  }
 
  if ( !PyObject_HasAttrString(self, "_dcd_C_ptr") ) {
    // Raise exception
    PyErr_SetString(PyExc_AttributeError, "_dcd_C_ptr is not an attribute");
    return NULL;
  }

  if ((temp = PyObject_GetAttrString(self, "_dcd_C_ptr")) == NULL) { // This gives me a New Reference
    // Raise exception
    PyErr_SetString(PyExc_AttributeError, "_dcd_C_ptr is not an attribute");
    return NULL;
  }

  dcd = (dcdhandle*)PyCObject_AsVoidPtr(temp);
	
  close_dcd_read(dcd->freeind, dcd->fixedcoords);
  free(dcd);
  Py_DECREF(temp);
  Py_INCREF(Py_None);
  return Py_None;
}                                                                                                                      
int jump_to_frame(dcdhandle *dcd, int frame)
{
  int rc;
  if (frame > dcd->nsets) {
    return -1;
  }
  // Calculate file offset
  {
    off_t extrablocksize, ndims, firstframesize, framesize;
    off_t pos;
    extrablocksize = dcd->charmm & DCD_HAS_EXTRA_BLOCK ? 48 + 8 : 0;
    ndims = dcd->charmm & DCD_HAS_4DIMS ? 4 : 3;
    firstframesize = (dcd->natoms+2) * ndims * sizeof(float) + extrablocksize;
    framesize = (dcd->natoms-dcd->nfixed+2) * ndims * sizeof(float)
      + extrablocksize;
    // Use zero indexing
    if (frame == 0) {
      pos = dcd->header_size;
      dcd->first = 1;
    }
    else {
      pos = dcd->header_size + firstframesize + framesize * (frame-1);
      dcd->first = 0;
    }
    rc = fio_fseek(dcd->fd, pos, FIO_SEEK_SET);
  }
  dcd->setsread = frame;
  return rc;
}

static PyObject *
__jump_to_frame(PyObject *self, PyObject *args)
{
  PyObject* temp;
  dcdhandle *dcd;
  int frame;

  if (! self) {
    /* we were in fact called as a module function, try to retrieve 
       a matching object from args */
    if( !PyArg_ParseTuple(args, "Oi", &self, &frame) )
      return NULL;
  } else {
    /* we were obviously called as an object method so args should 
       only have the int value. */
    if( !PyArg_ParseTuple(args, "i", &frame) )
      return NULL;
  }             
  
  if ( !PyObject_HasAttrString(self, "_dcd_C_ptr") ) {
    // Raise exception
    PyErr_SetString(PyExc_AttributeError, "_dcd_C_ptr is not an attribute");
    return NULL;
  } 
  
  if ((temp = PyObject_GetAttrString(self, "_dcd_C_ptr")) == NULL) { // This gives me a New Reference
    // Raise exception
    PyErr_SetString(PyExc_AttributeError, "_dcd_C_ptr is not an attribute");
    return NULL;
  }
  dcd = (dcdhandle*)PyCObject_AsVoidPtr(temp);
  Py_DECREF(temp);
	
  /*if (frame > dcd->nsets) {
    PyErr_SetString(PyExc_IndexError, "Invalid frame");
    return NULL;
    }

    // Calculate file offset
    {
    off_t extrablocksize, ndims, firstframesize, framesize;
    off_t pos;
    extrablocksize = dcd->charmm & DCD_HAS_EXTRA_BLOCK ? 48 + 8 : 0;
    ndims = dcd->charmm & DCD_HAS_4DIMS ? 4 : 3;
    firstframesize = (dcd->natoms+2) * ndims * sizeof(float) + extrablocksize;
    framesize = (dcd->natoms-dcd->nfixed+2) * ndims * sizeof(float) 
    + extrablocksize;
    // Use zero indexing
    if (frame == 0) {
    pos = dcd->header_size;
    }
    else {
    pos = dcd->header_size + firstframesize + framesize * (frame-1);
    }
    rc = fio_fseek(dcd->fd, pos, FIO_SEEK_SET);
    }
    dcd->setsread = frame;
    dcd->first = 0;*/
  jump_to_frame(dcd, frame);

  temp = Py_BuildValue("i", dcd->setsread);
  return temp;
}

static PyObject *
__reset_dcd_read(PyObject *self, PyObject *args)
{
  PyObject* temp;
  dcdhandle *dcd;
  int rc;
	
  if (! self) {
    /* we were in fact called as a module function, try to retrieve 
       a matching object from args */
    if( !PyArg_ParseTuple(args, "O", &self) )
      return NULL;
  } else {
    /* we were obviously called as an object method so args should 
       only have the int value. */
    if( !PyArg_ParseTuple(args, "") )
      return NULL;
  }                                                                                                                    

  if ( !PyObject_HasAttrString(self, "_dcd_C_ptr") ) {
    // Raise exception
    PyErr_SetString(PyExc_AttributeError, "_dcd_C_ptr is not an attribute");
    return NULL;
  }

  if ((temp = PyObject_GetAttrString(self, "_dcd_C_ptr")) == NULL) { // This gives me a New Reference
    // Raise exception
    PyErr_SetString(PyExc_AttributeError, "_dcd_C_ptr is not an attribute");
    return NULL;
  }

  dcd = (dcdhandle*)PyCObject_AsVoidPtr(temp);
  rc = fio_fseek(dcd->fd, dcd->header_size, FIO_SEEK_SET);
  dcd->setsread = 0;
  dcd->first = 1;
  Py_DECREF(temp);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyMethodDef DCDMethods[] = {
  {"__write_dcd_header", __write_dcd_header, METH_VARARGS, "Write a DCD header."},
  {"__write_next_frame", __write_next_frame, METH_VARARGS, "Write the next timestep."},
  {"__finish_dcd_write", __finish_dcd_write, METH_VARARGS, "Clean up and data for handling dcd file."},
  {"__read_dcd_header", __read_dcd_header, METH_VARARGS, "Read in a DCD header."},
  {"__read_next_frame", __read_next_frame, METH_VARARGS, "Read in the next timestep."},
  {"__jump_to_frame", __jump_to_frame, METH_VARARGS, "Jump to specified timestep."},
  {"__reset_dcd_read", __reset_dcd_read, METH_VARARGS, "Reset dcd file reading."},
  {"__finish_dcd_read", __finish_dcd_read, METH_VARARGS, "Clean up any data for handling dcd file."},
  {"__read_timeseries", __read_timeseries, METH_VARARGS, "Return a Numpy array of the coordinates for a set of atoms for the whole trajectory."},
  {NULL, NULL, 0, NULL}	/* Sentinel */
};

PyMODINIT_FUNC
init_dcdmodule(void)
{
  (void) Py_InitModule("_dcdmodule", DCDMethods);
  import_array();
}

