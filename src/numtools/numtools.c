
#include <math.h>
#include <Python.h>
#include <Numeric/arrayobject.h>

static PyObject *
length(PyObject *self, PyObject *args)
{
  PyArrayObject *array = NULL, *res = NULL;
	int *dim;
	int nd,i,j;
  double x, y, z, r;
	
	// Like the numeric ufuncs, can pass in a destination array for in place
  if( !PyArg_ParseTuple(args, "O!|O!", &PyArray_Type, &array, &PyArray_Type, &res) ) {
		return NULL;
	}

	// Incr ref count if not null over here, that way if we have an exception
	// we can clean it up properly
	Py_XINCREF(res);

	nd = array->nd;
	dim = (int*)malloc(sizeof(int)*nd);
	memcpy(dim, array->dimensions, sizeof(int)*nd);
	
  
	// Check to make sure the last dimension is of size 3	
	if (nd != 3) {
		PyErr_SetString(PyExc_Exception, "input array requires 3 dimensions");
		goto error;
	}
	if (dim[nd-1] != 3) {
		PyErr_SetString(PyExc_Exception, "input array does not have the proper shape (...,3)");
		goto error;
	}
	if (array->descr->type_num != PyArray_DOUBLE) {
		PyErr_SetString(PyExc_TypeError, "input array has incorrect type");
		goto error;
	}

	// If a destination array is given, make sure it is the right typecode/size
	dim[nd-1] = 1;
	if (res != NULL) {
		if (res->descr->type_num != PyArray_DOUBLE) {
			PyErr_SetString(PyExc_TypeError, "result array has incorrect type");
			goto error;
		}
	  if (res->nd == nd) {
			for (i=0;i<nd;i++) {
				if (res->dimensions[i] != dim[i]){
					PyErr_SetString(PyExc_Exception, "result array does not have the proper shape");
					goto error;
				}
			}
		}	else {
			PyErr_SetString(PyExc_Exception, "result array does not have the proper dimensions");
			goto error;
		}
	} else { // Create a new array
		res = (PyArrayObject*) PyArray_FromDims(nd, dim, PyArray_DOUBLE);
		if (res == NULL) goto error;
	}

	for (i=0;i<dim[0];i++) {
		for (j=0;j<dim[1];j++) {
			x = *(double*)(array->data+i*array->strides[0]+j*array->strides[1]+0*array->strides[2]);
			y = *(double*)(array->data+i*array->strides[0]+j*array->strides[1]+1*array->strides[2]);
			z = *(double*)(array->data+i*array->strides[0]+j*array->strides[1]+2*array->strides[2]);
			r = sqrt((x*x)+(y*y)+(z*z));
			*(double*)(res->data+i*res->strides[0]+j*res->strides[1]+0*res->strides[2]) = (double)r;
		}
	}
	free(dim);
	return PyArray_Return(res);

  error:
	Py_XDECREF(res);
	free(dim);
	return NULL;
}

static PyObject *
normalize(PyObject *self, PyObject *args)
{
  PyArrayObject *array = NULL, *res = NULL;
	int *dim;
	int nd,i,j;
  double x, y, z, r;
	
	// Like the numeric ufuncs, can pass in a destination array for in place
  if( !PyArg_ParseTuple(args, "O!|O!", &PyArray_Type, &array, &PyArray_Type, &res) ) {
		return NULL;
	}

  // Incr ref count if not null over here, that way if we have an exception
  // we can clean it up properly
  Py_XINCREF(res);                                                                                            

	nd = array->nd;
	dim = (int*)malloc(sizeof(int)*nd);
	memcpy(dim, array->dimensions, sizeof(int)*nd);
  // Check to make sure the last dimension is of size 3	
	if (dim[nd-1] != 3) {
		PyErr_SetString(PyExc_Exception, "input array does not have the proper shape (...,3)");
		goto error;
	}
  if (array->descr->type_num != PyArray_DOUBLE) {
    PyErr_SetString(PyExc_TypeError, "input array has incorrect type");
    goto error;
  }	
	
	// If a destination array is given, make sure it is the right typecode/size
	if (res != NULL) {
    if (res->descr->type_num != PyArray_DOUBLE) {
      PyErr_SetString(PyExc_TypeError, "result array has incorrect type");
      goto error;
    }
		if (res->nd == nd) {
			for (i=0;i<nd;i++) {
				if (res->dimensions[i] != dim[i]){
					PyErr_SetString(PyExc_Exception, "result array does not have the proper shape");
					goto error;
				}
			}
		}	else {
			PyErr_SetString(PyExc_Exception, "result array does not have the proper dimensions");
			goto error;
		}
	} else { // Create a new array
		res = (PyArrayObject*) PyArray_FromDims(nd, dim, PyArray_DOUBLE);
		if (res == NULL) goto error;
	}

	// Now do the calculation
	for (i=0;i<dim[0];i++) {
		for (j=0;j<dim[1];j++) {
			x = *(double*)(array->data+i*array->strides[0]+j*array->strides[1]+0*array->strides[2]);
			y = *(double*)(array->data+i*array->strides[0]+j*array->strides[1]+1*array->strides[2]);
			z = *(double*)(array->data+i*array->strides[0]+j*array->strides[1]+2*array->strides[2]);
			r = (double)sqrt((x*x)+(y*y)+(z*z));
      *(double*)(res->data+i*res->strides[0]+j*res->strides[1]+0*res->strides[2]) = (double)x/r;
      *(double*)(res->data+i*res->strides[0]+j*res->strides[1]+1*res->strides[2]) = (double)y/r;
      *(double*)(res->data+i*res->strides[0]+j*res->strides[1]+2*res->strides[2]) = (double)z/r;
		}
	}
	free(dim);
	return PyArray_Return(res);

  error:
	Py_XDECREF(res);
	free(dim);
	return NULL;
}

static PyObject *
distance(PyObject *self, PyObject *args)
{
  PyArrayObject *array = NULL, *res = NULL;
	int *dim;
	int nd,i,j;
  double x, y, z, d;
  double xpos, ypos, zpos;
	
	// Like the numeric ufuncs, can pass in a destination array for in place
  if( !PyArg_ParseTuple(args, "O!ddd|O!", &PyArray_Type, &array, &xpos, &ypos, &zpos, &PyArray_Type, &res) ) {
		return NULL;
	}

  // Incr ref count if not null over here, that way if we have an exception
  // we can clean it up properly
  Py_XINCREF(res);                                                                                            

	nd = array->nd;
	dim = (int*)malloc(sizeof(int)*nd);
	memcpy(dim, array->dimensions, sizeof(int)*nd);
  // Check to make sure the last dimension is of size 3	
	if (dim[nd-1] != 3) {
		PyErr_SetString(PyExc_Exception, "input array does not have the proper shape (...,3)");
		goto error;
	}
  if (array->descr->type_num != PyArray_DOUBLE) {
    PyErr_SetString(PyExc_TypeError, "input array has incorrect type");
    goto error;
  }	
	
	// If a destination array is given, make sure it is the right typecode/size
	dim[nd-1] = 1;
	if (res != NULL) {
    if (res->descr->type_num != PyArray_DOUBLE) {
      PyErr_SetString(PyExc_TypeError, "result array has incorrect type");
      goto error;
    }
		if (res->nd == nd) {
			for (i=0;i<nd;i++) {
				if (res->dimensions[i] != dim[i]){
					PyErr_SetString(PyExc_Exception, "result array does not have the proper shape");
					goto error;
				}
			}
		}	else {
			PyErr_SetString(PyExc_Exception, "result array does not have the proper dimensions");
			goto error;
		}
	} else { // Create a new array
		res = (PyArrayObject*) PyArray_FromDims(nd, dim, PyArray_DOUBLE);
		if (res == NULL) goto error;
	}

	// Now do the calculation
	for (i=0;i<dim[0];i++) {
		for (j=0;j<dim[1];j++) {
			x = *(double*)(array->data+i*array->strides[0]+j*array->strides[1]+0*array->strides[2]);
			y = *(double*)(array->data+i*array->strides[0]+j*array->strides[1]+1*array->strides[2]);
			z = *(double*)(array->data+i*array->strides[0]+j*array->strides[1]+2*array->strides[2]);
			d = (double)sqrt((x-xpos)*(x-xpos)+(y-ypos)*(y-ypos)+(z-zpos)*(z-zpos));
      *(double*)(res->data+i*res->strides[0]+j*res->strides[1]+0*res->strides[2]) = (double)d;
		}
	}
	free(dim);
	return PyArray_Return(res);

  error:
	Py_XDECREF(res);
	free(dim);
	return NULL;
}


static PyMethodDef NumericMethods[] = {
	{"normalize", normalize, METH_VARARGS, "Given an array of vectors returns a normalized array"},
	{"length", length, METH_VARARGS, "Computes |v| for array[:,:,3]"},
	{"distance", distance, METH_VARARGS, "Computes the distance between an array of vectors and a point (x,y,z)"},
	{NULL, NULL, 0, NULL}	/* Sentinel */
};

PyMODINIT_FUNC
init_numtools(void)
{
	(void) Py_InitModule("_numtools", NumericMethods);
	import_array();
}

