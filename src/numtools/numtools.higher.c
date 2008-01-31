
#include <math.h>
#include <Python.h>
#include <Numeric/arrayobject.h>

typedef struct {
	void (*func)(PyArrayObject*,PyArrayObject*,int*);
	PyArrayObject* array;
	PyArrayObject* res;
} callback;

void call_on_indices(int* list, int size, void*c);

void __length(PyArrayObject* array, PyArrayObject* res, int* ind)
{
	double x, y, z, r;
	void *apos, *rpos = NULL;
	int nd, i;
	// calculate the position of the data to access
	nd = array->nd;
	apos = array->data;
	rpos = res->data;
	for (i=0;i<nd-1;i++) {
		apos += ind[i]*array->strides[i];
		rpos += ind[i]*res->strides[i];
	}
	x = *(double*)(apos+0*array->strides[nd-1]);
	y = *(double*)(apos+1*array->strides[nd-1]);
	z = *(double*)(apos+2*array->strides[nd-1]);
	r = sqrt((x*x)+(y*y)+(z*z));
	*(double*)(rpos+0*res->strides[nd-1]) = (double)r;
}

void __norm(PyArrayObject* array, PyArrayObject* res, int* ind)
{
	double x, y, z, r;
	void *apos, *rpos = NULL;
	int nd, i;
	// calculate the position of the data to access
	nd = array->nd;
	apos = array->data;
	rpos = res->data;
	for (i=0;i<nd-1;i++) {
		apos += ind[i]*array->strides[i];
		rpos += ind[i]*res->strides[i];
	}
	x = *(double*)(apos+0*array->strides[nd-1]);
	y = *(double*)(apos+1*array->strides[nd-1]);
	z = *(double*)(apos+2*array->strides[nd-1]);
	r = sqrt((x*x)+(y*y)+(z*z));
	*(double*)(rpos+0*res->strides[nd-1]) = (double)x/r;
	*(double*)(rpos+1*res->strides[nd-1]) = (double)y/r;
	*(double*)(rpos+2*res->strides[nd-1]) = (double)z/r;
}

// Nothing like trying to do higher-order control structures in c
// this generalizes a nested for loop, nesting 'size' number of times
void __rec_indices(int* ind, int* list, int size, int pos, void* c) {
  int i;
	callback* cb;
  if (pos == size-1) {
    // call function on indices
		cb = (callback*) c;
		//printf("calling func - ");
		//for (i=0;i<size;i++)
		//	printf("%d ", ind[i]);
		//printf("\n");
		cb->func(cb->array, cb->res, ind);
  }
  else {
    for (i=0;i<list[0];i++) {
      ind[pos] = i;
      __rec_indices(ind, &list[1], size, pos+1, c);
    }
  }
}

void call_on_indices(int* list, int size, void*c) {
	int i;
	callback* cb;
	int* indices;

	indices = (int*)malloc(sizeof(int)*size);
	memset(indices, -1, sizeof(int)*size);
	
	// Special case
	if (size == 1) {
		indices[0] = list[0];
		cb = (callback*) c;
		cb->func(cb->array, cb->res, indices);
	}
	// Otherwise calculate the indices recursively
	// Unfortunately this is tree-recursive and can't be
	// converted into a tail recursive version
	for (i=0;i<list[0];i++) {
		indices[0] = i;
		__rec_indices(indices, &list[1], size, 1, c);
	}
	free(indices);
}

static PyObject *
length(PyObject *self, PyObject *args)
{
  PyArrayObject *array = NULL, *res = NULL;
  int *dim;
  int nd,i;
  callback cb;
	
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

	cb.func = __length;
	cb.array = array;
	cb.res = res;
	call_on_indices(array->dimensions, array->nd, &cb);
	
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
	int nd,i;
  callback cb;
	
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

	cb.func = __norm;
	cb.array = array;
	cb.res = res;
	call_on_indices(array->dimensions, array->nd, &cb);

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
	{NULL, NULL, 0, NULL}	/* Sentinel */
};

PyMODINIT_FUNC
init_numtools(void)
{
	(void) Py_InitModule("_numtools", NumericMethods);
	import_array();
}

