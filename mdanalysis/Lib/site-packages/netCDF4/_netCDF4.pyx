"""
Version 1.6.3
-------------

# Introduction

netcdf4-python is a Python interface to the netCDF C library.

[netCDF](http://www.unidata.ucar.edu/software/netcdf/) version 4 has many features
not found in earlier versions of the library and is implemented on top of
[HDF5](http://www.hdfgroup.org/HDF5). This module can read and write
files in both the new netCDF 4 and the old netCDF 3 format, and can create
files that are readable by HDF5 clients. The API modelled after
[Scientific.IO.NetCDF](http://dirac.cnrs-orleans.fr/ScientificPython/),
and should be familiar to users of that module.

Most new features of netCDF 4 are implemented, such as multiple
unlimited dimensions, groups and data compression.  All the new
numeric data types (such as 64 bit and unsigned integer types) are
implemented. Compound (struct), variable length (vlen) and
enumerated (enum) data types are supported, but not the opaque data type.
Mixtures of compound, vlen and enum data types (such as
compound types containing enums, or vlens containing compound
types) are not supported.

## Quick Install

 - the easiest way to get going is to install via `pip install netCDF4`.
   (or if you use the [conda](http://conda.io) package manager `conda install -c conda-forge netCDF4`).

## Developer Install

 - Clone the
   [github repository](http://github.com/Unidata/netcdf4-python).
 - Make sure the dependencies are satisfied (Python 3.7 or later,
   [numpy](http://numpy.scipy.org), 
   [Cython](http://cython.org),
   [cftime](https://github.com/Unidata/cftime),
   [setuptools](https://pypi.python.org/pypi/setuptools),
   the [HDF5 C library](https://www.hdfgroup.org/solutions/hdf5/),
   and the [netCDF C library](https://www.unidata.ucar.edu/software/netcdf/)).
   For MPI parallel IO support, an MPI-enabled versions of the netcdf library
   is required, as is [mpi4py](http://mpi4py.scipy.org).
   Parallel IO further depends on the existence of MPI-enabled HDF5 or the
   [PnetCDF](https://parallel-netcdf.github.io/) library.
 - By default, the utility `nc-config` (installed with netcdf-c)
   will be run used to determine where all the dependencies live.
 - If `nc-config` is not in your default `PATH`, you can set the `NETCDF4_DIR`
   environment variable and `setup.py` will look in `$NETCDF4_DIR/bin`.
   You can also use the file `setup.cfg` to set the path to `nc-config`, or
   enter the paths to the libraries and include files manually. Just  edit the `setup.cfg` file
   in a text editor and follow the instructions in the comments.
   To disable the use of `nc-config`, set the env var `USE_NCCONFIG` to 0.
   To disable the use of `setup.cfg`, set `USE_SETUPCFG` to 0.
   As a last resort, the library and include paths can be set via environment variables.
   If you go this route, set `USE_NCCONFIG` and `USE_SETUPCFG` to 0, and specify
   `NETCDF4_LIBDIR`, `NETCDF4_INCDIR`, `HDF5_LIBDIR` and `HDF5_INCDIR`.
   If the dependencies are not found
   in any of the paths specified by environment variables, then standard locations
   (such as `/usr` and `/usr/local`) are searched.
 - if the env var `NETCDF_PLUGIN_DIR` is set to point to the location of the netcdf-c compression
   plugins built by netcdf >= 4.9.0, they will be installed inside the package.  In this
   case `HDF5_PLUGIN_PATH` will be set to the package installation path on import,
   so the extra compression algorithms available in netcdf-c >= 4.9.0 will automatically
   be available.  Otherwise, the user will have to set `HDF5_PLUGIN_PATH` explicitly
   to have access to the extra compression plugins.
 - run `python setup.py build`, then `python setup.py install` (as root if
   necessary).
 - run the tests in the 'test' directory by running `python run_all.py`.

# Tutorial

- [Creating/Opening/Closing a netCDF file](#creatingopeningclosing-a-netcdf-file)
- [Groups in a netCDF file](#groups-in-a-netcdf-file)
- [Dimensions in a netCDF file](#dimensions-in-a-netcdf-file)
- [Variables in a netCDF file](#variables-in-a-netcdf-file)
- [Attributes in a netCDF file](#attributes-in-a-netcdf-file)
- [Dealing with time coordinates](#dealing-with-time-coordinates)
- [Writing data to and retrieving data from a netCDF variable](#writing-data-to-and-retrieving-data-from-a-netcdf-variable)
- [Reading data from a multi-file netCDF dataset](#reading-data-from-a-multi-file-netcdf-dataset)
- [Efficient compression of netCDF variables](#efficient-compression-of-netcdf-variables)
- [Beyond homogeneous arrays of a fixed type - compound data types](#beyond-homogeneous-arrays-of-a-fixed-type-compound-data-types)
- [Variable-length (vlen) data types](#variable-length-vlen-data-types)
- [Enum data type](#enum-data-type)
- [Parallel IO](#parallel-io)
- [Dealing with strings](#dealing-with-strings)
- [In-memory (diskless) Datasets](#in-memory-diskless-datasets)

## Creating/Opening/Closing a netCDF file

To create a netCDF file from python, you simply call the `Dataset`
constructor. This is also the method used to open an existing netCDF
file.  If the file is open for write access (`mode='w', 'r+'` or `'a'`), you may
write any type of data including new dimensions, groups, variables and
attributes.  netCDF files come in five flavors (`NETCDF3_CLASSIC,
NETCDF3_64BIT_OFFSET, NETCDF3_64BIT_DATA, NETCDF4_CLASSIC`, and `NETCDF4`).
`NETCDF3_CLASSIC` was the original netcdf binary format, and was limited
to file sizes less than 2 Gb. `NETCDF3_64BIT_OFFSET` was introduced
in version 3.6.0 of the library, and extended the original binary format
to allow for file sizes greater than 2 Gb.
`NETCDF3_64BIT_DATA` is a new format that requires version 4.4.0 of
the C library - it extends the `NETCDF3_64BIT_OFFSET` binary format to
allow for unsigned/64 bit integer data types and 64-bit dimension sizes.
`NETCDF3_64BIT` is an alias for `NETCDF3_64BIT_OFFSET`.
`NETCDF4_CLASSIC` files use the version 4 disk format (HDF5), but omits features
not found in the version 3 API. They can be read by netCDF 3 clients
only if they have been relinked against the netCDF 4 library. They can
also be read by HDF5 clients. `NETCDF4` files use the version 4 disk
format (HDF5) and use the new features of the version 4 API.  The
netCDF4 module can read and write files in any of these formats. When
creating a new file, the format may be specified using the `format`
keyword in the `Dataset` constructor.  The default format is
`NETCDF4`. To see how a given file is formatted, you can examine the
`data_model` attribute.  Closing the netCDF file is
accomplished via the `Dataset.close` method of the `Dataset`
instance.

Here's an example:

```python
>>> from netCDF4 import Dataset
>>> rootgrp = Dataset("test.nc", "w", format="NETCDF4")
>>> print(rootgrp.data_model)
NETCDF4
>>> rootgrp.close()
```

Remote [OPeNDAP](http://opendap.org)-hosted datasets can be accessed for
reading over http if a URL is provided to the `Dataset` constructor instead of a
filename.  However, this requires that the netCDF library be built with
OPenDAP support, via the `--enable-dap` configure option (added in
version 4.0.1).


## Groups in a netCDF file

netCDF version 4 added support for organizing data in hierarchical
groups, which are analogous to directories in a filesystem. Groups serve
as containers for variables, dimensions and attributes, as well as other
groups.  A `Dataset` creates a special group, called
the 'root group', which is similar to the root directory in a unix
filesystem.  To create `Group` instances, use the
`Dataset.createGroup` method of a `Dataset` or `Group`
instance. `Dataset.createGroup` takes a single argument, a
python string containing the name of the new group. The new `Group`
instances contained within the root group can be accessed by name using
the `groups` dictionary attribute of the `Dataset` instance.  Only
`NETCDF4` formatted files support Groups, if you try to create a Group
in a netCDF 3 file you will get an error message.

```python
>>> rootgrp = Dataset("test.nc", "a")
>>> fcstgrp = rootgrp.createGroup("forecasts")
>>> analgrp = rootgrp.createGroup("analyses")
>>> print(rootgrp.groups)
{'forecasts': <class 'netCDF4._netCDF4.Group'>
group /forecasts:
    dimensions(sizes): 
    variables(dimensions): 
    groups: , 'analyses': <class 'netCDF4._netCDF4.Group'>
group /analyses:
    dimensions(sizes): 
    variables(dimensions): 
    groups: }
>>>
```

Groups can exist within groups in a `Dataset`, just as directories
exist within directories in a unix filesystem. Each `Group` instance
has a `groups` attribute dictionary containing all of the group
instances contained within that group. Each `Group` instance also has a
`path` attribute that contains a simulated unix directory path to
that group.  To simplify the creation of nested groups, you can
use a unix-like path as an argument to `Dataset.createGroup`.

```python
>>> fcstgrp1 = rootgrp.createGroup("/forecasts/model1")
>>> fcstgrp2 = rootgrp.createGroup("/forecasts/model2")
```

If any of the intermediate elements of the path do not exist, they are created,
just as with the unix command `'mkdir -p'`. If you try to create a group
that already exists, no error will be raised, and the existing group will be
returned.

Here's an example that shows how to navigate all the groups in a
`Dataset`. The function `walktree` is a Python generator that is used
to walk the directory tree. Note that printing the `Dataset` or `Group`
object yields summary information about it's contents.

```python
>>> def walktree(top):
...     yield top.groups.values()
...     for value in top.groups.values():
...         yield from walktree(value)
>>> print(rootgrp)
<class 'netCDF4._netCDF4.Dataset'>
root group (NETCDF4 data model, file format HDF5):
    dimensions(sizes): 
    variables(dimensions): 
    groups: forecasts, analyses
>>> for children in walktree(rootgrp):
...     for child in children:
...         print(child)
<class 'netCDF4._netCDF4.Group'>
group /forecasts:
    dimensions(sizes): 
    variables(dimensions): 
    groups: model1, model2
<class 'netCDF4._netCDF4.Group'>
group /analyses:
    dimensions(sizes): 
    variables(dimensions): 
    groups: 
<class 'netCDF4._netCDF4.Group'>
group /forecasts/model1:
    dimensions(sizes): 
    variables(dimensions): 
    groups: 
<class 'netCDF4._netCDF4.Group'>
group /forecasts/model2:
    dimensions(sizes): 
    variables(dimensions): 
    groups: 
```

## Dimensions in a netCDF file

netCDF defines the sizes of all variables in terms of dimensions, so
before any variables can be created the dimensions they use must be
created first. A special case, not often used in practice, is that of a
scalar variable, which has no dimensions. A dimension is created using
the `Dataset.createDimension` method of a `Dataset`
or `Group` instance. A Python string is used to set the name of the
dimension, and an integer value is used to set the size. To create an
unlimited dimension (a dimension that can be appended to), the size
value is set to `None` or 0. In this example, there both the `time` and
`level` dimensions are unlimited.  Having more than one unlimited
dimension is a new netCDF 4 feature, in netCDF 3 files there may be only
one, and it must be the first (leftmost) dimension of the variable.

```python
>>> level = rootgrp.createDimension("level", None)
>>> time = rootgrp.createDimension("time", None)
>>> lat = rootgrp.createDimension("lat", 73)
>>> lon = rootgrp.createDimension("lon", 144)
```


All of the `Dimension` instances are stored in a python dictionary.

```python
>>> print(rootgrp.dimensions)
{'level': <class 'netCDF4._netCDF4.Dimension'> (unlimited): name = 'level', size = 0, 'time': <class 'netCDF4._netCDF4.Dimension'> (unlimited): name = 'time', size = 0, 'lat': <class 'netCDF4._netCDF4.Dimension'>: name = 'lat', size = 73, 'lon': <class 'netCDF4._netCDF4.Dimension'>: name = 'lon', size = 144}
```

Using the python `len` function with a `Dimension` instance returns
current size of that dimension.
`Dimension.isunlimited` method of a `Dimension` instance
be used to determine if the dimensions is unlimited, or appendable.

```python
>>> print(len(lon))
144
>>> print(lon.isunlimited())
False
>>> print(time.isunlimited())
True
```

Printing the `Dimension` object
provides useful summary info, including the name and length of the dimension,
and whether it is unlimited.

```python
>>> for dimobj in rootgrp.dimensions.values():
...     print(dimobj)
<class 'netCDF4._netCDF4.Dimension'> (unlimited): name = 'level', size = 0
<class 'netCDF4._netCDF4.Dimension'> (unlimited): name = 'time', size = 0
<class 'netCDF4._netCDF4.Dimension'>: name = 'lat', size = 73
<class 'netCDF4._netCDF4.Dimension'>: name = 'lon', size = 144
```

`Dimension` names can be changed using the
`Datatset.renameDimension` method of a `Dataset` or
`Group` instance.

## Variables in a netCDF file

netCDF variables behave much like python multidimensional array objects
supplied by the [numpy module](http://numpy.scipy.org). However,
unlike numpy arrays, netCDF4 variables can be appended to along one or
more 'unlimited' dimensions. To create a netCDF variable, use the
`Dataset.createVariable` method of a `Dataset` or
`Group` instance. The `Dataset.createVariable`j method
has two mandatory arguments, the variable name (a Python string), and
the variable datatype. The variable's dimensions are given by a tuple
containing the dimension names (defined previously with
`Dataset.createDimension`). To create a scalar
variable, simply leave out the dimensions keyword. The variable
primitive datatypes correspond to the dtype attribute of a numpy array.
You can specify the datatype as a numpy dtype object, or anything that
can be converted to a numpy dtype object.  Valid datatype specifiers
include: `'f4'` (32-bit floating point), `'f8'` (64-bit floating
point), `'i4'` (32-bit signed integer), `'i2'` (16-bit signed
integer), `'i8'` (64-bit signed integer), `'i1'` (8-bit signed
integer), `'u1'` (8-bit unsigned integer), `'u2'` (16-bit unsigned
integer), `'u4'` (32-bit unsigned integer), `'u8'` (64-bit unsigned
integer), or `'S1'` (single-character string).  The old Numeric
single-character typecodes (`'f'`,`'d'`,`'h'`,
`'s'`,`'b'`,`'B'`,`'c'`,`'i'`,`'l'`), corresponding to
(`'f4'`,`'f8'`,`'i2'`,`'i2'`,`'i1'`,`'i1'`,`'S1'`,`'i4'`,`'i4'`),
will also work. The unsigned integer types and the 64-bit integer type
can only be used if the file format is `NETCDF4`.

The dimensions themselves are usually also defined as variables, called
coordinate variables. The `Dataset.createVariable`
method returns an instance of the `Variable` class whose methods can be
used later to access and set variable data and attributes.

```python
>>> times = rootgrp.createVariable("time","f8",("time",))
>>> levels = rootgrp.createVariable("level","i4",("level",))
>>> latitudes = rootgrp.createVariable("lat","f4",("lat",))
>>> longitudes = rootgrp.createVariable("lon","f4",("lon",))
>>> # two dimensions unlimited
>>> temp = rootgrp.createVariable("temp","f4",("time","level","lat","lon",))
>>> temp.units = "K"
```

To get summary info on a `Variable` instance in an interactive session,
just print it.

```python
>>> print(temp)
<class 'netCDF4._netCDF4.Variable'>
float32 temp(time, level, lat, lon)
    units: K
unlimited dimensions: time, level
current shape = (0, 0, 73, 144)
filling on, default _FillValue of 9.969209968386869e+36 used
```

You can use a path to create a Variable inside a hierarchy of groups.

```python
>>> ftemp = rootgrp.createVariable("/forecasts/model1/temp","f4",("time","level","lat","lon",))
```

If the intermediate groups do not yet exist, they will be created.

You can also query a `Dataset` or `Group` instance directly to obtain `Group` or
`Variable` instances using paths.

```python
>>> print(rootgrp["/forecasts/model1"])  # a Group instance
<class 'netCDF4._netCDF4.Group'>
group /forecasts/model1:
    dimensions(sizes): 
    variables(dimensions): float32 temp(time,level,lat,lon)
    groups: 
>>> print(rootgrp["/forecasts/model1/temp"])  # a Variable instance
<class 'netCDF4._netCDF4.Variable'>
float32 temp(time, level, lat, lon)
path = /forecasts/model1
unlimited dimensions: time, level
current shape = (0, 0, 73, 144)
filling on, default _FillValue of 9.969209968386869e+36 used
```


All of the variables in the `Dataset` or `Group` are stored in a
Python dictionary, in the same way as the dimensions:

```python
>>> print(rootgrp.variables)
{'time': <class 'netCDF4._netCDF4.Variable'>
float64 time(time)
unlimited dimensions: time
current shape = (0,)
filling on, default _FillValue of 9.969209968386869e+36 used, 'level': <class 'netCDF4._netCDF4.Variable'>
int32 level(level)
unlimited dimensions: level
current shape = (0,)
filling on, default _FillValue of -2147483647 used, 'lat': <class 'netCDF4._netCDF4.Variable'>
float32 lat(lat)
unlimited dimensions: 
current shape = (73,)
filling on, default _FillValue of 9.969209968386869e+36 used, 'lon': <class 'netCDF4._netCDF4.Variable'>
float32 lon(lon)
unlimited dimensions: 
current shape = (144,)
filling on, default _FillValue of 9.969209968386869e+36 used, 'temp': <class 'netCDF4._netCDF4.Variable'>
float32 temp(time, level, lat, lon)
    units: K
unlimited dimensions: time, level
current shape = (0, 0, 73, 144)
filling on, default _FillValue of 9.969209968386869e+36 used}
```

`Variable` names can be changed using the
`Dataset.renameVariable` method of a `Dataset`
instance.

Variables can be sliced similar to numpy arrays, but there are some differences.  See
[Writing data to and retrieving data from a netCDF variable](#writing-data-to-and-retrieving-data-from-a-netcdf-variable) for more details.


## Attributes in a netCDF file

There are two types of attributes in a netCDF file, global and variable.
Global attributes provide information about a group, or the entire
dataset, as a whole. `Variable` attributes provide information about
one of the variables in a group. Global attributes are set by assigning
values to `Dataset` or `Group` instance variables. `Variable`
attributes are set by assigning values to `Variable` instances
variables. Attributes can be strings, numbers or sequences. Returning to
our example,

```python
>>> import time
>>> rootgrp.description = "bogus example script"
>>> rootgrp.history = "Created " + time.ctime(time.time())
>>> rootgrp.source = "netCDF4 python module tutorial"
>>> latitudes.units = "degrees north"
>>> longitudes.units = "degrees east"
>>> levels.units = "hPa"
>>> temp.units = "K"
>>> times.units = "hours since 0001-01-01 00:00:00.0"
>>> times.calendar = "gregorian"
```

The `Dataset.ncattrs` method of a `Dataset`, `Group` or
`Variable` instance can be used to retrieve the names of all the netCDF
attributes. This method is provided as a convenience, since using the
built-in `dir` Python function will return a bunch of private methods
and attributes that cannot (or should not) be modified by the user.

```python
>>> for name in rootgrp.ncattrs():
...     print("Global attr {} = {}".format(name, getattr(rootgrp, name)))
Global attr description = bogus example script
Global attr history = Created Mon Jul  8 14:19:41 2019
Global attr source = netCDF4 python module tutorial
```

The `__dict__` attribute of a `Dataset`, `Group` or `Variable`
instance provides all the netCDF attribute name/value pairs in a python
dictionary:

```python
>>> print(rootgrp.__dict__)
{'description': 'bogus example script', 'history': 'Created Mon Jul  8 14:19:41 2019', 'source': 'netCDF4 python module tutorial'}
```

Attributes can be deleted from a netCDF `Dataset`, `Group` or
`Variable` using the python `del` statement (i.e. `del grp.foo`
removes the attribute `foo` the the group `grp`).

## Writing data to and retrieving data from a netCDF variable

Now that you have a netCDF `Variable` instance, how do you put data
into it? You can just treat it like an array and assign data to a slice.

```python
>>> import numpy as np
>>> lats =  np.arange(-90,91,2.5)
>>> lons =  np.arange(-180,180,2.5)
>>> latitudes[:] = lats
>>> longitudes[:] = lons
>>> print("latitudes =\\n{}".format(latitudes[:]))
latitudes =
[-90.  -87.5 -85.  -82.5 -80.  -77.5 -75.  -72.5 -70.  -67.5 -65.  -62.5
 -60.  -57.5 -55.  -52.5 -50.  -47.5 -45.  -42.5 -40.  -37.5 -35.  -32.5
 -30.  -27.5 -25.  -22.5 -20.  -17.5 -15.  -12.5 -10.   -7.5  -5.   -2.5
   0.    2.5   5.    7.5  10.   12.5  15.   17.5  20.   22.5  25.   27.5
  30.   32.5  35.   37.5  40.   42.5  45.   47.5  50.   52.5  55.   57.5
  60.   62.5  65.   67.5  70.   72.5  75.   77.5  80.   82.5  85.   87.5
  90. ]
```

Unlike NumPy's array objects, netCDF `Variable`
objects with unlimited dimensions will grow along those dimensions if you
assign data outside the currently defined range of indices.

```python
>>> # append along two unlimited dimensions by assigning to slice.
>>> nlats = len(rootgrp.dimensions["lat"])
>>> nlons = len(rootgrp.dimensions["lon"])
>>> print("temp shape before adding data = {}".format(temp.shape))
temp shape before adding data = (0, 0, 73, 144)
>>>
>>> from numpy.random import uniform
>>> temp[0:5, 0:10, :, :] = uniform(size=(5, 10, nlats, nlons))
>>> print("temp shape after adding data = {}".format(temp.shape))
temp shape after adding data = (5, 10, 73, 144)
>>>
>>> # levels have grown, but no values yet assigned.
>>> print("levels shape after adding pressure data = {}".format(levels.shape))
levels shape after adding pressure data = (10,)
```

Note that the size of the levels variable grows when data is appended
along the `level` dimension of the variable `temp`, even though no
data has yet been assigned to levels.

```python
>>> # now, assign data to levels dimension variable.
>>> levels[:] =  [1000.,850.,700.,500.,300.,250.,200.,150.,100.,50.]
```

However, that there are some differences between NumPy and netCDF
variable slicing rules. Slices behave as usual, being specified as a
`start:stop:step` triplet. Using a scalar integer index `i` takes the ith
element and reduces the rank of the output array by one. Boolean array and
integer sequence indexing behaves differently for netCDF variables
than for numpy arrays.  Only 1-d boolean arrays and integer sequences are
allowed, and these indices work independently along each dimension (similar
to the way vector subscripts work in fortran).  This means that

```python
>>> temp[0, 0, [0,1,2,3], [0,1,2,3]].shape
(4, 4)
```

returns an array of shape (4,4) when slicing a netCDF variable, but for a
numpy array it returns an array of shape (4,).
Similarly, a netCDF variable of shape `(2,3,4,5)` indexed
with `[0, array([True, False, True]), array([False, True, True, True]), :]`
would return a `(2, 3, 5)` array. In NumPy, this would raise an error since
it would be equivalent to `[0, [0,1], [1,2,3], :]`. When slicing with integer
sequences, the indices ***need not be sorted*** and ***may contain
duplicates*** (both of these are new features in version 1.2.1).
While this behaviour may cause some confusion for those used to NumPy's 'fancy indexing' rules,
it provides a very powerful way to extract data from multidimensional netCDF
variables by using logical operations on the dimension arrays to create slices.

For example,

```python
>>> tempdat = temp[::2, [1,3,6], lats>0, lons>0]
```

will extract time indices 0,2 and 4, pressure levels
850, 500 and 200 hPa, all Northern Hemisphere latitudes and Eastern
Hemisphere longitudes, resulting in a numpy array of shape  (3, 3, 36, 71).

```python
>>> print("shape of fancy temp slice = {}".format(tempdat.shape))
shape of fancy temp slice = (3, 3, 36, 71)
```

***Special note for scalar variables***: To extract data from a scalar variable
`v` with no associated dimensions, use `numpy.asarray(v)` or `v[...]`.
The result will be a numpy scalar array.

By default, netcdf4-python returns numpy masked arrays with values equal to the
`missing_value` or `_FillValue` variable attributes masked for primitive and
enum data types.
The `Dataset.set_auto_mask` `Dataset` and `Variable` methods
can be used to disable this feature so that
numpy arrays are always returned, with the missing values included. Prior to
version 1.4.0 the default behavior was to only return masked arrays when the
requested slice contained missing values.  This behavior can be recovered
using the `Dataset.set_always_mask` method. If a masked array is
written to a netCDF variable, the masked elements are filled with the
value specified by the `missing_value` attribute.  If the variable has
no `missing_value`, the `_FillValue` is used instead.

## Dealing with time coordinates

Time coordinate values pose a special challenge to netCDF users.  Most
metadata standards (such as CF) specify that time should be
measure relative to a fixed date using a certain calendar, with units
specified like `hours since YY-MM-DD hh:mm:ss`.  These units can be
awkward to deal with, without a utility to convert the values to and
from calendar dates.  The functions [num2date](https://unidata.github.io/cftime/api.html)
and [date2num](https://unidata.github.io/cftime/api.html) are
provided by [cftime](https://unidata.github.io/cftime) to do just that.
Here's an example of how they can be used:

```python
>>> # fill in times.
>>> from datetime import datetime, timedelta
>>> from cftime import num2date, date2num
>>> dates = [datetime(2001,3,1)+n*timedelta(hours=12) for n in range(temp.shape[0])]
>>> times[:] = date2num(dates,units=times.units,calendar=times.calendar)
>>> print("time values (in units {}):\\n{}".format(times.units, times[:]))
time values (in units hours since 0001-01-01 00:00:00.0):
[17533104. 17533116. 17533128. 17533140. 17533152.]
>>> dates = num2date(times[:],units=times.units,calendar=times.calendar)
>>> print("dates corresponding to time values:\\n{}".format(dates))
 [cftime.DatetimeGregorian(2001, 3, 1, 0, 0, 0, 0, has_year_zero=False)
  cftime.DatetimeGregorian(2001, 3, 1, 12, 0, 0, 0, has_year_zero=False)
  cftime.DatetimeGregorian(2001, 3, 2, 0, 0, 0, 0, has_year_zero=False)
  cftime.DatetimeGregorian(2001, 3, 2, 12, 0, 0, 0, has_year_zero=False)
  cftime.DatetimeGregorian(2001, 3, 3, 0, 0, 0, 0, has_year_zero=False)]
```

`num2date` converts numeric values of time in the specified `units`
and `calendar` to datetime objects, and `date2num` does the reverse.
All the calendars currently defined in the
[CF metadata convention](http://cfconventions.org) are supported.
A function called `date2index` is also provided which returns the indices
of a netCDF time variable corresponding to a sequence of datetime instances.


## Reading data from a multi-file netCDF dataset

If you want to read data from a variable that spans multiple netCDF files,
you can use the `MFDataset` class to read the data as if it were
contained in a single file. Instead of using a single filename to create
a `Dataset` instance, create a `MFDataset` instance with either a list
of filenames, or a string with a wildcard (which is then converted to
a sorted list of files using the python glob module).
Variables in the list of files that share the same unlimited
dimension are aggregated together, and can be sliced across multiple
files.  To illustrate this, let's first create a bunch of netCDF files with
the same variable (with the same unlimited dimension).  The files
must in be in `NETCDF3_64BIT_OFFSET`, `NETCDF3_64BIT_DATA`, `NETCDF3_CLASSIC` or
`NETCDF4_CLASSIC` format (`NETCDF4` formatted multi-file
datasets are not supported).

```python
>>> for nf in range(10):
...     with Dataset("mftest%s.nc" % nf, "w", format="NETCDF4_CLASSIC") as f:
...         _ = f.createDimension("x",None)
...         x = f.createVariable("x","i",("x",))
...         x[0:10] = np.arange(nf*10,10*(nf+1))
```

Now read all the files back in at once with `MFDataset`

```python
>>> from netCDF4 import MFDataset
>>> f = MFDataset("mftest*nc")
>>> print(f.variables["x"][:])
[ 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47
 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71
 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95
 96 97 98 99]
```

Note that `MFDataset` can only be used to read, not write, multi-file
datasets.

## Efficient compression of netCDF variables

Data stored in netCDF `Variable` objects can be compressed and
decompressed on the fly. The compression algorithm used is determined
by the `compression` keyword argument to the `Dataset.createVariable` method.
`zlib` compression is always available, `szip` is available if the linked HDF5
library supports it, and `zstd`, `bzip2`, `blosc_lz`,`blosc_lz4`,`blosc_lz4hc`,
`blosc_zlib` and `blosc_zstd` are available via optional external plugins.
The `complevel` keyword regulates the
speed and efficiency of the compression for `zlib`, `bzip` and `zstd` (1 being fastest, but lowest
compression ratio, 9 being slowest but best compression ratio). The
default value of `complevel` is 4. Setting `shuffle=False` will turn
off the HDF5 shuffle filter, which de-interlaces a block of data before
`zlib` compression by reordering the bytes.  The shuffle filter can
significantly improve compression ratios, and is on by default if `compression=zlib`.  Setting
`fletcher32` keyword argument to
`Dataset.createVariable` to `True` (it's `False` by
default) enables the Fletcher32 checksum algorithm for error detection.
It's also possible to set the HDF5 chunking parameters and endian-ness
of the binary data stored in the HDF5 file with the `chunksizes`
and `endian` keyword arguments to
`Dataset.createVariable`.  These keyword arguments only
are relevant for `NETCDF4` and `NETCDF4_CLASSIC` files (where the
underlying file format is HDF5) and are silently ignored if the file
format is `NETCDF3_CLASSIC`, `NETCDF3_64BIT_OFFSET` or `NETCDF3_64BIT_DATA`.
If the HDF5 library is built with szip support, compression=`szip` can also
be used (in conjunction with the `szip_coding` and `szip_pixels_per_block` keyword
arguments).  

If your data only has a certain number of digits of precision (say for
example, it is temperature data that was measured with a precision of
0.1 degrees), you can dramatically improve compression by
quantizing (or truncating) the data. There are two methods supplied for
doing this.  You can use the `least_significant_digit`
keyword argument to `Dataset.createVariable` to specify
the power of ten of the smallest decimal place in
the data that is a reliable value. For example if the data has a
precision of 0.1, then setting `least_significant_digit=1` will cause
data the data to be quantized using `numpy.around(scale*data)/scale`, where
scale = 2**bits, and bits is determined so that a precision of 0.1 is
retained (in this case bits=4).  This is done at the python level and is
not a part of the underlying C library.  Starting with netcdf-c version 4.9.0,
a quantization capability is provided in the library.  This can be
used via the `significant_digits` `Dataset.createVariable` kwarg (new in
version 1.6.0).
The interpretation of `significant_digits` is different than `least_signficant_digit`
in that it specifies the absolute number of significant digits independent
of the magnitude of the variable (the floating point exponent).
Either of these approaches makes the compression
'lossy' instead of 'lossless', that is some precision in the data is
sacrificed for the sake of disk space.

In our example, try replacing the line

```python
>>> temp = rootgrp.createVariable("temp","f4",("time","level","lat","lon",))
```

with

```python
>>> temp = rootgrp.createVariable("temp","f4",("time","level","lat","lon",),compression='zlib')
```

and then

```python
>>> temp = rootgrp.createVariable("temp","f4",("time","level","lat","lon",),compression='zlib',least_significant_digit=3)
```

or with netcdf-c >= 4.9.0

```python
>>> temp = rootgrp.createVariable("temp","f4",("time","level","lat","lon",),compression='zlib',significant_digits=4)
```

and see how much smaller the resulting files are.

## Beyond homogeneous arrays of a fixed type - compound data types

Compound data types map directly to numpy structured (a.k.a 'record')
arrays.  Structured arrays are akin to C structs, or derived types
in Fortran. They allow for the construction of table-like structures
composed of combinations of other data types, including other
compound types. Compound types might be useful for representing multiple
parameter values at each point on a grid, or at each time and space
location for scattered (point) data. You can then access all the
information for a point by reading one variable, instead of reading
different parameters from different variables.  Compound data types
are created from the corresponding numpy data type using the
`Dataset.createCompoundType` method of a `Dataset` or `Group` instance.
Since there is no native complex data type in netcdf, compound types are handy
for storing numpy complex arrays.  Here's an example:

```python
>>> f = Dataset("complex.nc","w")
>>> size = 3 # length of 1-d complex array
>>> # create sample complex data.
>>> datac = np.exp(1j*(1.+np.linspace(0, np.pi, size)))
>>> # create complex128 compound data type.
>>> complex128 = np.dtype([("real",np.float64),("imag",np.float64)])
>>> complex128_t = f.createCompoundType(complex128,"complex128")
>>> # create a variable with this data type, write some data to it.
>>> x_dim = f.createDimension("x_dim",None)
>>> v = f.createVariable("cmplx_var",complex128_t,"x_dim")
>>> data = np.empty(size,complex128) # numpy structured array
>>> data["real"] = datac.real; data["imag"] = datac.imag
>>> v[:] = data # write numpy structured array to netcdf compound var
>>> # close and reopen the file, check the contents.
>>> f.close(); f = Dataset("complex.nc")
>>> v = f.variables["cmplx_var"]
>>> datain = v[:] # read in all the data into a numpy structured array
>>> # create an empty numpy complex array
>>> datac2 = np.empty(datain.shape,np.complex128)
>>> # .. fill it with contents of structured array.
>>> datac2.real = datain["real"]; datac2.imag = datain["imag"]
>>> print('{}: {}'.format(datac.dtype, datac)) # original data
complex128: [ 0.54030231+0.84147098j -0.84147098+0.54030231j -0.54030231-0.84147098j]
>>>
>>> print('{}: {}'.format(datac2.dtype, datac2)) # data from file
complex128: [ 0.54030231+0.84147098j -0.84147098+0.54030231j -0.54030231-0.84147098j]
```

Compound types can be nested, but you must create the 'inner'
ones first. All possible numpy structured arrays cannot be
represented as Compound variables - an error message will be
raise if you try to create one that is not supported.
All of the compound types defined for a `Dataset` or `Group` are stored
in a Python dictionary, just like variables and dimensions. As always, printing
objects gives useful summary information in an interactive session:

```python
>>> print(f)
<class 'netCDF4._netCDF4.Dataset'>
root group (NETCDF4 data model, file format HDF5):
    dimensions(sizes): x_dim(3)
    variables(dimensions): {'names':['real','imag'], 'formats':['<f8','<f8'], 'offsets':[0,8], 'itemsize':16, 'aligned':True} cmplx_var(x_dim)
    groups: 
>>> print(f.variables["cmplx_var"])
<class 'netCDF4._netCDF4.Variable'>
compound cmplx_var(x_dim)
compound data type: {'names':['real','imag'], 'formats':['<f8','<f8'], 'offsets':[0,8], 'itemsize':16, 'aligned':True}
unlimited dimensions: x_dim
current shape = (3,)
>>> print(f.cmptypes)
{'complex128': <class 'netCDF4._netCDF4.CompoundType'>: name = 'complex128', numpy dtype = {'names':['real','imag'], 'formats':['<f8','<f8'], 'offsets':[0,8], 'itemsize':16, 'aligned':True}}
>>> print(f.cmptypes["complex128"])
<class 'netCDF4._netCDF4.CompoundType'>: name = 'complex128', numpy dtype = {'names':['real','imag'], 'formats':['<f8','<f8'], 'offsets':[0,8], 'itemsize':16, 'aligned':True}
```

## Variable-length (vlen) data types

NetCDF 4 has support for variable-length or "ragged" arrays.  These are arrays
of variable length sequences having the same type. To create a variable-length
data type, use the `Dataset.createVLType` method
method of a `Dataset` or `Group` instance.

```python
>>> f = Dataset("tst_vlen.nc","w")
>>> vlen_t = f.createVLType(np.int32, "phony_vlen")
```

The numpy datatype of the variable-length sequences and the name of the
new datatype must be specified. Any of the primitive datatypes can be
used (signed and unsigned integers, 32 and 64 bit floats, and characters),
but compound data types cannot.
A new variable can then be created using this datatype.

```python
>>> x = f.createDimension("x",3)
>>> y = f.createDimension("y",4)
>>> vlvar = f.createVariable("phony_vlen_var", vlen_t, ("y","x"))
```

Since there is no native vlen datatype in numpy, vlen arrays are represented
in python as object arrays (arrays of dtype `object`). These are arrays whose
elements are Python object pointers, and can contain any type of python object.
For this application, they must contain 1-D numpy arrays all of the same type
but of varying length.
In this case, they contain 1-D numpy `int32` arrays of random length between
1 and 10.

```python
>>> import random
>>> random.seed(54321)
>>> data = np.empty(len(y)*len(x),object)
>>> for n in range(len(y)*len(x)):
...     data[n] = np.arange(random.randint(1,10),dtype="int32")+1
>>> data = np.reshape(data,(len(y),len(x)))
>>> vlvar[:] = data
>>> print("vlen variable =\\n{}".format(vlvar[:]))
vlen variable =
[[array([1, 2, 3, 4, 5, 6, 7, 8], dtype=int32) array([1, 2], dtype=int32)
  array([1, 2, 3, 4], dtype=int32)]
 [array([1, 2, 3], dtype=int32)
  array([1, 2, 3, 4, 5, 6, 7, 8, 9], dtype=int32)
  array([1, 2, 3, 4, 5, 6, 7, 8, 9], dtype=int32)]
 [array([1, 2, 3, 4, 5, 6, 7], dtype=int32) array([1, 2, 3], dtype=int32)
  array([1, 2, 3, 4, 5, 6], dtype=int32)]
 [array([1, 2, 3, 4, 5, 6, 7, 8, 9], dtype=int32)
  array([1, 2, 3, 4, 5], dtype=int32) array([1, 2], dtype=int32)]]
>>> print(f)
<class 'netCDF4._netCDF4.Dataset'>
root group (NETCDF4 data model, file format HDF5):
    dimensions(sizes): x(3), y(4)
    variables(dimensions): int32 phony_vlen_var(y,x)
    groups: 
>>> print(f.variables["phony_vlen_var"])
<class 'netCDF4._netCDF4.Variable'>
vlen phony_vlen_var(y, x)
vlen data type: int32
unlimited dimensions: 
current shape = (4, 3)
>>> print(f.vltypes["phony_vlen"])
<class 'netCDF4._netCDF4.VLType'>: name = 'phony_vlen', numpy dtype = int32
```

Numpy object arrays containing python strings can also be written as vlen
variables,  For vlen strings, you don't need to create a vlen data type.
Instead, simply use the python `str` builtin (or a numpy string datatype
with fixed length greater than 1) when calling the
`Dataset.createVariable` method.

```python
>>> z = f.createDimension("z",10)
>>> strvar = f.createVariable("strvar", str, "z")
```

In this example, an object array is filled with random python strings with
random lengths between 2 and 12 characters, and the data in the object
array is assigned to the vlen string variable.

```python
>>> chars = "1234567890aabcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
>>> data = np.empty(10,"O")
>>> for n in range(10):
...     stringlen = random.randint(2,12)
...     data[n] = "".join([random.choice(chars) for i in range(stringlen)])
>>> strvar[:] = data
>>> print("variable-length string variable:\\n{}".format(strvar[:]))
variable-length string variable:
['Lh' '25F8wBbMI' '53rmM' 'vvjnb3t63ao' 'qjRBQk6w' 'aJh' 'QF'
 'jtIJbJACaQk4' '3Z5' 'bftIIq']
>>> print(f)
<class 'netCDF4._netCDF4.Dataset'>
root group (NETCDF4 data model, file format HDF5):
    dimensions(sizes): x(3), y(4), z(10)
    variables(dimensions): int32 phony_vlen_var(y,x), <class 'str'> strvar(z)
    groups: 
>>> print(f.variables["strvar"])
<class 'netCDF4._netCDF4.Variable'>
vlen strvar(z)
vlen data type: <class 'str'>
unlimited dimensions: 
current shape = (10,)
```

It is also possible to set contents of vlen string variables with numpy arrays
of any string or unicode data type. Note, however, that accessing the contents
of such variables will always return numpy arrays with dtype `object`.

## Enum data type

netCDF4 has an enumerated data type, which is an integer datatype that is
restricted to certain named values. Since Enums don't map directly to
a numpy data type, they are read and written as integer arrays.

Here's an example of using an Enum type to hold cloud type data.
The base integer data type and a python dictionary describing the allowed
values and their names are used to define an Enum data type using
`Dataset.createEnumType`.

```python
>>> nc = Dataset('clouds.nc','w')
>>> # python dict with allowed values and their names.
>>> enum_dict = {'Altocumulus': 7, 'Missing': 255,
... 'Stratus': 2, 'Clear': 0,
... 'Nimbostratus': 6, 'Cumulus': 4, 'Altostratus': 5,
... 'Cumulonimbus': 1, 'Stratocumulus': 3}
>>> # create the Enum type called 'cloud_t'.
>>> cloud_type = nc.createEnumType(np.uint8,'cloud_t',enum_dict)
>>> print(cloud_type)
<class 'netCDF4._netCDF4.EnumType'>: name = 'cloud_t', numpy dtype = uint8, fields/values ={'Altocumulus': 7, 'Missing': 255, 'Stratus': 2, 'Clear': 0, 'Nimbostratus': 6, 'Cumulus': 4, 'Altostratus': 5, 'Cumulonimbus': 1, 'Stratocumulus': 3}
```

A new variable can be created in the usual way using this data type.
Integer data is written to the variable that represents the named
cloud types in enum_dict. A `ValueError` will be raised if an attempt
is made to write an integer value not associated with one of the
specified names.

```python
>>> time = nc.createDimension('time',None)
>>> # create a 1d variable of type 'cloud_type'.
>>> # The fill_value is set to the 'Missing' named value.
>>> cloud_var = nc.createVariable('primary_cloud',cloud_type,'time',
...                               fill_value=enum_dict['Missing'])
>>> # write some data to the variable.
>>> cloud_var[:] = [enum_dict[k] for k in ['Clear', 'Stratus', 'Cumulus',
...                                        'Missing', 'Cumulonimbus']]
>>> nc.close()
>>> # reopen the file, read the data.
>>> nc = Dataset('clouds.nc')
>>> cloud_var = nc.variables['primary_cloud']
>>> print(cloud_var)
<class 'netCDF4._netCDF4.Variable'>
enum primary_cloud(time)
    _FillValue: 255
enum data type: uint8
unlimited dimensions: time
current shape = (5,)
>>> print(cloud_var.datatype.enum_dict)
{'Altocumulus': 7, 'Missing': 255, 'Stratus': 2, 'Clear': 0, 'Nimbostratus': 6, 'Cumulus': 4, 'Altostratus': 5, 'Cumulonimbus': 1, 'Stratocumulus': 3}
>>> print(cloud_var[:])
[0 2 4 -- 1]
>>> nc.close()
```

## Parallel IO

If MPI parallel enabled versions of netcdf and hdf5 or pnetcdf are detected,
and [mpi4py](https://mpi4py.scipy.org) is installed, netcdf4-python will
be built with parallel IO capabilities enabled. Parallel IO of NETCDF4 or 
NETCDF4_CLASSIC formatted files is only available if the MPI parallel HDF5
library is available. Parallel IO of classic netcdf-3 file formats is only
available if the [PnetCDF](https://parallel-netcdf.github.io/) library is
available. To use parallel IO, your program must be running in an MPI
environment using [mpi4py](https://mpi4py.scipy.org).

```python
>>> from mpi4py import MPI
>>> import numpy as np
>>> from netCDF4 import Dataset
>>> rank = MPI.COMM_WORLD.rank  # The process ID (integer 0-3 for 4-process run)
```

To run an MPI-based parallel program like this, you must use `mpiexec` to launch several
parallel instances of Python (for example, using `mpiexec -np 4 python mpi_example.py`).
The parallel features of netcdf4-python are mostly transparent -
when a new dataset is created or an existing dataset is opened,
use the `parallel` keyword to enable parallel access.

```python
>>> nc = Dataset('parallel_test.nc','w',parallel=True)
```

The optional `comm` keyword may be used to specify a particular
MPI communicator (`MPI_COMM_WORLD` is used by default).  Each process (or rank)
can now write to the file indepedently.  In this example the process rank is
written to a different variable index on each task

```python
>>> d = nc.createDimension('dim',4)
>>> v = nc.createVariable('var', np.int64, 'dim')
>>> v[rank] = rank
>>> nc.close()

% ncdump parallel_test.nc
netcdf parallel_test {
dimensions:
    dim = 4 ;
variables:
    int64 var(dim) ;
data:

    var = 0, 1, 2, 3 ;
}
```

There are two types of parallel IO, independent (the default) and collective.
Independent IO means that each process can do IO independently. It should not
depend on or be affected by other processes. Collective IO is a way of doing
IO defined in the MPI-IO standard; unlike independent IO, all processes must
participate in doing IO. To toggle back and forth between
the two types of IO, use the `Variable.set_collective`
`Variable` method. All metadata
operations (such as creation of groups, types, variables, dimensions, or attributes)
are collective.  There are a couple of important limitations of parallel IO:

 - parallel IO for NETCDF4 or NETCDF4_CLASSIC formatted files is only available
   if the netcdf library was compiled with MPI enabled HDF5.
 - parallel IO for all classic netcdf-3 file formats is only available if the
   netcdf library was compiled with [PnetCDF](https://parallel-netcdf.github.io).
 - If a variable has an unlimited dimension, appending data must be done in collective mode.
   If the write is done in independent mode, the operation will fail with a
   a generic "HDF Error".
 - You can write compressed data in parallel only with netcdf-c >= 4.7.4
   and hdf5 >= 1.10.3 (although you can read in parallel with earlier versions). To write
   compressed data in parallel, the variable must be in 'collective IO mode'.  This is done
   automatically on variable creation if compression is turned on, but if you are appending
   to a variable in an existing file, you must use `Variable.set_collective(True)` before attempting
   to write to it.
 - You cannot use variable-length (VLEN) data types.

## Dealing with strings

The most flexible way to store arrays of strings is with the
[Variable-length (vlen) string data type](#variable-length-vlen-data-type). However, this requires
the use of the NETCDF4 data model, and the vlen type does not map very well
numpy arrays (you have to use numpy arrays of dtype=`object`, which are arrays of
arbitrary python objects). numpy does have a fixed-width string array
data type, but unfortunately the netCDF data model does not.
Instead fixed-width byte strings are typically stored as [arrays of 8-bit
characters](https://www.unidata.ucar.edu/software/netcdf/docs/BestPractices.html#bp_Strings-and-Variables-of-type-char).
To perform the conversion to and from character arrays to fixed-width numpy string arrays, the
following convention is followed by the python interface.
If the `_Encoding` special attribute is set for a character array
(dtype `S1`) variable, the `chartostring` utility function is used to convert the array of
characters to an array of strings with one less dimension (the last dimension is
interpreted as the length of each string) when reading the data. The character
set (usually ascii) is specified by the `_Encoding` attribute. If `_Encoding`
is 'none' or 'bytes', then the character array is converted to a numpy
fixed-width byte string array (dtype `S#`), otherwise a numpy unicode (dtype
`U#`) array is created.  When writing the data,
`stringtochar` is used to convert the numpy string array to an array of
characters with one more dimension. For example,

```python
>>> from netCDF4 import stringtochar
>>> nc = Dataset('stringtest.nc','w',format='NETCDF4_CLASSIC')
>>> _ = nc.createDimension('nchars',3)
>>> _ = nc.createDimension('nstrings',None)
>>> v = nc.createVariable('strings','S1',('nstrings','nchars'))
>>> datain = np.array(['foo','bar'],dtype='S3')
>>> v[:] = stringtochar(datain) # manual conversion to char array
>>> print(v[:]) # data returned as char array
[[b'f' b'o' b'o']
 [b'b' b'a' b'r']]
>>> v._Encoding = 'ascii' # this enables automatic conversion
>>> v[:] = datain # conversion to char array done internally
>>> print(v[:])  # data returned in numpy string array
['foo' 'bar']
>>> nc.close()
```

Even if the `_Encoding` attribute is set, the automatic conversion of char
arrays to/from string arrays can be disabled with
`Variable.set_auto_chartostring`.

A similar situation is often encountered with numpy structured arrays with subdtypes
containing fixed-wdith byte strings (dtype=`S#`). Since there is no native fixed-length string
netCDF datatype, these numpy structure arrays are mapped onto netCDF compound
types with character array elements.  In this case the string <-> char array
conversion is handled automatically (without the need to set the `_Encoding`
attribute) using [numpy
views](https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.view.html).
The structured array dtype (including the string elements) can even be used to
define the compound data type - the string dtype will be converted to
character array dtype under the hood when creating the netcdf compound type.
Here's an example:

```python
>>> nc = Dataset('compoundstring_example.nc','w')
>>> dtype = np.dtype([('observation', 'f4'),
...                      ('station_name','S10')])
>>> station_data_t = nc.createCompoundType(dtype,'station_data')
>>> _ = nc.createDimension('station',None)
>>> statdat = nc.createVariable('station_obs', station_data_t, ('station',))
>>> data = np.empty(2,dtype)
>>> data['observation'][:] = (123.,3.14)
>>> data['station_name'][:] = ('Boulder','New York')
>>> print(statdat.dtype) # strings actually stored as character arrays
{'names':['observation','station_name'], 'formats':['<f4',('S1', (10,))], 'offsets':[0,4], 'itemsize':16, 'aligned':True}
>>> statdat[:] = data # strings converted to character arrays internally
>>> print(statdat[:])  # character arrays converted back to strings
[(123.  , b'Boulder') (  3.14, b'New York')]
>>> print(statdat[:].dtype)
{'names':['observation','station_name'], 'formats':['<f4','S10'], 'offsets':[0,4], 'itemsize':16, 'aligned':True}
>>> statdat.set_auto_chartostring(False) # turn off auto-conversion
>>> statdat[:] = data.view(dtype=[('observation', 'f4'),('station_name','S1',10)])
>>> print(statdat[:])  # now structured array with char array subtype is returned
[(123.  , [b'B', b'o', b'u', b'l', b'd', b'e', b'r', b'', b'', b''])
 (  3.14, [b'N', b'e', b'w', b' ', b'Y', b'o', b'r', b'k', b'', b''])]
>>> nc.close()
```

Note that there is currently no support for mapping numpy structured arrays with
unicode elements (dtype `U#`) onto netCDF compound types, nor is there support
for netCDF compound types with vlen string components.

## In-memory (diskless) Datasets

You can create netCDF Datasets whose content is held in memory
instead of in a disk file.  There are two ways to do this.  If you
don't need access to the memory buffer containing the Dataset from
within python, the best way is to use the `diskless=True` keyword
argument when creating the Dataset.  If you want to save the Dataset
to disk when you close it, also set `persist=True`.  If you want to
create a new read-only Dataset from an existing python memory buffer, use the
`memory` keyword argument to pass the memory buffer when creating the Dataset.
If you want to create a new in-memory Dataset, and then access the memory buffer
directly from Python, use the `memory` keyword argument to specify the
estimated size of the Dataset in bytes when creating the Dataset with
`mode='w'`.  Then, the `Dataset.close` method will return a python memoryview
object representing the Dataset. Below are examples illustrating both
approaches.

```python
>>> # create a diskless (in-memory) Dataset,
>>> # and persist the file to disk when it is closed.
>>> nc = Dataset('diskless_example.nc','w',diskless=True,persist=True)
>>> d = nc.createDimension('x',None)
>>> v = nc.createVariable('v',np.int32,'x')
>>> v[0:5] = np.arange(5)
>>> print(nc)
<class 'netCDF4._netCDF4.Dataset'>
root group (NETCDF4 data model, file format HDF5):
    dimensions(sizes): x(5)
    variables(dimensions): int32 v(x)
    groups: 
>>> print(nc['v'][:])
[0 1 2 3 4]
>>> nc.close() # file saved to disk
>>> # create an in-memory dataset from an existing python
>>> # python memory buffer.
>>> # read the newly created netcdf file into a python
>>> # bytes object.
>>> with open('diskless_example.nc', 'rb') as f:
...     nc_bytes = f.read()
>>> # create a netCDF in-memory dataset from the bytes object.
>>> nc = Dataset('inmemory.nc', memory=nc_bytes)
>>> print(nc)
<class 'netCDF4._netCDF4.Dataset'>
root group (NETCDF4 data model, file format HDF5):
    dimensions(sizes): x(5)
    variables(dimensions): int32 v(x)
    groups: 
>>> print(nc['v'][:])
[0 1 2 3 4]
>>> nc.close()
>>> # create an in-memory Dataset and retrieve memory buffer
>>> # estimated size is 1028 bytes - this is actually only
>>> # used if format is NETCDF3
>>> # (ignored for NETCDF4/HDF5 files).
>>> nc = Dataset('inmemory.nc', mode='w',memory=1028)
>>> d = nc.createDimension('x',None)
>>> v = nc.createVariable('v',np.int32,'x')
>>> v[0:5] = np.arange(5)
>>> nc_buf = nc.close() # close returns memoryview
>>> print(type(nc_buf))
<class 'memoryview'>
>>> # save nc_buf to disk, read it back in and check.
>>> with open('inmemory.nc', 'wb') as f:
...     f.write(nc_buf)
>>> nc = Dataset('inmemory.nc')
>>> print(nc)
<class 'netCDF4._netCDF4.Dataset'>
root group (NETCDF4 data model, file format HDF5):
    dimensions(sizes): x(5)
    variables(dimensions): int32 v(x)
    groups:
>>> print(nc['v'][:])
[0 1 2 3 4]
>>> nc.close()
```


All of the code in this tutorial is available in `examples/tutorial.py`, except
the parallel IO example, which is in `examples/mpi_example.py`.
Unit tests are in the `test` directory.

**contact**: Jeffrey Whitaker <jeffrey.s.whitaker@noaa.gov>

**copyright**: 2008 by Jeffrey Whitaker.

**license**: Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

# Make changes to this file, not the c-wrappers that Cython generates.
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from cpython.buffer cimport PyObject_GetBuffer, PyBuffer_Release, PyBUF_SIMPLE, PyBUF_ANY_CONTIGUOUS
from cpython.bytes cimport PyBytes_FromStringAndSize

# pure python utilities
from .utils import (_StartCountStride, _quantize, _find_dim, _walk_grps,
                    _out_array_shape, _sortbylist, _tostr, _safecast, _is_int)
import sys

__version__ = "1.6.3"

# Initialize numpy
import posixpath
from cftime import date2num, num2date, date2index
import numpy
cimport numpy
import weakref
import warnings
import subprocess
import pathlib
import os
from glob import glob
from numpy import ma
from libc.string cimport memcpy, memset
from libc.stdlib cimport malloc, free
numpy.import_array()
include "constants.pyx"
include "membuf.pyx"
include "netCDF4.pxi"
IF HAS_PARALLEL4_SUPPORT or HAS_PNETCDF_SUPPORT:
    cimport mpi4py.MPI as MPI
    from mpi4py.libmpi cimport MPI_Comm, MPI_Info, MPI_Comm_dup, MPI_Info_dup, \
                               MPI_Comm_free, MPI_Info_free, MPI_INFO_NULL,\
                               MPI_COMM_WORLD
    ctypedef MPI.Comm Comm
    ctypedef MPI.Info Info
ELSE:
    ctypedef object Comm
    ctypedef object Info

# check for required version of netcdf-4 and hdf5.

def _gethdf5libversion():
    cdef unsigned int majorvers, minorvers, releasevers
    cdef herr_t ierr
    with nogil:
        ierr = H5get_libversion( &majorvers, &minorvers, &releasevers)
    if ierr < 0:
        raise RuntimeError('error getting HDF5 library version info')
    return '%d.%d.%d' % (majorvers,minorvers,releasevers)

def getlibversion():
    """
**`getlibversion()`**

returns a string describing the version of the netcdf library
used to build the module, and when it was built.
    """
    return (<char *>nc_inq_libvers()).decode('ascii')

def get_chunk_cache():
    """
**`get_chunk_cache()`**

return current netCDF chunk cache information in a tuple (size,nelems,preemption).
See netcdf C library documentation for `nc_get_chunk_cache` for
details. Values can be reset with `set_chunk_cache`."""
    cdef int ierr
    cdef size_t sizep, nelemsp
    cdef float preemptionp
    with nogil:
        ierr = nc_get_chunk_cache(&sizep, &nelemsp, &preemptionp)
    _ensure_nc_success(ierr)
    size = sizep; nelems = nelemsp; preemption = preemptionp
    return (size,nelems,preemption)

def set_chunk_cache(size=None,nelems=None,preemption=None):
    """
**`set_chunk_cache(self,size=None,nelems=None,preemption=None)`**

change netCDF4 chunk cache settings.
See netcdf C library documentation for `nc_set_chunk_cache` for
details."""
    cdef int ierr
    cdef size_t sizep, nelemsp
    cdef float preemptionp
    # reset chunk cache size, leave other parameters unchanged.
    size_orig, nelems_orig, preemption_orig = get_chunk_cache()
    if size is not None:
        sizep = size
    else:
        sizep = size_orig
    if nelems is not None:
        nelemsp = nelems
    else:
        nelemsp = nelems_orig
    if preemption is not None:
        preemptionp = preemption
    else:
        preemptionp = preemption_orig
    with nogil:
        ierr = nc_set_chunk_cache(sizep,nelemsp, preemptionp)
    _ensure_nc_success(ierr)

IF HAS_SET_ALIGNMENT:
    def get_alignment():
        """
    **`get_alignment()`**

    return current netCDF alignment within HDF5 files in a tuple
    (threshold,alignment). See netcdf C library documentation for
    `nc_get_alignment` for details. Values can be reset with
    `set_alignment`.

    This function was added in netcdf 4.9.0."""
        cdef int ierr
        cdef int thresholdp, alignmentp
        ierr = nc_get_alignment(&thresholdp, &alignmentp)
        _ensure_nc_success(ierr)
        threshold = thresholdp
        alignment = alignmentp
        return (threshold,alignment)

    def set_alignment(threshold, alignment):
        """
    **`set_alignment(threshold,alignment)`**

    Change the HDF5 file alignment.
    See netcdf C library documentation for `nc_set_alignment` for
    details.

    This function was added in netcdf 4.9.0."""
        cdef int ierr
        cdef int thresholdp, alignmentp
        thresholdp = threshold
        alignmentp = alignment

        ierr = nc_set_alignment(thresholdp, alignmentp)
        _ensure_nc_success(ierr)
ELSE:
    def get_alignment():
        raise RuntimeError(
            "This function requires netcdf4 4.9.0+ to be used at compile time"
        )

    def set_alignment(threshold, alignment):
        raise RuntimeError(
            "This function requires netcdf4 4.9.0+ to be used at compile time"
        )

__netcdf4libversion__ = getlibversion().split()[0]
__hdf5libversion__ = _gethdf5libversion()
__has_rename_grp__ = HAS_RENAME_GRP
__has_nc_inq_path__ = HAS_NC_INQ_PATH
__has_nc_inq_format_extended__ = HAS_NC_INQ_FORMAT_EXTENDED
__has_cdf5_format__ = HAS_CDF5_FORMAT
__has_nc_open_mem__ = HAS_NC_OPEN_MEM
__has_nc_create_mem__ = HAS_NC_CREATE_MEM
__has_parallel4_support__ = HAS_PARALLEL4_SUPPORT
__has_pnetcdf_support__ = HAS_PNETCDF_SUPPORT
__has_quantization_support__ = HAS_QUANTIZATION_SUPPORT
__has_zstandard_support__ = HAS_ZSTANDARD_SUPPORT
__has_bzip2_support__ = HAS_BZIP2_SUPPORT
__has_blosc_support__ = HAS_BLOSC_SUPPORT
__has_szip_support__ = HAS_SZIP_SUPPORT
__has_set_alignment__ = HAS_SET_ALIGNMENT
_needsworkaround_issue485 = __netcdf4libversion__ < "4.4.0" or \
               (__netcdf4libversion__.startswith("4.4.0") and \
                "-development" in __netcdf4libversion__)

# issue warning for hdf5 1.10 (issue #549)
if __netcdf4libversion__[0:5] < "4.4.1" and\
   __hdf5libversion__.startswith("1.10"):
    msg = """
WARNING: Backwards incompatible files will be created with HDF5 1.10.x 
and netCDF < 4.4.1. Upgrading to netCDF4 >= 4.4.1 or downgrading to 
to HDF5 version 1.8.x is highly recommended 
(see https://github.com/Unidata/netcdf-c/issues/250)."""
    warnings.warn(msg)

# numpy data type <--> netCDF 4 data type mapping.
_nptonctype  = {'S1' : NC_CHAR,
                'i1' : NC_BYTE,
                'u1' : NC_UBYTE,
                'i2' : NC_SHORT,
                'u2' : NC_USHORT,
                'i4' : NC_INT,
                'u4' : NC_UINT,
                'i8' : NC_INT64,
                'u8' : NC_UINT64,
                'f4' : NC_FLOAT,
                'f8' : NC_DOUBLE}

# just integer types.
_intnptonctype  = {'i1' : NC_BYTE,
                   'u1' : NC_UBYTE,
                   'i2' : NC_SHORT,
                   'u2' : NC_USHORT,
                   'i4' : NC_INT,
                   'u4' : NC_UINT,
                   'i8' : NC_INT64,
                   'u8' : NC_UINT64}

# create dictionary mapping string identifiers to netcdf format codes
_format_dict  = {'NETCDF3_CLASSIC' : NC_FORMAT_CLASSIC,
                 'NETCDF4_CLASSIC' : NC_FORMAT_NETCDF4_CLASSIC,
                 'NETCDF4'         : NC_FORMAT_NETCDF4}
# create dictionary mapping string identifiers to netcdf create format codes
_cmode_dict  = {'NETCDF3_CLASSIC' : NC_CLASSIC_MODEL,
                'NETCDF4_CLASSIC' : NC_CLASSIC_MODEL | NC_NETCDF4,
                'NETCDF4'         : NC_NETCDF4}
# dicts for blosc, szip compressors.
_blosc_dict={'blosc_lz':0,'blosc_lz4':1,'blosc_lz4hc':2,'blosc_snappy':3,'blosc_zlib':4,'blosc_zstd':5}
_blosc_dict_inv = {v: k for k, v in _blosc_dict.items()}
_szip_dict = {'ec': 4, 'nn': 32}
_szip_dict_inv = {v: k for k, v in _szip_dict.items()}
IF HAS_CDF5_FORMAT:
    # NETCDF3_64BIT deprecated, saved for compatibility.
    # use NETCDF3_64BIT_OFFSET instead.
    _format_dict['NETCDF3_64BIT_OFFSET'] = NC_FORMAT_64BIT_OFFSET
    _format_dict['NETCDF3_64BIT_DATA'] = NC_FORMAT_64BIT_DATA
    _cmode_dict['NETCDF3_64BIT_OFFSET'] = NC_64BIT_OFFSET
    _cmode_dict['NETCDF3_64BIT_DATA'] = NC_64BIT_DATA
ELSE:
    _format_dict['NETCDF3_64BIT'] = NC_FORMAT_64BIT
    _cmode_dict['NETCDF3_64BIT'] = NC_64BIT_OFFSET
# invert dictionary mapping
_reverse_format_dict = dict((v, k) for k, v in _format_dict.iteritems())
# add duplicate entry (NETCDF3_64BIT == NETCDF3_64BIT_OFFSET)
IF HAS_CDF5_FORMAT:
    _format_dict['NETCDF3_64BIT'] = NC_FORMAT_64BIT_OFFSET
    _cmode_dict['NETCDF3_64BIT'] = NC_64BIT_OFFSET
ELSE:
    _format_dict['NETCDF3_64BIT_OFFSET'] = NC_FORMAT_64BIT
    _cmode_dict['NETCDF3_64BIT_OFFSET'] = NC_64BIT_OFFSET

# default fill_value to numpy datatype mapping.
default_fillvals = {#'S1':NC_FILL_CHAR,
                     'S1':'\0',
                     'i1':NC_FILL_BYTE,
                     'u1':NC_FILL_UBYTE,
                     'i2':NC_FILL_SHORT,
                     'u2':NC_FILL_USHORT,
                     'i4':NC_FILL_INT,
                     'u4':NC_FILL_UINT,
                     'i8':NC_FILL_INT64,
                     'u8':NC_FILL_UINT64,
                     'f4':NC_FILL_FLOAT,
                     'f8':NC_FILL_DOUBLE}

# logical for native endian type.
is_native_little = numpy.dtype('<f4').byteorder == c'='
is_native_big = numpy.dtype('>f4').byteorder == c'='

# hard code these here, instead of importing from netcdf.h
# so it will compile with versions <= 4.2.
NC_DISKLESS = 0x0008
# introduced in 4.6.2
if __netcdf4libversion__[0:5] >= "4.6.2":
    NC_PERSIST = 0x4000
else:  # prior to 4.6.2 this flag doesn't work, so make the same as NC_DISKLESS
    NC_PERSIST = NC_DISKLESS

# next two lines do nothing, preserved for backwards compatibility.
default_encoding = 'utf-8'
unicode_error = 'replace'

_nctonptype = {}
for _key,_value in _nptonctype.items():
    _nctonptype[_value] = _key
_supportedtypes = _nptonctype.keys()
# make sure NC_CHAR points to S1
_nctonptype[NC_CHAR]='S1'

# internal C functions.

cdef _get_att_names(int grpid, int varid):
    # Private function to get all the attribute names in a group
    cdef int ierr, numatts, n
    cdef char namstring[NC_MAX_NAME+1]
    if varid == NC_GLOBAL:
        with nogil:
            ierr = nc_inq_natts(grpid, &numatts)
    else:
        with nogil:
            ierr = nc_inq_varnatts(grpid, varid, &numatts)
    _ensure_nc_success(ierr, err_cls=AttributeError)
    attslist = []
    for n from 0 <= n < numatts:
        with nogil:
            ierr = nc_inq_attname(grpid, varid, n, namstring)
        _ensure_nc_success(ierr, err_cls=AttributeError)
        # attribute names are assumed to be utf-8
        attslist.append(namstring.decode('utf-8'))
    return attslist

cdef _get_att(grp, int varid, name, encoding='utf-8'):
    # Private function to get an attribute value given its name
    cdef int ierr, n, _grpid
    cdef size_t att_len
    cdef char *attname
    cdef nc_type att_type
    cdef ndarray value_arr
    # attribute names are assumed to be utf-8
    bytestr = _strencode(name,encoding='utf-8')
    attname = bytestr
    _grpid = grp._grpid
    with nogil:
        ierr = nc_inq_att(_grpid, varid, attname, &att_type, &att_len)
    _ensure_nc_success(ierr, err_cls=AttributeError)
    # attribute is a character or string ...
    if att_type == NC_CHAR:
        value_arr = numpy.empty(att_len,'S1')
        with nogil:
            ierr = nc_get_att_text(_grpid, varid, attname,
                    PyArray_BYTES(value_arr))
        _ensure_nc_success(ierr, err_cls=AttributeError)
        if name == '_FillValue':
            # make sure _FillValue for character arrays is a byte on python 3
            # (issue 271).
            pstring = value_arr.tobytes()
        else:
            pstring =\
            value_arr.tobytes().decode(encoding,errors='replace').replace('\x00','')
        return pstring
    elif att_type == NC_STRING:
        values = <char**>PyMem_Malloc(sizeof(char*) * att_len)
        if not values:
            raise MemoryError()
        try:
            with nogil:
                ierr = nc_get_att_string(_grpid, varid, attname, values)
            _ensure_nc_success(ierr, err_cls=AttributeError)
            try:
                result = [values[j].decode(encoding,errors='replace').replace('\x00','')
                          if values[j] else "" for j in range(att_len)]
            finally:
                with nogil:
                    ierr = nc_free_string(att_len, values) # free memory in netcdf C lib
        finally:
            PyMem_Free(values)

        if len(result) == 1:
            return result[0]
        else:
            return result
    else:
    # a regular numeric or compound type.
        if att_type == NC_LONG:
            att_type = NC_INT
        try:
            type_att = _nctonptype[att_type] # see if it is a primitive type
            value_arr = numpy.empty(att_len,type_att)
        except KeyError:
            # check if it's a compound
            try:
                type_att = _read_compound(grp, att_type)
                value_arr = numpy.empty(att_len,type_att.dtype_view)
            except:
                # check if it's an enum
                try:
                    type_att = _read_enum(grp, att_type)
                    value_arr = numpy.empty(att_len,type_att.dtype)
                except:
                    raise KeyError('attribute %s has unsupported datatype' % attname)
        with nogil:
            ierr = nc_get_att(_grpid, varid, attname, PyArray_BYTES(value_arr))
        _ensure_nc_success(ierr, err_cls=AttributeError)
        if value_arr.shape == ():
            # return a scalar for a scalar array
            return value_arr.item()
        elif att_len == 1:
            # return a scalar for a single element array
            return value_arr[0]
        else:
            return value_arr

def _set_default_format(object format='NETCDF4'):
    # Private function to set the netCDF file format
    cdef int ierr, formatid
    if format not in _format_dict:
        raise ValueError("unrecognized format requested")
    formatid = _format_dict[format]
    with nogil:
        ierr = nc_set_default_format(formatid, NULL)
    _ensure_nc_success(ierr)

cdef _get_format(int grpid):
    # Private function to get the netCDF file format
    cdef int ierr, formatp
    with nogil:
        ierr = nc_inq_format(grpid, &formatp)
    _ensure_nc_success(ierr)
    if formatp not in _reverse_format_dict:
        raise ValueError('format not supported by python interface')
    return _reverse_format_dict[formatp]

cdef _get_full_format(int grpid):
    # Private function to get the underlying disk format
    cdef int ierr, formatp, modep
    IF HAS_NC_INQ_FORMAT_EXTENDED:
        with nogil:
            ierr = nc_inq_format_extended(grpid, &formatp, &modep)
        _ensure_nc_success(ierr)
        if formatp == NC_FORMAT_NC3:
            return 'NETCDF3'
        elif formatp == NC_FORMAT_NC_HDF5:
            return 'HDF5'
        elif formatp == NC_FORMAT_NC_HDF4:
            return 'HDF4'
        elif formatp == NC_FORMAT_PNETCDF:
            return 'PNETCDF'
        elif formatp == NC_FORMAT_DAP2:
            return 'DAP2'
        elif formatp == NC_FORMAT_DAP4:
            return 'DAP4'
        elif formatp == NC_FORMAT_UNDEFINED:
            return 'UNDEFINED'
    ELSE:
        return 'UNDEFINED'

cdef issue485_workaround(int grpid, int varid, char* attname):
    # check to see if attribute already exists
    # and is NC_CHAR, if so delete it and re-create it
    # (workaround for issue #485). Fixed in C library
    # with commit 473259b7728120bb281c52359b1af50cca2fcb72,
    # which was included in 4.4.0-RC5.
    cdef nc_type att_type
    cdef size_t att_len

    if not _needsworkaround_issue485:
        return
    with nogil:
        ierr = nc_inq_att(grpid, varid, attname, &att_type, &att_len)
    if ierr == NC_NOERR and att_type == NC_CHAR:
        with nogil:
            ierr = nc_del_att(grpid, varid, attname)
        _ensure_nc_success(ierr)


cdef _set_att(grp, int varid, name, value,\
              nc_type xtype=-99, force_ncstring=False):
    # Private function to set an attribute name/value pair
    cdef int ierr, lenarr, N, grpid
    cdef char *attname
    cdef char *datstring
    cdef char **string_ptrs
    cdef ndarray value_arr
    bytestr = _strencode(name)
    attname = bytestr
    grpid = grp._grpid
    # put attribute value into a numpy array.
    value_arr = numpy.array(value)
    if value_arr.ndim > 1: # issue #841
        if __version__ > "1.4.2":
            raise ValueError('multi-dimensional array attributes not supported')
        else:
            msg = """
Multi-dimensional array attributes are now deprecated.
Instead of silently flattening the array, an error will
be raised in the next release."""
            warnings.warn(msg,FutureWarning)
    # if array is 64 bit integers or
    # if 64-bit datatype not supported, cast to 32 bit integers.
    fmt = _get_format(grpid)
    is_netcdf3 = fmt.startswith('NETCDF3') or fmt == 'NETCDF4_CLASSIC'
    if value_arr.dtype.str[1:] == 'i8' and ('i8' not in _supportedtypes or\
       (is_netcdf3 and fmt != 'NETCDF3_64BIT_DATA')):
        value_arr = value_arr.astype('i4')
    # if array contains ascii strings, write a text attribute (stored as bytes).
    # if array contains unicode strings, and data model is NETCDF4,
    # write as a string.
    if value_arr.dtype.char in ['S','U']:
        # force array of strings if array has multiple elements (issue #770)
        N = value_arr.size
        if N > 1: force_ncstring=True
        if not is_netcdf3 and force_ncstring and N > 1:
            string_ptrs = <char**>PyMem_Malloc(N * sizeof(char*))
            if not string_ptrs:
                raise MemoryError()
            try:
                strings = [_strencode(s) for s in value_arr.flat]
                for j in range(N):
                    if len(strings[j]) == 0:
                        strings[j] = _strencode('\x00')
                    string_ptrs[j] = strings[j]
                issue485_workaround(grpid, varid, attname)
                with nogil:
                    ierr = nc_put_att_string(grpid, varid, attname, N, string_ptrs)
            finally:
                PyMem_Free(string_ptrs)
        else:
            # don't allow string array attributes in NETCDF3 files.
            if is_netcdf3 and N > 1:
                msg='array string attributes can only be written with NETCDF4'
                raise OSError(msg)
            if not value_arr.shape:
                dats = _strencode(value_arr.item())
            else:
                value_arr1 = value_arr.ravel()
                dats = _strencode(''.join(value_arr1.tolist()))
            lenarr = len(dats)
            datstring = dats
            if lenarr == 0:
                # write null byte
                lenarr=1; datstring = '\x00'
            if (force_ncstring or value_arr.dtype.char == 'U') and not is_netcdf3:
                # try to convert to ascii string, write as NC_CHAR
                # else it's a unicode string, write as NC_STRING (if NETCDF4)
                try:
                    if force_ncstring: raise UnicodeError
                    dats_ascii = _to_ascii(dats) # try to encode bytes as ascii string
                    with nogil:
                        ierr = nc_put_att_text(grpid, varid, attname, lenarr, datstring)
                except UnicodeError:
                    issue485_workaround(grpid, varid, attname)
                    with nogil:
                        ierr = nc_put_att_string(grpid, varid, attname, 1, &datstring)
            else:
                with nogil:
                    ierr = nc_put_att_text(grpid, varid, attname, lenarr, datstring)
        _ensure_nc_success(ierr, err_cls=AttributeError)
    # a 'regular' array type ('f4','i4','f8' etc)
    else:
        if value_arr.dtype.kind == 'V': # compound attribute.
            xtype = _find_cmptype(grp,value_arr.dtype)
        elif value_arr.dtype.str[1:] not in _supportedtypes:
            raise TypeError, 'illegal data type for attribute %r, must be one of %s, got %s' % (attname, _supportedtypes, value_arr.dtype.str[1:])
        elif xtype == -99: # if xtype is not passed in as kwarg.
            xtype = _nptonctype[value_arr.dtype.str[1:]]
        lenarr = PyArray_SIZE(value_arr)
        with nogil:
            ierr = nc_put_att(grpid, varid, attname, xtype, lenarr,
                PyArray_DATA(value_arr))
        _ensure_nc_success(ierr, err_cls=AttributeError)

cdef _get_types(group):
    # Private function to create `CompoundType`,
    # `VLType` or `EnumType` instances for all the
    # compound, VLEN or Enum types in a `Group` or `Dataset`.
    cdef int ierr, ntypes, classp, n, _grpid
    cdef nc_type xtype
    cdef nc_type *typeids
    cdef char namstring[NC_MAX_NAME+1]
    _grpid = group._grpid
    # get the number of user defined types in this group.
    with nogil:
        ierr = nc_inq_typeids(_grpid, &ntypes, NULL)
    _ensure_nc_success(ierr)
    if ntypes > 0:
        typeids = <nc_type *>malloc(sizeof(nc_type) * ntypes)
        with nogil:
            ierr = nc_inq_typeids(_grpid, &ntypes, typeids)
        _ensure_nc_success(ierr)
    # create empty dictionary for CompoundType instances.
    cmptypes = dict()
    vltypes = dict()
    enumtypes = dict()

    if ntypes > 0:
        for n from 0 <= n < ntypes:
            xtype = typeids[n]
            with nogil:
                ierr = nc_inq_user_type(_grpid, xtype, namstring,
                                        NULL,NULL,NULL,&classp)
            _ensure_nc_success(ierr)
            if classp == NC_COMPOUND: # a compound
                name = namstring.decode('utf-8')
                # read the compound type info from the file,
                # create a CompoundType instance from it.
                try:
                    cmptype = _read_compound(group, xtype)
                except KeyError:
                    msg='WARNING: unsupported Compound type, skipping...'
                    warnings.warn(msg)
                    continue
                cmptypes[name] = cmptype
            elif classp == NC_VLEN: # a vlen
                name = namstring.decode('utf-8')
                # read the VLEN type info from the file,
                # create a VLType instance from it.
                try:
                    vltype = _read_vlen(group, xtype)
                except KeyError:
                    msg='WARNING: unsupported VLEN type, skipping...'
                    warnings.warn(msg)
                    continue
                vltypes[name] = vltype
            elif classp == NC_ENUM: # an enum type
                name = namstring.decode('utf-8')
                # read the Enum type info from the file,
                # create a EnumType instance from it.
                try:
                    enumtype = _read_enum(group, xtype)
                except KeyError:
                    msg='WARNING: unsupported Enum type, skipping...'
                    warnings.warn(msg)
                    continue
                enumtypes[name] = enumtype
        free(typeids)
    return cmptypes, vltypes, enumtypes

cdef _get_dims(group):
    # Private function to create `Dimension` instances for all the
    # dimensions in a `Group` or Dataset
    cdef int ierr, numdims, n, _grpid
    cdef int *dimids
    cdef char namstring[NC_MAX_NAME+1]
    # get number of dimensions in this Group.
    _grpid = group._grpid
    with nogil:
        ierr = nc_inq_ndims(_grpid, &numdims)
    _ensure_nc_success(ierr)
    # create empty dictionary for dimensions.
    dimensions = dict()
    if numdims > 0:
        dimids = <int *>malloc(sizeof(int) * numdims)
        if group.data_model == 'NETCDF4':
            with nogil:
                ierr = nc_inq_dimids(_grpid, &numdims, dimids, 0)
            _ensure_nc_success(ierr)
        else:
            for n from 0 <= n < numdims:
                dimids[n] = n
        for n from 0 <= n < numdims:
            with nogil:
                ierr = nc_inq_dimname(_grpid, dimids[n], namstring)
            _ensure_nc_success(ierr)
            name = namstring.decode('utf-8')
            dimensions[name] = Dimension(group, name, id=dimids[n])
        free(dimids)
    return dimensions

cdef _get_grps(group):
    # Private function to create `Group` instances for all the
    # groups in a `Group` or Dataset
    cdef int ierr, numgrps, n, _grpid
    cdef int *grpids
    cdef char namstring[NC_MAX_NAME+1]
    # get number of groups in this Group.
    _grpid = group._grpid
    with nogil:
        ierr = nc_inq_grps(_grpid, &numgrps, NULL)
    _ensure_nc_success(ierr)
    # create dictionary containing `Group` instances for groups in this group
    groups = dict()
    if numgrps > 0:
        grpids = <int *>malloc(sizeof(int) * numgrps)
        with nogil:
            ierr = nc_inq_grps(_grpid, NULL, grpids)
        _ensure_nc_success(ierr)
        for n from 0 <= n < numgrps:
            with nogil:
                ierr = nc_inq_grpname(grpids[n], namstring)
            _ensure_nc_success(ierr)
            name = namstring.decode('utf-8')
            groups[name] = Group(group, name, id=grpids[n])
        free(grpids)
    return groups

cdef _get_vars(group):
    # Private function to create `Variable` instances for all the
    # variables in a `Group` or Dataset
    cdef int ierr, numvars, n, nn, numdims, varid, classp, iendian, _grpid
    cdef int *varids
    cdef int *dimids
    cdef nc_type xtype
    cdef char namstring[NC_MAX_NAME+1]
    cdef char namstring_cmp[NC_MAX_NAME+1]
    # get number of variables in this Group.
    _grpid = group._grpid
    with nogil:
        ierr = nc_inq_nvars(_grpid, &numvars)
    _ensure_nc_success(ierr, err_cls=AttributeError)
    # create empty dictionary for variables.
    variables = dict()
    if numvars > 0:
        # get variable ids.
        varids = <int *>malloc(sizeof(int) * numvars)
        if group.data_model == 'NETCDF4':
            with nogil:
                ierr = nc_inq_varids(_grpid, &numvars, varids)
            _ensure_nc_success(ierr)
        else:
            for n from 0 <= n < numvars:
                varids[n] = n
        # loop over variables.
        for n from 0 <= n < numvars:
            varid = varids[n]
            # get variable name.
            with nogil:
                ierr = nc_inq_varname(_grpid, varid, namstring)
            _ensure_nc_success(ierr)
            name = namstring.decode('utf-8')
            # get variable type.
            with nogil:
                ierr = nc_inq_vartype(_grpid, varid, &xtype)
            _ensure_nc_success(ierr)
            # get endian-ness of variable.
            endianness = None
            with nogil:
                ierr = nc_inq_var_endian(_grpid, varid, &iendian)
            if ierr == NC_NOERR and iendian == NC_ENDIAN_LITTLE:
                endianness = '<'
            elif iendian == NC_ENDIAN_BIG:
                endianness = '>'
            # check to see if it is a supported user-defined type.
            try:
                datatype = _nctonptype[xtype]
                if endianness is not None:
                    datatype = endianness + datatype
            except KeyError:
                if xtype == NC_STRING:
                    datatype = str
                else:
                    with nogil:
                        ierr = nc_inq_user_type(_grpid, xtype, namstring_cmp,
                                                NULL, NULL, NULL, &classp)
                    _ensure_nc_success(ierr)
                    if classp == NC_COMPOUND: # a compound type
                        # create CompoundType instance describing this compound type.
                        try:
                            datatype = _read_compound(group, xtype, endian=endianness)
                        except KeyError:
                            msg="WARNING: variable '%s' has unsupported compound datatype, skipping .." % name
                            warnings.warn(msg)
                            continue
                    elif classp == NC_VLEN: # a compound type
                        # create VLType instance describing this compound type.
                        try:
                            datatype = _read_vlen(group, xtype, endian=endianness)
                        except KeyError:
                            msg="WARNING: variable '%s' has unsupported VLEN datatype, skipping .." % name
                            warnings.warn(msg)
                            continue
                    elif classp == NC_ENUM:
                        # create EnumType instance describing this compound type.
                        try:
                            datatype = _read_enum(group, xtype, endian=endianness)
                        except KeyError:
                            msg="WARNING: variable '%s' has unsupported Enum datatype, skipping .." % name
                            warnings.warn(msg)
                            continue
                    else:
                        msg="WARNING: variable '%s' has unsupported datatype, skipping .." % name
                        warnings.warn(msg)
                        continue
            # get number of dimensions.
            with nogil:
                ierr = nc_inq_varndims(_grpid, varid, &numdims)
            _ensure_nc_success(ierr)
            dimids = <int *>malloc(sizeof(int) * numdims)
            # get dimension ids.
            with nogil:
                ierr = nc_inq_vardimid(_grpid, varid, dimids)
            _ensure_nc_success(ierr)
            # loop over dimensions, retrieve names.
            # if not found in current group, look in parents.
            # QUESTION:  what if grp1 has a dimension named 'foo'
            # and so does it's parent - can a variable in grp1
            # use the 'foo' dimension from the parent?
            dimensions = []
            for nn from 0 <= nn < numdims:
                grp = group
                found = False
                while not found:
                    for key, value in grp.dimensions.items():
                        if value._dimid == dimids[nn]:
                            dimensions.append(key)
                            found = True
                            break
                    grp = grp.parent
            free(dimids)
            # create new variable instance.
            dimensions = tuple(_find_dim(group,d) for d in dimensions)
            if endianness == '>':
                variables[name] = Variable(group, name, datatype, dimensions, id=varid, endian='big')
            elif endianness == '<':
                variables[name] = Variable(group, name, datatype, dimensions, id=varid, endian='little')
            else:
                variables[name] = Variable(group, name, datatype, dimensions, id=varid)
        free(varids) # free pointer holding variable ids.
    return variables

cdef _ensure_nc_success(ierr, err_cls=RuntimeError, filename=None):
    # print netcdf error message, raise error.
    if ierr != NC_NOERR:
        err_str = (<char *>nc_strerror(ierr)).decode('ascii')
        if issubclass(err_cls, OSError):
            if isinstance(filename, bytes):
                filename = filename.decode()
            raise err_cls(ierr, err_str, filename)
        else:
            raise err_cls(err_str)

# these are class attributes that
# only exist at the python level (not in the netCDF file).

_private_atts = \
['_grpid','_grp','_varid','groups','dimensions','variables','dtype','data_model','disk_format',
 '_nunlimdim','path','parent','ndim','mask','scale','cmptypes','vltypes','enumtypes','_isprimitive',
 'file_format','_isvlen','_isenum','_iscompound','_cmptype','_vltype','_enumtype','name',
 '__orthogoral_indexing__','keepweakref','_has_lsd',
 '_buffer','chartostring','_use_get_vars','_ncstring_attrs__']

cdef class Dataset:
    """
A netCDF `Dataset` is a collection of dimensions, groups, variables and
attributes. Together they describe the meaning of data and relations among
data fields stored in a netCDF file. See `Dataset.__init__` for more
details.

A list of attribute names corresponding to global netCDF attributes
defined for the `Dataset` can be obtained with the
`Dataset.ncattrs` method.
These attributes can be created by assigning to an attribute of the
`Dataset` instance. A dictionary containing all the netCDF attribute
name/value pairs is provided by the `__dict__` attribute of a
`Dataset` instance.

The following class variables are read-only and should not be
modified by the user.

**`dimensions`**: The `dimensions` dictionary maps the names of
dimensions defined for the `Group` or `Dataset` to instances of the
`Dimension` class.

**`variables`**: The `variables` dictionary maps the names of variables
defined for this `Dataset` or `Group` to instances of the
`Variable` class.

**`groups`**: The groups dictionary maps the names of groups created for
this `Dataset` or `Group` to instances of the `Group` class (the
`Dataset` class is simply a special case of the `Group` class which
describes the root group in the netCDF4 file).

**`cmptypes`**: The `cmptypes` dictionary maps the names of
compound types defined for the `Group` or `Dataset` to instances of the
`CompoundType` class.

**`vltypes`**: The `vltypes` dictionary maps the names of
variable-length types defined for the `Group` or `Dataset` to instances
of the `VLType` class.

**`enumtypes`**: The `enumtypes` dictionary maps the names of
Enum types defined for the `Group` or `Dataset` to instances
of the `EnumType` class.

**`data_model`**: `data_model` describes the netCDF
data model version, one of `NETCDF3_CLASSIC`, `NETCDF4`,
`NETCDF4_CLASSIC`, `NETCDF3_64BIT_OFFSET` or `NETCDF3_64BIT_DATA`.

**`file_format`**: same as `data_model`, retained for backwards compatibility.

**`disk_format`**: `disk_format` describes the underlying
file format, one of `NETCDF3`, `HDF5`, `HDF4`,
`PNETCDF`, `DAP2`, `DAP4` or `UNDEFINED`. Only available if using
netcdf C library version >= 4.3.1, otherwise will always return
`UNDEFINED`.

**`parent`**: `parent` is a reference to the parent
`Group` instance. `None` for the root group or `Dataset`
instance.

**`path`**: `path` shows the location of the `Group` in
the `Dataset` in a unix directory format (the names of groups in the
hierarchy separated by backslashes). A `Dataset` instance is the root
group, so the path is simply `'/'`.

**`keepweakref`**: If `True`, child Dimension and Variables objects only keep weak
references to the parent Dataset or Group.

**`_ncstring_attrs__`**: If `True`, all text attributes will be written as variable-length
strings.
    """
    cdef object __weakref__, _inmemory
    cdef public int _grpid
    cdef public int _isopen
    cdef Py_buffer _buffer
    cdef public groups, dimensions, variables, disk_format, path, parent,\
    file_format, data_model, cmptypes, vltypes, enumtypes,  __orthogonal_indexing__, \
    keepweakref, _ncstring_attrs__

    def __init__(self, filename, mode='r', clobber=True, format='NETCDF4',
                     diskless=False, persist=False, keepweakref=False,
                     memory=None, encoding=None, parallel=False,
                     Comm comm=None, Info info=None, **kwargs):
        """
        **`__init__(self, filename, mode="r", clobber=True, diskless=False,
        persist=False, keepweakref=False, memory=None, encoding=None,
        parallel=False, comm=None, info=None, format='NETCDF4')`**

        `Dataset` constructor.

        **`filename`**: Name of netCDF file to hold dataset. Can also
        be a python 3 pathlib instance or the URL of an OpenDAP dataset.  When memory is
        set this is just used to set the `filepath()`.

        **`mode`**: access mode. `r` means read-only; no data can be
        modified. `w` means write; a new file is created, an existing file with
        the same name is deleted. `x` means write, but fail if an existing
        file with the same name already exists. `a` and `r+` mean append; 
        an existing file is opened for reading and writing, if 
        file does not exist already, one is created.
        Appending `s` to modes `r`, `w`, `r+` or `a` will enable unbuffered shared
        access to `NETCDF3_CLASSIC`, `NETCDF3_64BIT_OFFSET` or
        `NETCDF3_64BIT_DATA` formatted files.
        Unbuffered access may be useful even if you don't need shared
        access, since it may be faster for programs that don't access data
        sequentially. This option is ignored for `NETCDF4` and `NETCDF4_CLASSIC`
        formatted files.

        **`clobber`**: if `True` (default), opening a file with `mode='w'`
        will clobber an existing file with the same name.  if `False`, an
        exception will be raised if a file with the same name already exists.
        mode=`x` is identical to mode=`w` with clobber=False.

        **`format`**: underlying file format (one of `'NETCDF4',
        'NETCDF4_CLASSIC', 'NETCDF3_CLASSIC'`, `'NETCDF3_64BIT_OFFSET'` or
        `'NETCDF3_64BIT_DATA'`.
        Only relevant if `mode = 'w'` (if `mode = 'r','a'` or `'r+'` the file format
        is automatically detected). Default `'NETCDF4'`, which means the data is
        stored in an HDF5 file, using netCDF 4 API features.  Setting
        `format='NETCDF4_CLASSIC'` will create an HDF5 file, using only netCDF 3
        compatible API features. netCDF 3 clients must be recompiled and linked
        against the netCDF 4 library to read files in `NETCDF4_CLASSIC` format.
        `'NETCDF3_CLASSIC'` is the classic netCDF 3 file format that does not
        handle 2+ Gb files. `'NETCDF3_64BIT_OFFSET'` is the 64-bit offset
        version of the netCDF 3 file format, which fully supports 2+ GB files, but
        is only compatible with clients linked against netCDF version 3.6.0 or
        later. `'NETCDF3_64BIT_DATA'` is the 64-bit data version of the netCDF 3
        file format, which supports 64-bit dimension sizes plus unsigned and
        64 bit integer data types, but is only compatible with clients linked against
        netCDF version 4.4.0 or later.

        **`diskless`**: If `True`, create diskless (in-core) file.
        This is a feature added to the C library after the
        netcdf-4.2 release. If you need to access the memory buffer directly,
        use the in-memory feature instead (see `memory` kwarg).

        **`persist`**: if `diskless=True`, persist file to disk when closed
        (default `False`).

        **`keepweakref`**: if `True`, child Dimension and Variable instances will keep weak
        references to the parent Dataset or Group object.  Default is `False`, which
        means strong references will be kept.  Having Dimension and Variable instances
        keep a strong reference to the parent Dataset instance, which in turn keeps a
        reference to child Dimension and Variable instances, creates circular references.
        Circular references complicate garbage collection, which may mean increased
        memory usage for programs that create may Dataset instances with lots of
        Variables. It also will result in the Dataset object never being deleted, which
        means it may keep open files alive as well. Setting `keepweakref=True` allows
        Dataset instances to be garbage collected as soon as they go out of scope, potentially
        reducing memory usage and open file handles.  However, in many cases this is not
        desirable, since the associated Variable instances may still be needed, but are
        rendered unusable when the parent Dataset instance is garbage collected.

        **`memory`**: if not `None`, create or open an in-memory Dataset.
        If mode = `r`, the memory kwarg must contain a memory buffer object
        (an object that supports the python buffer interface).
        The Dataset will then be created with contents taken from this block of memory.
        If mode = `w`, the memory kwarg should contain the anticipated size
        of the Dataset in bytes (used only for NETCDF3 files).  A memory
        buffer containing a copy of the Dataset is returned by the
        `Dataset.close` method. Requires netcdf-c version 4.4.1 for mode=`r`
        netcdf-c 4.6.2 for mode=`w`. To persist the file to disk, the raw
        bytes from the returned buffer can be written into a binary file.
        The Dataset can also be re-opened using this memory buffer.

        **`encoding`**: encoding used to encode filename string into bytes.
        Default is None (`sys.getdefaultfileencoding()` is used).

        **`parallel`**: open for parallel access using MPI (requires mpi4py and
        parallel-enabled netcdf-c and hdf5 libraries).  Default is `False`. If
        `True`, `comm` and `info` kwargs may also be specified.

        **`comm`**: MPI_Comm object for parallel access. Default `None`, which
        means MPI_COMM_WORLD will be used.  Ignored if `parallel=False`.

        **`info`**: MPI_Info object for parallel access. Default `None`, which
        means MPI_INFO_NULL will be used.  Ignored if `parallel=False`.
        """
        cdef int grpid, ierr, numgrps, numdims, numvars,
        cdef size_t initialsize
        cdef char *path
        cdef char namstring[NC_MAX_NAME+1]
        cdef int cmode, parmode
        IF HAS_PARALLEL4_SUPPORT or HAS_PNETCDF_SUPPORT:
            cdef MPI_Comm mpicomm
            cdef MPI_Info mpiinfo

        memset(&self._buffer, 0, sizeof(self._buffer))

        # flag to indicate that Variables in this Dataset support orthogonal indexing.
        self.__orthogonal_indexing__ = True
        if diskless and __netcdf4libversion__ < '4.2.1':
            #diskless = False # don't raise error, instead silently ignore
            raise ValueError('diskless mode requires netcdf lib >= 4.2.1, you have %s' % __netcdf4libversion__)
        # convert filename into string (from os.path object for example),
        # encode into bytes.
        if encoding is None:
            encoding = sys.getfilesystemencoding()
        bytestr = _strencode(_tostr(filename), encoding=encoding)
        path = bytestr

        if memory is not None and mode not in ['r','w']:
            msg='if memory kwarg specified, mode must be \'r\' or \'w\''
            raise ValueError(msg)

        if parallel:
            IF HAS_PARALLEL4_SUPPORT != 1 and HAS_PNETCDF_SUPPORT != 1:
                msg='parallel mode requires MPI enabled netcdf-c'
                raise ValueError(msg)
            ELSE:
                parallel_formats = []
                IF HAS_PARALLEL4_SUPPORT:
                    parallel_formats += ['NETCDF4','NETCDF4_CLASSIC']
                IF HAS_PNETCDF_SUPPORT:
                    parallel_formats += ['NETCDF3_CLASSIC',
                                         'NETCDF3_64BIT_OFFSET',
                                         'NETCDF3_64BIT_DATA',
                                         'NETCDF3_64BIT']
                if format not in parallel_formats:
                    msg='parallel mode only works with the following formats: ' + ' '.join(parallel_formats)
                    raise ValueError(msg)
                if comm is not None:
                    mpicomm = comm.ob_mpi
                else:
                    mpicomm = MPI_COMM_WORLD
                if info is not None:
                    mpiinfo = info.ob_mpi
                else:
                    mpiinfo = MPI_INFO_NULL
                parmode = NC_MPIIO | _cmode_dict[format]

        self._inmemory = False

        # mode='x' is the same as mode='w' with clobber=False
        if mode == 'x':
            mode = 'w'; clobber = False

        if mode == 'w' or (mode in ['a','r+'] and not os.path.exists(filename)):
            _set_default_format(format=format)
            if memory is not None:
                # if memory is not None and mode='w', memory
                # kwarg is interpreted as advisory size.
                IF HAS_NC_CREATE_MEM:
                   initialsize = <size_t>memory
                   with nogil:
                       ierr = nc_create_mem(path, 0, initialsize, &grpid)
                   self._inmemory = True # checked in close method
                ELSE:
                    msg = """
        nc_create_mem functionality not enabled.  To enable, install Cython, make sure you have
        version 4.6.2 or higher of the netcdf C lib, and rebuild netcdf4-python."""
                    raise ValueError(msg)
            else:
                if clobber:
                    if parallel:
                        IF HAS_PARALLEL4_SUPPORT or HAS_PNETCDF_SUPPORT:
                            cmode = NC_CLOBBER | parmode
                            with nogil:
                                ierr = nc_create_par(path, cmode, \
                                       mpicomm, mpiinfo, &grpid)
                        ELSE:
                            pass
                    elif diskless:
                        if persist:
                            cmode = NC_WRITE | NC_CLOBBER | NC_DISKLESS | NC_PERSIST
                            with nogil:
                                ierr = nc_create(path, cmode, &grpid)
                        else:
                            cmode = NC_CLOBBER | NC_DISKLESS
                            with nogil:
                                ierr = nc_create(path, cmode , &grpid)
                    else:
                        with nogil:
                            ierr = nc_create(path, NC_CLOBBER, &grpid)
                else:
                    if parallel:
                        IF HAS_PARALLEL4_SUPPORT or HAS_PNETCDF_SUPPORT:
                            cmode = NC_NOCLOBBER | parmode
                            with nogil:
                                ierr = nc_create_par(path, cmode, \
                                       mpicomm, mpiinfo, &grpid)
                        ELSE:
                            pass
                    elif diskless:
                        if persist:
                            cmode = NC_WRITE | NC_NOCLOBBER | NC_DISKLESS | NC_PERSIST
                            with nogil:
                                ierr = nc_create(path, cmode, &grpid)
                        else:
                            cmode = NC_NOCLOBBER | NC_DISKLESS
                            with nogil:
                                ierr = nc_create(path, cmode , &grpid)
                    else:
                        with nogil:
                            ierr = nc_create(path, NC_NOCLOBBER, &grpid)
            # reset default format to netcdf3 - this is a workaround
            # for issue 170 (nc_open'ing a DAP dataset after switching
            # format to NETCDF4). This bug should be fixed in version
            # 4.3.0 of the netcdf library (add a version check here?).
            # **this causes parallel mode to fail when both hdf5-parallel and
            # pnetcdf are enabled - issue #820 **
            #_set_default_format(format='NETCDF3_64BIT_OFFSET')
        elif mode in ('r', 'rs'):
            if memory is not None:
                IF HAS_NC_OPEN_MEM:
                    # Store reference to memory
                    result = PyObject_GetBuffer(memory, &self._buffer, PyBUF_SIMPLE | PyBUF_ANY_CONTIGUOUS)
                    if result != 0:
                        raise ValueError("Unable to retrieve Buffer from %s" % (memory,))

                    with nogil:
                        ierr = nc_open_mem(<char *>path, 0, self._buffer.len, <void *>self._buffer.buf, &grpid)
                ELSE:
                    msg = """
        nc_open_mem functionality not enabled.  To enable, install Cython, make sure you have
        version 4.4.1 or higher of the netcdf C lib, and rebuild netcdf4-python."""
                    raise ValueError(msg)
            elif parallel:
                IF HAS_PARALLEL4_SUPPORT or HAS_PNETCDF_SUPPORT:
                    cmode = NC_NOWRITE | NC_MPIIO
                    with nogil:
                        ierr = nc_open_par(path, cmode, \
                               mpicomm, mpiinfo, &grpid)
                ELSE:
                    pass
            elif diskless:
                cmode = NC_NOWRITE | NC_DISKLESS
                with nogil:
                    ierr = nc_open(path, cmode, &grpid)
            else:
                if mode == 'rs':
                    # NC_SHARE is very important for speed reading
                    # large netcdf3 files with a record dimension
                    # (pull request #902).
                    cmode = NC_NOWRITE | NC_SHARE
                    with nogil:
                        ierr = nc_open(path, cmode, &grpid)
                else:
                    with nogil:
                        ierr = nc_open(path, NC_NOWRITE, &grpid)
        elif mode in ['a','r+'] and os.path.exists(filename):
            if parallel:
                IF HAS_PARALLEL4_SUPPORT or HAS_PNETCDF_SUPPORT:
                    cmode = NC_WRITE | NC_MPIIO
                    with nogil:
                        ierr = nc_open_par(path, cmode, \
                               mpicomm, mpiinfo, &grpid)
                ELSE:
                    pass
            elif diskless:
                cmode = NC_WRITE | NC_DISKLESS
                with nogil:
                    ierr = nc_open(path, cmode, &grpid)
            else:
                with nogil:
                    ierr = nc_open(path, NC_WRITE, &grpid)
        elif mode in ['as','r+s'] and os.path.exists(filename):
            if parallel:
                # NC_SHARE ignored
                IF HAS_PARALLEL4_SUPPORT or HAS_PNETCDF_SUPPORT:
                    cmode =  NC_WRITE | NC_MPIIO
                    with nogil:
                        ierr = nc_open_par(path, cmode, \
                               mpicomm, mpiinfo, &grpid)
                ELSE:
                    pass
            elif diskless:
                cmode = NC_SHARE | NC_DISKLESS
                with nogil:
                    ierr = nc_open(path, cmode, &grpid)
            else:
                with nogil:
                    ierr = nc_open(path, NC_SHARE, &grpid)
        elif mode == 'ws' or (mode in ['as','r+s'] and not os.path.exists(filename)):
            _set_default_format(format=format)
            if clobber:
                if parallel:
                    # NC_SHARE ignored
                    IF HAS_PARALLEL4_SUPPORT or HAS_PNETCDF_SUPPORT:
                        cmode = NC_CLOBBER | parmode 
                        with nogil:
                            ierr = nc_create_par(path, NC_CLOBBER | cmode, \
                                   mpicomm, mpiinfo, &grpid)
                    ELSE:
                        pass
                elif diskless:
                    if persist:
                        cmode = NC_WRITE | NC_SHARE | NC_CLOBBER | NC_DISKLESS
                        with nogil:
                            ierr = nc_create(path, cmode, &grpid)
                    else:
                        cmode = NC_SHARE | NC_CLOBBER | NC_DISKLESS
                        with nogil:
                            ierr = nc_create(path, cmode , &grpid)
                else:
                    cmode = NC_SHARE | NC_CLOBBER
                    with nogil:
                        ierr = nc_create(path, cmode, &grpid)
            else:
                if parallel:
                    # NC_SHARE ignored
                    IF HAS_PARALLEL4_SUPPORT or HAS_PNETCDF_SUPPORT:
                        cmode = NC_NOCLOBBER | parmode
                        with nogil:
                            ierr = nc_create_par(path, cmode, \
                                   mpicomm, mpiinfo, &grpid)
                    ELSE:
                        pass
                elif diskless:
                    if persist:
                        cmode = NC_WRITE | NC_SHARE | NC_NOCLOBBER | NC_DISKLESS
                        with nogil:
                            ierr = nc_create(path, cmode , &grpid)
                    else:
                        cmode = NC_SHARE | NC_NOCLOBBER | NC_DISKLESS
                        with nogil:
                            ierr = nc_create(path, cmode , &grpid)
                else:
                    cmode = NC_SHARE | NC_NOCLOBBER
                    with nogil:
                        ierr = nc_create(path, cmode, &grpid)
        else:
            raise ValueError("mode must be 'w', 'x', 'r', 'a' or 'r+', got '%s'" % mode)

        _ensure_nc_success(ierr, err_cls=OSError, filename=path)

        # data model and file format attributes
        self.data_model = _get_format(grpid)
        # data_model attribute used to be file_format (versions < 1.0.8), retain
        # file_format for backwards compatibility.
        self.file_format = self.data_model
        self.disk_format = _get_full_format(grpid)
        # diskless read access only works with NETCDF_CLASSIC (for now)
        #ncopen = mode.startswith('a') or mode.startswith('r')
        #if diskless and self.data_model != 'NETCDF3_CLASSIC' and ncopen:
        #    raise ValueError("diskless access only supported for NETCDF3_CLASSIC format")
        self._grpid = grpid
        self._isopen = 1
        self.path = '/'
        self.parent = None
        self.keepweakref = keepweakref
        self._ncstring_attrs__ = False
        # get compound, vlen and enum types in the root Group.
        self.cmptypes, self.vltypes, self.enumtypes = _get_types(self)
        # get dimensions in the root group.
        self.dimensions = _get_dims(self)
        # get variables in the root Group.
        self.variables = _get_vars(self)
        # get groups in the root Group.
        if self.data_model == 'NETCDF4':
            self.groups = _get_grps(self)
        else:
            self.groups = dict()

    # these allow Dataset objects to be used via a "with" statement.
    def __enter__(self):
        return self
    def __exit__(self,atype,value,traceback):
        self.close()

    def __getitem__(self, elem):
        # return variable or group defined in relative path.
        # split out group names in unix path.
        elem = posixpath.normpath(elem)
        # last name in path, could be a variable or group
        dirname, lastname = posixpath.split(elem)
        nestedgroups = dirname.split('/')
        group = self
        # iterate over groups in path.
        for g in nestedgroups:
            if g: group = group.groups[g]
        # return last one, either a group or a variable.
        if lastname in group.groups:
            return group.groups[lastname]
        elif lastname in group.variables:
            return group.variables[lastname]
        else:
            raise IndexError('%s not found in %s' % (lastname,group.path))

    def filepath(self,encoding=None):
        """
**`filepath(self,encoding=None)`**

Get the file system path (or the opendap URL) which was used to
open/create the Dataset. Requires netcdf >= 4.1.2.  The path
is decoded into a string using `sys.getfilesystemencoding()` by default, this can be
changed using the `encoding` kwarg."""
        cdef int ierr
        cdef size_t pathlen
        cdef char *c_path
        if encoding is None:
            encoding = sys.getfilesystemencoding()
        IF HAS_NC_INQ_PATH:
            with nogil:
                ierr = nc_inq_path(self._grpid, &pathlen, NULL)
            _ensure_nc_success(ierr)

            c_path = <char *>malloc(sizeof(char) * (pathlen + 1))
            if not c_path:
                raise MemoryError()
            try:
                with nogil:
                    ierr = nc_inq_path(self._grpid, &pathlen, c_path)
                _ensure_nc_success(ierr)

                py_path = c_path[:pathlen] # makes a copy of pathlen bytes from c_string
            finally:
                free(c_path)
            return py_path.decode(encoding)
        ELSE:
            msg = """
filepath method not enabled.  To enable, install Cython, make sure you have
version 4.1.2 or higher of the netcdf C lib, and rebuild netcdf4-python."""
            raise ValueError(msg)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        ncdump = [repr(type(self))]
        dimnames = tuple(_tostr(dimname)+'(%s)'%len(self.dimensions[dimname])\
        for dimname in self.dimensions.keys())
        varnames = tuple(\
        [_tostr(self.variables[varname].dtype)+' '+_tostr(varname)+
         ((_tostr(self.variables[varname].dimensions)).replace(",)",")")).replace("'","")
         for varname in self.variables.keys()])
        grpnames = tuple(_tostr(grpname) for grpname in self.groups.keys())
        if self.path == '/':
            ncdump.append('root group (%s data model, file format %s):' %
                    (self.data_model, self.disk_format))
        else:
            ncdump.append('group %s:' % self.path)
        for name in self.ncattrs():
            ncdump.append('    %s: %s' % (name, self.getncattr(name)))
        ncdump.append('    dimensions(sizes): %s' % ', '.join(dimnames))
        ncdump.append('    variables(dimensions): %s' % ', '.join(varnames))
        ncdump.append('    groups: %s' % ', '.join(grpnames))
        return '\n'.join(ncdump)

    def _close(self, check_err):
        cdef int ierr
        with nogil:
            ierr = nc_close(self._grpid)

        if check_err:
            _ensure_nc_success(ierr)

        self._isopen = 0 # indicates file already closed, checked by __dealloc__

        # Only release buffer if close succeeded
        # per impl of PyBuffer_Release: https://github.com/python/cpython/blob/master/Objects/abstract.c#L667
        # view.obj is checked, ref on obj is decremented and obj will be null'd out
        PyBuffer_Release(&self._buffer)

    IF HAS_NC_CREATE_MEM:
        def _close_mem(self, check_err):
            cdef int ierr
            cdef NC_memio memio
            with nogil:
                ierr = nc_close_memio(self._grpid, &memio)

            if check_err:
                _ensure_nc_success(ierr)

            self._isopen = 0
            PyBuffer_Release(&self._buffer)

            # membuf_fromptr from membuf.pyx - creates a python memoryview
            # from a raw pointer without making a copy.
            return memview_fromptr(<char *>memio.memory, memio.size)


    def close(self):
        """
**`close(self)`**

Close the Dataset.
        """
        IF HAS_NC_CREATE_MEM:
            if self._inmemory:
                return self._close_mem(True)
            else:
                self._close(True)
        ELSE:
            self._close(True)

    def isopen(self):
        """
**`isopen(self)`**

Is the Dataset open or closed?
        """
        return bool(self._isopen)

    def __dealloc__(self):
        # close file when there are no references to object left
        if self._isopen:
           self._close(False)

    def __reduce__(self):
        # raise error is user tries to pickle a Dataset object.
        raise NotImplementedError('Dataset is not picklable')

    def sync(self):
        """
**`sync(self)`**

Writes all buffered data in the `Dataset` to the disk file."""
        cdef int ierr
        with nogil:
            ierr = nc_sync(self._grpid)
        _ensure_nc_success(ierr)

    def _redef(self):
        cdef int ierr
        with nogil:
            ierr = nc_redef(self._grpid)

    def _enddef(self):
        cdef int ierr
        with nogil:
            ierr = nc_enddef(self._grpid)

    def set_fill_on(self):
        """
**`set_fill_on(self)`**

Sets the fill mode for a `Dataset` open for writing to `on`.

This causes data to be pre-filled with fill values. The fill values can be
controlled by the variable's `_Fill_Value` attribute, but is usually
sufficient to the use the netCDF default `_Fill_Value` (defined
separately for each variable type). The default behavior of the netCDF
library corresponds to `set_fill_on`.  Data which are equal to the
`_Fill_Value` indicate that the variable was created, but never written
to."""
        cdef int oldmode, ierr
        with nogil:
            ierr = nc_set_fill(self._grpid, NC_FILL, &oldmode)
        _ensure_nc_success(ierr)

    def set_fill_off(self):
        """
**`set_fill_off(self)`**

Sets the fill mode for a `Dataset` open for writing to `off`.

This will prevent the data from being pre-filled with fill values, which
may result in some performance improvements. However, you must then make
sure the data is actually written before being read."""
        cdef int oldmode, ierr
        with nogil:
            ierr = nc_set_fill(self._grpid, NC_NOFILL, &oldmode)
        _ensure_nc_success(ierr)

    def createDimension(self, dimname, size=None):
        """
**`createDimension(self, dimname, size=None)`**

Creates a new dimension with the given `dimname` and `size`.

`size` must be a positive integer or `None`, which stands for
"unlimited" (default is `None`). Specifying a size of 0 also
results in an unlimited dimension. The return value is the `Dimension`
class instance describing the new dimension.  To determine the current
maximum size of the dimension, use the `len` function on the `Dimension`
instance. To determine if a dimension is 'unlimited', use the
`Dimension.isunlimited` method of the `Dimension` instance."""
        self.dimensions[dimname] = Dimension(self, dimname, size=size)
        return self.dimensions[dimname]

    def renameDimension(self, oldname, newname):
        """
**`renameDimension(self, oldname, newname)`**

rename a `Dimension` named `oldname` to `newname`."""
        cdef char *namstring
        cdef Dimension dim
        bytestr = _strencode(newname)
        namstring = bytestr
        if self.data_model != 'NETCDF4': self._redef()
        try:
            dim = self.dimensions[oldname]
        except KeyError:
            raise KeyError('%s not a valid dimension name' % oldname)
        with nogil:
            ierr = nc_rename_dim(self._grpid, dim._dimid, namstring)
        if self.data_model != 'NETCDF4': self._enddef()

        _ensure_nc_success(ierr)
        # remove old key from dimensions dict.
        self.dimensions.pop(oldname)
        # add new key.
        self.dimensions[newname] = dim
        # Variable.dimensions is determined by a method that
        # looks in the file, so no need to manually update.

    def createCompoundType(self, datatype, datatype_name):
        """
**`createCompoundType(self, datatype, datatype_name)`**

Creates a new compound data type named `datatype_name` from the numpy
dtype object `datatype`.

***Note***: If the new compound data type contains other compound data types
(i.e. it is a 'nested' compound type, where not all of the elements
are homogeneous numeric data types), then the 'inner' compound types **must** be
created first.

The return value is the `CompoundType` class instance describing the new
datatype."""
        self.cmptypes[datatype_name] = CompoundType(self, datatype,\
                datatype_name)
        return self.cmptypes[datatype_name]

    def createVLType(self, datatype, datatype_name):
        """
**`createVLType(self, datatype, datatype_name)`**

Creates a new VLEN data type named `datatype_name` from a numpy
dtype object `datatype`.

The return value is the `VLType` class instance describing the new
datatype."""
        self.vltypes[datatype_name] = VLType(self, datatype, datatype_name)
        return self.vltypes[datatype_name]

    def createEnumType(self, datatype, datatype_name, enum_dict):
        """
**`createEnumType(self, datatype, datatype_name, enum_dict)`**

Creates a new Enum data type named `datatype_name` from a numpy
integer dtype object `datatype`, and a python dictionary
defining the enum fields and values.

The return value is the `EnumType` class instance describing the new
datatype."""
        self.enumtypes[datatype_name] = EnumType(self, datatype, datatype_name,
                enum_dict)
        return self.enumtypes[datatype_name]

    def createVariable(self, varname, datatype, dimensions=(), 
            compression=None, zlib=False,
            complevel=4, shuffle=True,
            szip_coding='nn',szip_pixels_per_block=8,
            blosc_shuffle=1,fletcher32=False, contiguous=False,
            chunksizes=None, endian='native', least_significant_digit=None,
            significant_digits=None,quantize_mode='BitGroom',fill_value=None, chunk_cache=None):
        """
**`createVariable(self, varname, datatype, dimensions=(), compression=None, zlib=False,
complevel=4, shuffle=True, fletcher32=False, contiguous=False, chunksizes=None,
szip_coding='nn', szip_pixels_per_block=8, blosc_shuffle=1,
endian='native', least_significant_digit=None, significant_digits=None, quantize_mode='BitGroom',
fill_value=None, chunk_cache=None)`**

Creates a new variable with the given `varname`, `datatype`, and
`dimensions`. If dimensions are not given, the variable is assumed to be
a scalar.

If `varname` is specified as a path, using forward slashes as in unix to
separate components, then intermediate groups will be created as necessary
For example, `createVariable('/GroupA/GroupB/VarC', float, ('x','y'))` will create groups `GroupA`
and `GroupA/GroupB`, plus the variable `GroupA/GroupB/VarC`, if the preceding
groups don't already exist.

The `datatype` can be a numpy datatype object, or a string that describes
a numpy dtype object (like the `dtype.str` attribute of a numpy array).
Supported specifiers include: `'S1' or 'c' (NC_CHAR), 'i1' or 'b' or 'B'
(NC_BYTE), 'u1' (NC_UBYTE), 'i2' or 'h' or 's' (NC_SHORT), 'u2'
(NC_USHORT), 'i4' or 'i' or 'l' (NC_INT), 'u4' (NC_UINT), 'i8' (NC_INT64),
'u8' (NC_UINT64), 'f4' or 'f' (NC_FLOAT), 'f8' or 'd' (NC_DOUBLE)`.
`datatype` can also be a `CompoundType` instance
(for a structured, or compound array), a `VLType` instance
(for a variable-length array), or the python `str` builtin
(for a variable-length string array). Numpy string and unicode datatypes with
length greater than one are aliases for `str`.

Data from netCDF variables is presented to python as numpy arrays with
the corresponding data type.

`dimensions` must be a tuple containing `Dimension` instances and/or
dimension names (strings) that have been defined
previously using `Dataset.createDimension`. The default value
is an empty tuple, which means the variable is a scalar.

If the optional keyword argument `compression` is set, the data will be
compressed in the netCDF file using the specified compression algorithm.
Currently `zlib`,`szip`,`zstd`,`bzip2`,`blosc_lz`,`blosc_lz4`,`blosc_lz4hc`,
`blosc_zlib` and `blosc_zstd` are supported.
Default is `None` (no compression).  All of the compressors except
`zlib` and `szip` use the HDF5 plugin architecture.

If the optional keyword `zlib` is `True`, the data will be compressed in
the netCDF file using zlib compression (default `False`).  The use of this option is 
deprecated in favor of `compression='zlib'`.

The optional keyword `complevel` is an integer between 0 and 9 describing
the level of compression desired (default 4). Ignored if `compression=None`.
A value of zero disables compression.

If the optional keyword `shuffle` is `True`, the HDF5 shuffle filter
will be applied before compressing the data with zlib (default `True`).  This
significantly improves compression. Default is `True`. Ignored if
`zlib=False`.

The optional kwarg `blosc_shuffle`is  ignored
unless the blosc compressor is used. `blosc_shuffle` can be 0 (no shuffle),
1 (byte-wise shuffle) or 2 (bit-wise shuffle). Default is 1.

The optional kwargs `szip_coding` and `szip_pixels_per_block` are ignored
unless the szip compressor is used. `szip_coding` can be `ec` (entropy coding)
or `nn` (nearest neighbor coding). Default is `nn`. `szip_pixels_per_block`
can be 4, 8, 16 or 32 (default 8).

If the optional keyword `fletcher32` is `True`, the Fletcher32 HDF5
checksum algorithm is activated to detect errors. Default `False`.

If the optional keyword `contiguous` is `True`, the variable data is
stored contiguously on disk.  Default `False`. Setting to `True` for
a variable with an unlimited dimension will trigger an error.
Fixed size variables (with no unlimited dimension) with no compression filters
are contiguous by default.

The optional keyword `chunksizes` can be used to manually specify the
HDF5 chunksizes for each dimension of the variable.
A detailed discussion of HDF chunking and I/O performance is available
[here](https://support.hdfgroup.org/HDF5/doc/Advanced/Chunking).
The default chunking scheme in the netcdf-c library is discussed
[here](https://www.unidata.ucar.edu/software/netcdf/documentation/NUG/netcdf_perf_chunking.html).
Basically, you want the chunk size for each dimension to match as
closely as possible the size of the data block that users will read
from the file. `chunksizes` cannot be set if `contiguous=True`.

The optional keyword `endian` can be used to control whether the
data is stored in little or big endian format on disk. Possible
values are `little, big` or `native` (default). The library
will automatically handle endian conversions when the data is read,
but if the data is always going to be read on a computer with the
opposite format as the one used to create the file, there may be
some performance advantage to be gained by setting the endian-ness.

The optional keyword `fill_value` can be used to override the default
netCDF `_FillValue` (the value that the variable gets filled with before
any data is written to it, defaults given in the dict `netCDF4.default_fillvals`).
If fill_value is set to `False`, then the variable is not pre-filled.

If the optional keyword parameters `least_significant_digit` or `significant_digits` are
specified, variable data will be truncated (quantized). In conjunction
with `compression='zlib'` this produces 'lossy', but significantly more
efficient compression. For example, if `least_significant_digit=1`,
data will be quantized using `numpy.around(scale*data)/scale`, where
scale = 2**bits, and bits is determined so that a precision of 0.1 is
retained (in this case bits=4). From the
[PSL metadata conventions](http://www.esrl.noaa.gov/psl/data/gridded/conventions/cdc_netcdf_standard.shtml):
"least_significant_digit -- power of ten of the smallest decimal place
in unpacked data that is a reliable value." Default is `None`, or no
quantization, or 'lossless' compression.  If `significant_digits=3`
then the data will be quantized so that three significant digits are retained, independent
of the floating point exponent. The keyword argument `quantize_mode` controls
the quantization algorithm (default 'BitGroom', 'BitRound' and
'GranularBitRound' also available).  The 'GranularBitRound'
algorithm may result in better compression for typical geophysical datasets.
This `significant_digits` kwarg is only available  with netcdf-c >= 4.9.0, and
only works with `NETCDF4` or `NETCDF4_CLASSIC` formatted files.

When creating variables in a `NETCDF4` or `NETCDF4_CLASSIC` formatted file,
HDF5 creates something called a 'chunk cache' for each variable.  The
default size of the chunk cache may be large enough to completely fill
available memory when creating thousands of variables.  The optional
keyword `chunk_cache` allows you to reduce (or increase) the size of
the default chunk cache when creating a variable.  The setting only
persists as long as the Dataset is open - you can use the set_var_chunk_cache
method to change it the next time the Dataset is opened.
Warning - messing with this parameter can seriously degrade performance.

The return value is the `Variable` class instance describing the new
variable.

A list of names corresponding to netCDF variable attributes can be
obtained with the `Variable` method `Variable.ncattrs`. A dictionary
containing all the netCDF attribute name/value pairs is provided by
the `__dict__` attribute of a `Variable` instance.

`Variable` instances behave much like array objects. Data can be
assigned to or retrieved from a variable with indexing and slicing
operations on the `Variable` instance. A `Variable` instance has six
Dataset standard attributes: `dimensions, dtype, shape, ndim, name` and
`least_significant_digit`. Application programs should never modify
these attributes. The `dimensions` attribute is a tuple containing the
names of the dimensions associated with this variable. The `dtype`
attribute is a string describing the variable's data type (`i4, f8,
S1,` etc). The `shape` attribute is a tuple describing the current
sizes of all the variable's dimensions. The `name` attribute is a
string containing the name of the Variable instance.
The `least_significant_digit`
attributes describes the power of ten of the smallest decimal place in
the data the contains a reliable value.  assigned to the `Variable`
instance. The `ndim` attribute
is the number of variable dimensions."""
        # if varname specified as a path, split out group names.
        varname = posixpath.normpath(varname)
        dirname, varname = posixpath.split(varname) # varname is last.
        # create parent groups (like mkdir -p).
        if not dirname:
            group = self
        else:
            group = self.createGroup(dirname)
        # if dimensions is a single string or Dimension instance,
        # convert to a tuple.
        # This prevents a common error that occurs when
        # dimensions = 'lat' instead of ('lat',)
        if isinstance(dimensions, (str, bytes, Dimension)):
            dimensions = dimensions,
        # convert elements of dimensions tuple to Dimension
        # instances if they are strings.
        # _find_dim looks for dimension in this group, and if not
        # found there, looks in parent (and it's parent, etc, back to root).
        dimensions =\
        tuple(_find_dim(group,d) if isinstance(d,(str,bytes)) else d for d in dimensions)
        # create variable.
        group.variables[varname] = Variable(group, varname, datatype,
        dimensions=dimensions, compression=compression, zlib=zlib, complevel=complevel, shuffle=shuffle,
        szip_coding=szip_coding, szip_pixels_per_block=szip_pixels_per_block,
        blosc_shuffle=blosc_shuffle,
        fletcher32=fletcher32, contiguous=contiguous, chunksizes=chunksizes,
        endian=endian, least_significant_digit=least_significant_digit,
        significant_digits=significant_digits,quantize_mode=quantize_mode,fill_value=fill_value, chunk_cache=chunk_cache)
        return group.variables[varname]

    def renameVariable(self, oldname, newname):
        """
**`renameVariable(self, oldname, newname)`**

rename a `Variable` named `oldname` to `newname`"""
        cdef char *namstring
        cdef Variable var
        try:
            var = self.variables[oldname]
        except KeyError:
            raise KeyError('%s not a valid variable name' % oldname)
        bytestr = _strencode(newname)
        namstring = bytestr
        if self.data_model != 'NETCDF4': self._redef()
        with nogil:
            ierr = nc_rename_var(self._grpid, var._varid, namstring)
        if self.data_model != 'NETCDF4': self._enddef()

        _ensure_nc_success(ierr)
        # remove old key from dimensions dict.
        self.variables.pop(oldname)
        # add new key.
        self.variables[newname] = var

    def createGroup(self, groupname):
        """
**`createGroup(self, groupname)`**

Creates a new `Group` with the given `groupname`.

If `groupname` is specified as a path, using forward slashes as in unix to
separate components, then intermediate groups will be created as necessary
(analogous to `mkdir -p` in unix).  For example,
`createGroup('/GroupA/GroupB/GroupC')` will create `GroupA`,
`GroupA/GroupB`, and `GroupA/GroupB/GroupC`, if they don't already exist.
If the specified path describes a group that already exists, no error is
raised.

The return value is a `Group` class instance."""
        # if group specified as a path, split out group names
        groupname = posixpath.normpath(groupname)
        nestedgroups = groupname.split('/')
        group = self
        # loop over group names, create parent groups if they do not already
        # exist.
        for g in nestedgroups:
            if not g: continue
            if g not in group.groups:
                group.groups[g] = Group(group, g)
            group = group.groups[g]
        # if group already exists, just return the group
        # (prior to 1.1.8, this would have raised an error)
        return group

    def ncattrs(self):
        """
**`ncattrs(self)`**

return netCDF global attribute names for this `Dataset` or `Group` in a list."""
        return _get_att_names(self._grpid, NC_GLOBAL)

    def setncattr(self,name,value):
        """
**`setncattr(self,name,value)`**

set a netCDF dataset or group attribute using name,value pair.
Use if you need to set a netCDF attribute with the
with the same name as one of the reserved python attributes."""
        cdef nc_type xtype
        xtype=-99
        if self.data_model != 'NETCDF4': self._redef()
        _set_att(self, NC_GLOBAL, name, value, xtype=xtype, force_ncstring=self._ncstring_attrs__)
        if self.data_model !=  'NETCDF4': self._enddef()

    def setncattr_string(self,name,value):
        """
**`setncattr_string(self,name,value)`**

set a netCDF dataset or group string attribute using name,value pair.
Use if you need to ensure that a netCDF attribute is created with type
`NC_STRING` if the file format is `NETCDF4`."""
        cdef nc_type xtype
        xtype=-99
        if self.data_model != 'NETCDF4':
            msg='file format does not support NC_STRING attributes'
            raise OSError(msg)
        _set_att(self, NC_GLOBAL, name, value, xtype=xtype, force_ncstring=True)

    def setncatts(self,attdict):
        """
**`setncatts(self,attdict)`**

set a bunch of netCDF dataset or group attributes at once using a python dictionary.
This may be faster when setting a lot of attributes for a `NETCDF3`
formatted file, since nc_redef/nc_enddef is not called in between setting
each attribute"""
        if self.data_model != 'NETCDF4': self._redef()
        for name, value in attdict.items():
            _set_att(self, NC_GLOBAL, name, value)
        if self.data_model != 'NETCDF4': self._enddef()

    def getncattr(self,name,encoding='utf-8'):
        """
**`getncattr(self,name)`**

retrieve a netCDF dataset or group attribute.
Use if you need to get a netCDF attribute with the same
name as one of the reserved python attributes.

option kwarg `encoding` can be used to specify the
character encoding of a string attribute (default is `utf-8`)."""
        return _get_att(self, NC_GLOBAL, name, encoding=encoding)

    def __delattr__(self,name):
        # if it's a netCDF attribute, remove it
        if name not in _private_atts:
            self.delncattr(name)
        else:
            raise AttributeError(
            "'%s' is one of the reserved attributes %s, cannot delete. Use delncattr instead." % (name, tuple(_private_atts)))

    def delncattr(self, name):
        """
**`delncattr(self,name,value)`**

delete a netCDF dataset or group attribute.  Use if you need to delete a
netCDF attribute with the same name as one of the reserved python
attributes."""
        cdef char *attname
        cdef int ierr
        bytestr = _strencode(name)
        attname = bytestr
        if self.data_model != 'NETCDF4': self._redef()
        with nogil:
            ierr = nc_del_att(self._grpid, NC_GLOBAL, attname)
        if self.data_model != 'NETCDF4': self._enddef()
        _ensure_nc_success(ierr)

    def __setattr__(self,name,value):
        # if name in _private_atts, it is stored at the python
        # level and not in the netCDF file.
        if name not in _private_atts:
            self.setncattr(name, value)
        elif not name.endswith('__'):
            if hasattr(self,name):
                raise AttributeError(
            "'%s' is one of the reserved attributes %s, cannot rebind. Use setncattr instead." % (name, tuple(_private_atts)))
            else:
                self.__dict__[name]=value

    def __getattr__(self,name):
        # if name in _private_atts, it is stored at the python
        # level and not in the netCDF file.
        if name.startswith('__') and name.endswith('__'):
            # if __dict__ requested, return a dict with netCDF attributes.
            if name == '__dict__':
                names = self.ncattrs()
                values = []
                for name in names:
                    values.append(_get_att(self, NC_GLOBAL, name))
                return dict(zip(names, values))
            else:
                raise AttributeError
        elif name in _private_atts:
            return self.__dict__[name]
        else:
            return self.getncattr(name)

    def renameAttribute(self, oldname, newname):
        """
**`renameAttribute(self, oldname, newname)`**

rename a `Dataset` or `Group` attribute named `oldname` to `newname`."""
        cdef char *oldnamec
        cdef char *newnamec
        cdef int ierr
        bytestr = _strencode(oldname)
        oldnamec = bytestr
        bytestr = _strencode(newname)
        newnamec = bytestr
        with nogil:
            ierr = nc_rename_att(self._grpid, NC_GLOBAL, oldnamec, newnamec)
        _ensure_nc_success(ierr)

    def renameGroup(self, oldname, newname):
        """
**`renameGroup(self, oldname, newname)`**

rename a `Group` named `oldname` to `newname` (requires netcdf >= 4.3.1)."""
        cdef char *newnamec
        cdef int grpid
        IF HAS_RENAME_GRP:
            cdef int ierr
            bytestr = _strencode(newname)
            newnamec = bytestr
            try:
                grp = self.groups[oldname]
                grpid = grp._grpid
            except KeyError:
                raise KeyError('%s not a valid group name' % oldname)
            with nogil:
                ierr = nc_rename_grp(grpid, newnamec)
            _ensure_nc_success(ierr)
            # remove old key from groups dict.
            self.groups.pop(oldname)
            # add new key.
            self.groups[newname] = grp
        ELSE:
            msg = """
renameGroup method not enabled.  To enable, install Cython, make sure you have
version 4.3.1 or higher of the netcdf C lib, and rebuild netcdf4-python."""
            raise ValueError(msg)

    def set_auto_chartostring(self, value):
        """
**`set_auto_chartostring(self, True_or_False)`**

Call `Variable.set_auto_chartostring` for all variables contained in this `Dataset` or
`Group`, as well as for all variables in all its subgroups.

**`True_or_False`**: Boolean determining if automatic conversion of
all character arrays <--> string arrays should be performed for
character variables (variables of type `NC_CHAR` or `S1`) with the
`_Encoding` attribute set.

***Note***: Calling this function only affects existing variables. Variables created
after calling this function will follow the default behaviour.
        """

        # this is a hack to make inheritance work in MFDataset
        # (which stores variables in _vars)
        _vars = self.variables
        if _vars is None: _vars = self._vars
        for var in _vars.values():
            var.set_auto_chartostring(value)

        for groups in _walk_grps(self):
            for group in groups:
                for var in group.variables.values():
                    var.set_auto_chartostring(value)

    def set_auto_maskandscale(self, value):
        """
**`set_auto_maskandscale(self, True_or_False)`**

Call `Variable.set_auto_maskandscale` for all variables contained in this `Dataset` or
`Group`, as well as for all variables in all its subgroups.

**`True_or_False`**: Boolean determining if automatic conversion to masked arrays
and variable scaling shall be applied for all variables.

***Note***: Calling this function only affects existing variables. Variables created
after calling this function will follow the default behaviour.
        """

        # this is a hack to make inheritance work in MFDataset
        # (which stores variables in _vars)
        _vars = self.variables
        if _vars is None: _vars = self._vars
        for var in _vars.values():
            var.set_auto_maskandscale(value)

        for groups in _walk_grps(self):
            for group in groups:
                for var in group.variables.values():
                    var.set_auto_maskandscale(value)


    def set_auto_mask(self, value):
        """
**`set_auto_mask(self, True_or_False)`**

Call `Variable.set_auto_mask` for all variables contained in this `Dataset` or
`Group`, as well as for all variables in all its subgroups. Only affects
Variables with primitive or enum types (not compound or vlen Variables).

**`True_or_False`**: Boolean determining if automatic conversion to masked arrays
shall be applied for all variables.

***Note***: Calling this function only affects existing variables. Variables created
after calling this function will follow the default behaviour.
        """

        # this is a hack to make inheritance work in MFDataset
        # (which stores variables in _vars)
        _vars = self.variables
        if _vars is None: _vars = self._vars
        for var in _vars.values():
            var.set_auto_mask(value)

        for groups in _walk_grps(self):
            for group in groups:
                for var in group.variables.values():
                    var.set_auto_mask(value)

    def set_auto_scale(self, value):
        """
**`set_auto_scale(self, True_or_False)`**

Call `Variable.set_auto_scale` for all variables contained in this `Dataset` or
`Group`, as well as for all variables in all its subgroups.

**`True_or_False`**: Boolean determining if automatic variable scaling
shall be applied for all variables.

***Note***: Calling this function only affects existing variables. Variables created
after calling this function will follow the default behaviour.
        """

        # this is a hack to make inheritance work in MFDataset
        # (which stores variables in _vars)
        _vars = self.variables
        if _vars is None: _vars = self._vars
        for var in _vars.values():
            var.set_auto_scale(value)

        for groups in _walk_grps(self):
            for group in groups:
                for var in group.variables.values():
                    var.set_auto_scale(value)

    def set_always_mask(self, value):
        """
**`set_always_mask(self, True_or_False)`**

Call `Variable.set_always_mask` for all variables contained in
this `Dataset` or `Group`, as well as for all
variables in all its subgroups.

**`True_or_False`**: Boolean determining if automatic conversion of
masked arrays with no missing values to regular numpy arrays shall be
applied for all variables. Default True. Set to False to restore the default behaviour
in versions prior to 1.4.1 (numpy array returned unless missing values are present,
otherwise masked array returned).

***Note***: Calling this function only affects existing
variables. Variables created after calling this function will follow
the default behaviour.
        """

        # this is a hack to make inheritance work in MFDataset
        # (which stores variables in _vars)
        _vars = self.variables
        if _vars is None: _vars = self._vars
        for var in _vars.values():
            var.set_always_mask(value)

        for groups in _walk_grps(self):
            for group in groups:
                for var in group.variables.values():
                    var.set_always_mask(value)

    def set_ncstring_attrs(self, value):
        """
**`set_ncstring_attrs(self, True_or_False)`**

Call `Variable.set_ncstring_attrs` for all variables contained in
this `Dataset` or `Group`, as well as for all its
subgroups and their variables.

**`True_or_False`**: Boolean determining if all string attributes are
created as variable-length NC_STRINGs, (if True), or if ascii text
attributes are stored as NC_CHARs (if False; default)

***Note***: Calling this function only affects newly created attributes
of existing (sub-) groups and their variables.
        """

        self._ncstring_attrs__ = bool(value)

        # this is a hack to make inheritance work in MFDataset
        # (which stores variables in _vars)
        _vars = self.variables
        if _vars is None:
            _vars = self._vars
        for var in _vars.values():
            var.set_ncstring_attrs(value)

        for groups in _walk_grps(self):
            for group in groups:
                group.set_ncstring_attrs(value) # recurse into subgroups...

    def get_variables_by_attributes(self, **kwargs):
        """
**`get_variables_by_attribute(self, **kwargs)`**

Returns a list of variables that match specific conditions.

Can pass in key=value parameters and variables are returned that
contain all of the matches. For example,

```python
>>> # Get variables with x-axis attribute.
>>> vs = nc.get_variables_by_attributes(axis='X')
>>> # Get variables with matching "standard_name" attribute
>>> vs = nc.get_variables_by_attributes(standard_name='northward_sea_water_velocity')
```

Can pass in key=callable parameter and variables are returned if the
callable returns True.  The callable should accept a single parameter,
the attribute value.  None is given as the attribute value when the
attribute does not exist on the variable. For example,

```python
>>> # Get Axis variables
>>> vs = nc.get_variables_by_attributes(axis=lambda v: v in ['X', 'Y', 'Z', 'T'])
>>> # Get variables that don't have an "axis" attribute
>>> vs = nc.get_variables_by_attributes(axis=lambda v: v is None)
>>> # Get variables that have a "grid_mapping" attribute
>>> vs = nc.get_variables_by_attributes(grid_mapping=lambda v: v is not None)
```
"""
        vs = []

        has_value_flag  = False
        # this is a hack to make inheritance work in MFDataset
        # (which stores variables in _vars)
        _vars = self.variables
        if _vars is None: _vars = self._vars
        for vname in _vars:
            var = _vars[vname]
            for k, v in kwargs.items():
                if callable(v):
                    has_value_flag = v(getattr(var, k, None))
                    if has_value_flag is False:
                        break
                elif hasattr(var, k) and getattr(var, k) == v:
                    has_value_flag = True
                else:
                    has_value_flag = False
                    break

            if has_value_flag is True:
                vs.append(_vars[vname])

        return vs

    def _getname(self):
        # private method to get name associated with instance.
        cdef int ierr
        cdef char namstring[NC_MAX_NAME+1]
        with nogil:
            ierr = nc_inq_grpname(self._grpid, namstring)
        _ensure_nc_success(ierr)
        return namstring.decode('utf-8')

    property name:
        """string name of Group instance"""
        def __get__(self):
            return self._getname()
        def __set__(self,value):
            raise AttributeError("name cannot be altered")

    @staticmethod
    def fromcdl(cdlfilename,ncfilename=None,mode='a',format='NETCDF4'):
        """
**`fromcdl(cdlfilename, ncfilename=None, mode='a',format='NETCDF4')`**

call [ncgen][ncgen] via subprocess to create Dataset from [CDL][cdl]
text representation. Requires [ncgen][ncgen] to be installed and in `$PATH`.

**`cdlfilename`**:  CDL file.

**`ncfilename`**: netCDF file to create. If not given, CDL filename with
suffix replaced by `.nc` is used..

**`mode`**:  Access mode to open Dataset (Default `'a'`).

**`format`**: underlying file format to use (one of `'NETCDF4',
'NETCDF4_CLASSIC', 'NETCDF3_CLASSIC'`, `'NETCDF3_64BIT_OFFSET'` or
`'NETCDF3_64BIT_DATA'`. Default `'NETCDF4'`.
       
Dataset instance for `ncfilename` is returned.
 
[ncgen]: https://www.unidata.ucar.edu/software/netcdf/docs/netcdf_utilities_guide.html#ncgen_guide
[cdl]: https://www.unidata.ucar.edu/software/netcdf/docs/netcdf_utilities_guide.html#cdl_guide
        """
        if ncfilename is None:
            filepath = pathlib.Path(cdlfilename)
            ncfilename = filepath.with_suffix('.nc') 
        formatcodes = {'NETCDF4': 4,
                       'NETCDF4_CLASSIC': 7,
                       'NETCDF3_CLASSIC': 3,
                       'NETCDF3_64BIT': 6, # legacy
                       'NETCDF3_64BIT_OFFSET': 6,
                       'NETCDF3_64BIT_DATA': 5}
        if format not in formatcodes:
            raise ValueError('illegal format requested')
        ncgenargs="-knc%s" % formatcodes[format]
        subprocess.run(["ncgen", ncgenargs, "-o", ncfilename, cdlfilename], check=True)
        return Dataset(ncfilename, mode=mode)

    def tocdl(self,coordvars=False,data=False,outfile=None):
        """
**`tocdl(self, coordvars=False, data=False, outfile=None)`**

call [ncdump][ncdump] via subprocess to create [CDL][cdl]
text representation of Dataset. Requires [ncdump][ncdump]
to be installed and in `$PATH`.

**`coordvars`**: include coordinate variable data (via `ncdump -c`). Default False

**`data`**: if True, write out variable data (Default False).

**`outfile`**: If not None, file to output ncdump to. Default is to return a string.

[ncdump]: https://www.unidata.ucar.edu/software/netcdf/docs/netcdf_utilities_guide.html#ncdump_guide
[cdl]: https://www.unidata.ucar.edu/software/netcdf/docs/netcdf_utilities_guide.html#cdl_guide
        """
        self.sync()
        if coordvars:
            ncdumpargs = "-cs"
        else:
            ncdumpargs = "-s"
        if not data: ncdumpargs += "h"
        result=subprocess.run(["ncdump", ncdumpargs, self.filepath()],
               check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
               encoding='utf-8')
        if outfile is None:
            return result.stdout
        else:
            f = open(outfile,'w')
            f.write(result.stdout)
            f.close()
    def has_blosc_filter(self):
        """
**`has_blosc_filter(self)`**
returns True if blosc compression filter is available"""
        cdef int ierr
        IF HAS_BLOSC_SUPPORT:
            with nogil:
                ierr = nc_inq_filter_avail(self._grpid, H5Z_FILTER_BLOSC)
            if ierr:
                return False
            else:
                return True
        ELSE:
            return False
    def has_zstd_filter(self):
        """
**`has_zstd_filter(self)`**
returns True if zstd compression filter is available"""
        cdef int ierr
        IF HAS_ZSTANDARD_SUPPORT:
            with nogil:
                ierr = nc_inq_filter_avail(self._grpid, H5Z_FILTER_ZSTD)
            if ierr:
                return False
            else:
                return True
        ELSE:
            return False
    def has_bzip2_filter(self):
        """
**`has_bzip2_filter(self)`**
returns True if bzip2 compression filter is available"""
        cdef int ierr
        IF HAS_BZIP2_SUPPORT:
            with nogil:
                ierr = nc_inq_filter_avail(self._grpid, H5Z_FILTER_BZIP2)
            if ierr:
                return False
            else:
                return True
        ELSE:
            return False
    def has_szip_filter(self):
        """
**`has_szip_filter(self)`**
returns True if szip compression filter is available"""
        cdef int ierr
        IF HAS_NCFILTER:
            IF HAS_SZIP_SUPPORT:
                with nogil:
                    ierr = nc_inq_filter_avail(self._grpid, H5Z_FILTER_SZIP)
                if ierr:
                    return False
                else:
                    return True
            ELSE:
                return False
        ELSE:
             IF HAS_SZIP_SUPPORT:
                 return True
             ELSE:
                 return False

cdef class Group(Dataset):
    """
Groups define a hierarchical namespace within a netCDF file. They are
analogous to directories in a unix filesystem. Each `Group` behaves like
a `Dataset` within a Dataset, and can contain it's own variables,
dimensions and attributes (and other Groups). See `Group.__init__`
for more details.

`Group` inherits from `Dataset`, so all the
`Dataset` class methods and variables are available
to a `Group` instance (except the `close` method).

Additional read-only class variables:

**`name`**: String describing the group name.
    """
    def __init__(self, parent, name, **kwargs):
        """
        **`__init__(self, parent, name)`**
        `Group` constructor.

        **`parent`**: `Group` instance for the parent group.  If being created
        in the root group, use a `Dataset` instance.

        **`name`**: - Name of the group.

        ***Note***: `Group` instances should be created using the
        `Dataset.createGroup` method of a `Dataset` instance, or
        another `Group` instance, not using this class directly.
        """
        cdef char *groupname
        cdef int ierr, grpid
        # flag to indicate that Variables in this Group support orthogonal indexing.
        self.__orthogonal_indexing__ = True
        # set data_model and file_format attributes.
        self.data_model = parent.data_model
        self.file_format = parent.file_format
        # full path to Group.
        self.path = posixpath.join(parent.path, name)
        # parent group.
        self.parent = parent
        # propagate weak reference setting from parent.
        self.keepweakref = parent.keepweakref
        # propagate _ncstring_attrs__ setting from parent.
        self._ncstring_attrs__ = parent._ncstring_attrs__
        if 'id' in kwargs:
            self._grpid = kwargs['id']
            # get compound, vlen and enum types in this Group.
            self.cmptypes, self.vltypes, self.enumtypes = _get_types(self)
            # get dimensions in this Group.
            self.dimensions = _get_dims(self)
            # get variables in this Group.
            self.variables = _get_vars(self)
            # get groups in this Group.
            self.groups = _get_grps(self)
        else:
            bytestr = _strencode(name)
            groupname = bytestr
            grpid = parent._grpid
            with nogil:
                ierr = nc_def_grp(grpid, groupname, &self._grpid)
            _ensure_nc_success(ierr)
            self.cmptypes = dict()
            self.vltypes = dict()
            self.enumtypes = dict()
            self.dimensions = dict()
            self.variables = dict()
            self.groups = dict()


    def close(self):
        """
**`close(self)`**

overrides `Dataset` close method which does not apply to `Group`
instances, raises OSError."""
        raise OSError('cannot close a `Group` (only applies to Dataset)')


cdef class Dimension:
    """
A netCDF `Dimension` is used to describe the coordinates of a `Variable`.
See `Dimension.__init__` for more details.

The current maximum size of a `Dimension` instance can be obtained by
calling the python `len` function on the `Dimension` instance. The
`Dimension.isunlimited` method of a `Dimension` instance can be used to
determine if the dimension is unlimited.

Read-only class variables:

**`name`**: String name, used when creating a `Variable` with
`Dataset.createVariable`.

**`size`**: Current `Dimension` size (same as `len(d)`, where `d` is a
`Dimension` instance).
    """
    cdef public int _dimid, _grpid
    cdef public _data_model, _name, _grp

    def __init__(self, grp, name, size=None, **kwargs):
        """
        **`__init__(self, group, name, size=None)`**

        `Dimension` constructor.

        **`group`**: `Group` instance to associate with dimension.

        **`name`**: Name of the dimension.

        **`size`**: Size of the dimension. `None` or 0 means unlimited. (Default `None`).

        ***Note***: `Dimension` instances should be created using the
        `Dataset.createDimension` method of a `Group` or
        `Dataset` instance, not using `Dimension.__init__` directly.
        """
        cdef int ierr
        cdef char *dimname
        cdef size_t lendim
        self._grpid = grp._grpid
        # make a weakref to group to avoid circular ref (issue 218)
        # keep strong reference the default behaviour (issue 251)
        if grp.keepweakref:
            self._grp = weakref.proxy(grp)
        else:
            self._grp = grp
        self._data_model = grp.data_model
        self._name = name
        if 'id' in kwargs:
            self._dimid = kwargs['id']
        else:
            bytestr = _strencode(name)
            dimname = bytestr
            if size is not None:
                lendim = size
            else:
                lendim = NC_UNLIMITED
            if grp.data_model != 'NETCDF4': grp._redef()
            with nogil:
                ierr = nc_def_dim(self._grpid, dimname, lendim, &self._dimid)
            if grp.data_model != 'NETCDF4': grp._enddef()
            _ensure_nc_success(ierr)

    def _getname(self):
        # private method to get name associated with instance.
        cdef int err, _grpid
        cdef char namstring[NC_MAX_NAME+1]
        _grpid = self._grp._grpid
        with nogil:
            ierr = nc_inq_dimname(_grpid, self._dimid, namstring)
        _ensure_nc_success(ierr)
        return namstring.decode('utf-8')

    property name:
        """string name of Dimension instance"""
        def __get__(self):
            return self._getname()
        def __set__(self,value):
            raise AttributeError("name cannot be altered")

    property size:
        """current size of Dimension (calls `len` on Dimension instance)"""
        def __get__(self):
            return len(self)
        def __set__(self,value):
            raise AttributeError("size cannot be altered")

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        if not dir(self._grp):
            return 'Dimension object no longer valid'
        if self.isunlimited():
            return "%r (unlimited): name = '%s', size = %s" %\
                (type(self), self._name, len(self))
        else:
            return "%r: name = '%s', size = %s" %\
                (type(self), self._name, len(self))

    def __len__(self):
        # len(`Dimension` instance) returns current size of dimension
        cdef int ierr
        cdef size_t lengthp
        with nogil:
            ierr = nc_inq_dimlen(self._grpid, self._dimid, &lengthp)
        _ensure_nc_success(ierr)
        return lengthp

    def group(self):
        """
**`group(self)`**

return the group that this `Dimension` is a member of."""
        return self._grp

    def isunlimited(self):
        """
**`isunlimited(self)`**

returns `True` if the `Dimension` instance is unlimited, `False` otherwise."""
        cdef int ierr, n, numunlimdims, ndims, nvars, ngatts, xdimid
        cdef int *unlimdimids
        if self._data_model == 'NETCDF4':
            with nogil:
                ierr = nc_inq_unlimdims(self._grpid, &numunlimdims, NULL)
            _ensure_nc_success(ierr)
            if numunlimdims == 0:
                return False
            else:
                unlimdimids = <int *>malloc(sizeof(int) * numunlimdims)
                dimid = self._dimid
                with nogil:
                    ierr = nc_inq_unlimdims(self._grpid, &numunlimdims, unlimdimids)
                _ensure_nc_success(ierr)
                unlimdim_ids = []
                for n from 0 <= n < numunlimdims:
                    unlimdim_ids.append(unlimdimids[n])
                free(unlimdimids)
                if dimid in unlimdim_ids:
                    return True
                else:
                    return False
        else: # if not NETCDF4, there is only one unlimited dimension.
            # nc_inq_unlimdims only works for NETCDF4.
            with nogil:
                ierr = nc_inq(self._grpid, &ndims, &nvars, &ngatts, &xdimid)
            if self._dimid == xdimid:
                return True
            else:
                return False

cdef class Variable:
    """
A netCDF `Variable` is used to read and write netCDF data.  They are
analogous to numpy array objects. See `Variable.__init__` for more
details.

A list of attribute names corresponding to netCDF attributes defined for
the variable can be obtained with the `Variable.ncattrs` method. These
attributes can be created by assigning to an attribute of the
`Variable` instance. A dictionary containing all the netCDF attribute
name/value pairs is provided by the `__dict__` attribute of a
`Variable` instance.

The following class variables are read-only:

**`dimensions`**: A tuple containing the names of the
dimensions associated with this variable.

**`dtype`**: A numpy dtype object describing the
variable's data type.

**`ndim`**: The number of variable dimensions.

**`shape`**: A tuple with the current shape (length of all dimensions).

**`scale`**: If True, `scale_factor` and `add_offset` are
applied, and signed integer data is automatically converted to
unsigned integer data if the `_Unsigned` attribute is set to "true" or "True".
Default is `True`, can be reset using `Variable.set_auto_scale` and
`Variable.set_auto_maskandscale` methods.

**`mask`**: If True, data is automatically converted to/from masked
arrays when missing values or fill values are present. Default is `True`, can be
reset using `Variable.set_auto_mask` and `Variable.set_auto_maskandscale`
methods. Only relevant for Variables with primitive or enum types (ignored
for compound and vlen Variables).

**`chartostring`**: If True, data is automatically converted to/from character
arrays to string arrays when the `_Encoding` variable attribute is set.
Default is `True`, can be reset using
`Variable.set_auto_chartostring` method.

**`least_significant_digit`**: Describes the power of ten of the
smallest decimal place in the data the contains a reliable value.  Data is
truncated to this decimal place when it is assigned to the `Variable`
instance. If `None`, the data is not truncated.

**`significant_digits`**: New in version 1.6.0. Describes the number of significant
digits in the data the contains a reliable value.  Data is
truncated to retain this number of significant digits when it is assigned to the
`Variable` instance. If `None`, the data is not truncated.
Only available with netcdf-c >= 4.9.0,
and only works with `NETCDF4` or `NETCDF4_CLASSIC` formatted files.
The number of significant digits used in the quantization of variable data can be
obtained using the `Variable.significant_digits` method. Default `None` - 
no quantization done.

**`quantize_mode`**: New in version 1.6.0. Controls
the quantization algorithm (default 'BitGroom', 'BitRound' and
'GranularBitRound' also available).  The 'GranularBitRound'
algorithm may result in better compression for typical geophysical datasets. 
Ignored if `significant_digits` not specified. If 'BitRound' is used, then
`significant_digits` is interpreted as binary (not decimal) digits.

**`__orthogonal_indexing__`**: Always `True`.  Indicates to client code
that the object supports 'orthogonal indexing', which means that slices
that are 1d arrays or lists slice along each dimension independently.  This
behavior is similar to Fortran or Matlab, but different than numpy.

**`datatype`**: numpy data type (for primitive data types) or VLType/CompoundType
 instance (for compound or vlen data types).

**`name`**: String name.

**`size`**: The number of stored elements.
    """
    cdef public int _varid, _grpid, _nunlimdim
    cdef public _name, ndim, dtype, mask, scale, always_mask, chartostring,  _isprimitive, \
    _iscompound, _isvlen, _isenum, _grp, _cmptype, _vltype, _enumtype,\
    __orthogonal_indexing__, _has_lsd, _use_get_vars, _ncstring_attrs__

    def __init__(self, grp, name, datatype, dimensions=(),
            compression=None, zlib=False,
            complevel=4, shuffle=True, szip_coding='nn', szip_pixels_per_block=8,
            blosc_shuffle=1, 
            fletcher32=False, contiguous=False,
            chunksizes=None, endian='native', least_significant_digit=None,
            significant_digits=None,quantize_mode='BitGroom',fill_value=None, chunk_cache=None, **kwargs):
        """
        **`__init__(self, group, name, datatype, dimensions=(), compression=None, zlib=False,
        complevel=4, shuffle=True, szip_coding='nn', szip_pixels_per_block=8,
        blosc_shuffle=1, fletcher32=False, contiguous=False,
        chunksizes=None, endian='native',
        least_significant_digit=None,fill_value=None,chunk_cache=None)`**

        `Variable` constructor.

        **`group`**: `Group` or `Dataset` instance to associate with variable.

        **`name`**: Name of the variable.

        **`datatype`**: `Variable` data type. Can be specified by providing a
        numpy dtype object, or a string that describes a numpy dtype object.
        Supported values, corresponding to `str` attribute of numpy dtype
        objects, include `'f4'` (32-bit floating point), `'f8'` (64-bit floating
        point), `'i4'` (32-bit signed integer), `'i2'` (16-bit signed integer),
        `'i8'` (64-bit signed integer), `'i4'` (8-bit signed integer), `'i1'`
        (8-bit signed integer), `'u1'` (8-bit unsigned integer), `'u2'` (16-bit
        unsigned integer), `'u4'` (32-bit unsigned integer), `'u8'` (64-bit
        unsigned integer), or `'S1'` (single-character string).  From
        compatibility with Scientific.IO.NetCDF, the old Numeric single character
        typecodes can also be used (`'f'` instead of `'f4'`, `'d'` instead of
        `'f8'`, `'h'` or `'s'` instead of `'i2'`, `'b'` or `'B'` instead of
        `'i1'`, `'c'` instead of `'S1'`, and `'i'` or `'l'` instead of
        `'i4'`). `datatype` can also be a `CompoundType` instance
        (for a structured, or compound array), a `VLType` instance
        (for a variable-length array), or the python `str` builtin
        (for a variable-length string array). Numpy string and unicode datatypes with
        length greater than one are aliases for `str`.

        **`dimensions`**: a tuple containing the variable's Dimension instances
        (defined previously with `createDimension`). Default is an empty tuple
        which means the variable is a scalar (and therefore has no dimensions).

        **`compression`**: compression algorithm to use. 
        Currently `zlib`,`szip`,`zstd`,`bzip2`,`blosc_lz`,`blosc_lz4`,`blosc_lz4hc`,
        `blosc_zlib` and `blosc_zstd` are supported.
        Default is `None` (no compression).  All of the compressors except
        `zlib` and `szip` use the HDF5 plugin architecture.

        **`zlib`**: if `True`, data assigned to the `Variable`
        instance is compressed on disk. Default `False`. Deprecated - use
        `compression='zlib'` instead.

        **`complevel`**: the level of compression to use (1 is the fastest,
        but poorest compression, 9 is the slowest but best compression). Default 4.
        Ignored if `compression=None` or `szip`. A value of 0 disables compression.

        **`shuffle`**: if `True`, the HDF5 shuffle filter is applied
        to improve zlib compression. Default `True`. Ignored unless `compression = 'zlib'`.

        **`blosc_shuffle`**: shuffle filter inside blosc compressor (only
        relevant if compression kwarg set to one of the blosc compressors).
        Can be 0 (no blosc shuffle), 1 (bytewise shuffle) or 2 (bitwise
        shuffle)). Default is 1. Ignored if blosc compressor not used.

        **`szip_coding`**: szip coding method. Can be `ec` (entropy coding)
        or `nn` (nearest neighbor coding). Default is `nn`.
        Ignored if szip compressor not used.

        **`szip_pixels_per_block`**: Can be 4,8,16 or 32 (Default 8). 
        Ignored if szip compressor not used.

        **`fletcher32`**: if `True` (default `False`), the Fletcher32 checksum
        algorithm is used for error detection.

        **`contiguous`**: if `True` (default `False`), the variable data is
        stored contiguously on disk.  Default `False`. Setting to `True` for
        a variable with an unlimited dimension will trigger an error. Fixed
        size variables (with no unlimited dimension) with no compression
        filters are contiguous by default.

        **`chunksizes`**: Can be used to specify the HDF5 chunksizes for each
        dimension of the variable. A detailed discussion of HDF chunking and I/O
        performance is available
        [here](https://support.hdfgroup.org/HDF5/doc/Advanced/Chunking).
        The default chunking scheme in the netcdf-c library is discussed
        [here](https://www.unidata.ucar.edu/software/netcdf/documentation/NUG/netcdf_perf_chunking.html).
        Basically, you want the chunk size for each dimension to match as
        closely as possible the size of the data block that users will read
        from the file. `chunksizes` cannot be set if `contiguous=True`.

        **`endian`**: Can be used to control whether the
        data is stored in little or big endian format on disk. Possible
        values are `little, big` or `native` (default). The library
        will automatically handle endian conversions when the data is read,
        but if the data is always going to be read on a computer with the
        opposite format as the one used to create the file, there may be
        some performance advantage to be gained by setting the endian-ness.
        For netCDF 3 files (that don't use HDF5), only `endian='native'` is allowed.

        The `compression, zlib, complevel, shuffle, fletcher32, contiguous` and `chunksizes`
        keywords are silently ignored for netCDF 3 files that do not use HDF5.

        **`least_significant_digit`**: If this or `significant_digits` are specified,
        variable data will be truncated (quantized).
        In conjunction with `compression='zlib'` this produces
        'lossy', but significantly more efficient compression. For example, if
        `least_significant_digit=1`, data will be quantized using
        around(scale*data)/scale, where scale = 2**bits, and bits is determined
        so that a precision of 0.1 is retained (in this case bits=4). Default is
        `None`, or no quantization.

        **`significant_digits`**: New in version 1.6.0.
        As described for `least_significant_digit`
        except the number of significant digits retained is prescribed independent
        of the floating point exponent. Default `None` - no quantization done.

        **`quantize_mode`**: New in version 1.6.0. Controls
        the quantization algorithm (default 'BitGroom', 'BitRound' and
        'GranularBitRound' also available).  The 'GranularBitRound'
        algorithm may result in better compression for typical geophysical datasets.
        Ignored if `significant_digts` not specified. If 'BitRound' is used, then
        `significant_digits` is interpreted as binary (not decimal) digits.

        **`fill_value`**:  If specified, the default netCDF `_FillValue` (the
        value that the variable gets filled with before any data is written to it)
        is replaced with this value.  If fill_value is set to `False`, then
        the variable is not pre-filled. The default netCDF fill values can be found
        in the dictionary `netCDF4.default_fillvals`.

        **`chunk_cache`**: If specified, sets the chunk cache size for this variable.
        Persists as long as Dataset is open. Use `set_var_chunk_cache` to
        change it when Dataset is re-opened.

        ***Note***: `Variable` instances should be created using the
        `Dataset.createVariable` method of a `Dataset` or
        `Group` instance, not using this class directly.
        """
        cdef int ierr, ndims, icontiguous, icomplevel, numdims, _grpid, nsd,
        cdef unsigned int iblosc_complevel,iblosc_blocksize,iblosc_compressor,iblosc_shuffle
        cdef int iszip_coding, iszip_pixels_per_block
        cdef char namstring[NC_MAX_NAME+1]
        cdef char *varname
        cdef nc_type xtype
        cdef int *dimids
        cdef size_t sizep, nelemsp
        cdef size_t *chunksizesp
        cdef float preemptionp
        # flag to indicate that orthogonal indexing is supported
        self.__orthogonal_indexing__ = True
        # For backwards compatibility, deprecated zlib kwarg takes 
        # precedence if compression kwarg not set.
        if zlib and not compression:
            compression = 'zlib'
        # if complevel is set to zero, turn off compression
        if not complevel:
            compression = None
        zlib = False
        szip = False
        zstd = False
        bzip2 = False
        blosc_lz = False
        blosc_lz4 = False
        blosc_lz4hc = False
        #blosc_snappy = False
        blosc_zlib = False
        blosc_zstd = False
        if compression == 'zlib':
            zlib = True
        elif compression == 'szip':
            szip = True
        elif compression == 'zstd':
            zstd = True
        elif compression == 'bzip2':
            bzip2 = True
        elif compression == 'blosc_lz':
            blosc_lz = True
        elif compression == 'blosc_lz4':
            blosc_lz4 = True
        elif compression == 'blosc_lz4hc':
            blosc_lz4hc = True
        #elif compression == 'blosc_snappy':
        #    blosc_snappy = True
        elif compression == 'blosc_zlib':
            blosc_zlib = True
        elif compression == 'blosc_zstd':
            blosc_zstd = True
        elif not compression:
            compression = None # if compression evaluates to False, set to None.
            pass
        else:
            raise ValueError("Unsupported value for compression kwarg")
        self._grpid = grp._grpid
        # make a weakref to group to avoid circular ref (issue 218)
        # keep strong reference the default behaviour (issue 251)
        if grp.keepweakref:
            self._grp = weakref.proxy(grp)
        else:
            self._grp = grp
        user_type = isinstance(datatype, CompoundType) or \
                    isinstance(datatype, VLType) or \
                    isinstance(datatype, EnumType) or \
                    datatype == str
        # convert to a real numpy datatype object if necessary.
        if not user_type and type(datatype) != numpy.dtype:
            datatype = numpy.dtype(datatype)
        # convert numpy string dtype with length > 1
        # or any numpy unicode dtype into str
        if (isinstance(datatype, numpy.dtype) and
            ((datatype.kind == 'S' and datatype.itemsize > 1) or
              datatype.kind == 'U')):
            datatype = str
            user_type = True
        # check if endian keyword consistent with datatype specification.
        dtype_endian = getattr(datatype,'byteorder',None)
        if dtype_endian == '=': dtype_endian='native'
        if dtype_endian == '>': dtype_endian='big'
        if dtype_endian == '<': dtype_endian='little'
        if dtype_endian == '|': dtype_endian=None
        if dtype_endian is not None and dtype_endian != endian:
            if dtype_endian == 'native' and endian == sys.byteorder:
                pass
            else:
                # endian keyword prevails, issue warning
                msg = 'endian-ness of dtype and endian kwarg do not match, using endian kwarg'
                #msg = 'endian-ness of dtype and endian kwarg do not match, dtype over-riding endian kwarg'
                warnings.warn(msg)
                #endian = dtype_endian # dtype prevails
        # check validity of datatype.
        self._isprimitive = False
        self._iscompound = False
        self._isvlen = False
        self._isenum = False
        if user_type:
            if isinstance(datatype, CompoundType):
                self._iscompound = True
                self._cmptype = datatype
            if isinstance(datatype, VLType) or datatype==str:
                self._isvlen = True
                self._vltype = datatype
            if isinstance(datatype, EnumType):
                self._isenum = True
                self._enumtype = datatype
            if datatype==str:
                if grp.data_model != 'NETCDF4':
                    raise ValueError(
                        'Variable length strings are only supported for the '
                        'NETCDF4 format. For other formats, consider using '
                        'netCDF4.stringtochar to convert string arrays into '
                        'character arrays with an additional dimension.')
                datatype = VLType(self._grp, str, None)
                self._vltype = datatype
            xtype = datatype._nc_type
            # make sure this a valid user defined datatype defined in this Group
            with nogil:
                ierr = nc_inq_type(self._grpid, xtype, namstring, NULL)
            _ensure_nc_success(ierr)
            # dtype variable attribute is a numpy datatype object.
            self.dtype = datatype.dtype
        elif datatype.str[1:] in _supportedtypes:
            self._isprimitive = True
            # find netCDF primitive data type corresponding to
            # specified numpy data type.
            xtype = _nptonctype[datatype.str[1:]]
            # dtype variable attribute is a numpy datatype object.
            self.dtype = datatype
        else:
            raise TypeError('illegal primitive data type, must be one of %s, got %s' % (_supportedtypes,datatype))
        if 'id' in kwargs:
            self._varid = kwargs['id']
        else:
            bytestr = _strencode(name)
            varname = bytestr
            ndims = len(dimensions)
            # find dimension ids.
            if ndims:
                dimids = <int *>malloc(sizeof(int) * ndims)
                for n from 0 <= n < ndims:
                    dimids[n] = dimensions[n]._dimid
            # go into define mode if it's a netCDF 3 compatible
            # file format.  Be careful to exit define mode before
            # any exceptions are raised.
            if grp.data_model != 'NETCDF4': grp._redef()
            # define variable.
            if ndims:
                with nogil:
                    ierr = nc_def_var(self._grpid, varname, xtype, ndims,
                                      dimids, &self._varid)
                free(dimids)
            else: # a scalar variable.
                with nogil:
                    ierr = nc_def_var(self._grpid, varname, xtype, ndims,
                                      NULL, &self._varid)
            # set chunk cache size if desired
            # default is 1mb per var, can cause problems when many (1000's)
            # of vars are created.  This change only lasts as long as file is
            # open.
            if grp.data_model.startswith('NETCDF4') and chunk_cache is not None:
                with nogil:
                    ierr = nc_get_var_chunk_cache(self._grpid, self._varid, &sizep,
                           &nelemsp, &preemptionp)
                _ensure_nc_success(ierr)
                # reset chunk cache size, leave other parameters unchanged.
                sizep = chunk_cache
                with nogil:
                    ierr = nc_set_var_chunk_cache(self._grpid, self._varid, sizep,
                           nelemsp, preemptionp)
                _ensure_nc_success(ierr)
            if ierr != NC_NOERR:
                if grp.data_model != 'NETCDF4': grp._enddef()
                _ensure_nc_success(ierr)
            # set compression, shuffle, chunking, fletcher32 and endian
            # variable settings.
            # don't bother for NETCDF3* formats.
            # for NETCDF3* formats, the comopression,zlib,shuffle,chunking,
            # and fletcher32 flags are silently ignored. Only
            # endian='native' allowed for NETCDF3.
            if grp.data_model in ['NETCDF4','NETCDF4_CLASSIC']:
                # set compression and shuffle parameters.
                if compression is not None and ndims: # don't bother for scalar variable
                    if zlib:
                        icomplevel = complevel
                        if shuffle:
                            with nogil:
                                ierr = nc_def_var_deflate(self._grpid, self._varid, 1, 1, icomplevel)
                        else:
                            with nogil:
                                ierr = nc_def_var_deflate(self._grpid, self._varid, 0, 1, icomplevel)
                        if ierr != NC_NOERR:
                            if grp.data_model != 'NETCDF4': grp._enddef()
                            _ensure_nc_success(ierr)
                    if szip:
                        IF HAS_SZIP_SUPPORT:
                            try:
                                iszip_coding = _szip_dict[szip_coding]
                            except KeyError:
                                msg="unknown szip coding ('ec' or 'nn' supported)"
                                raise ValueError(msg)
                            iszip_pixels_per_block = szip_pixels_per_block
                            with nogil:
                                ierr = nc_def_var_szip(self._grpid, self._varid, iszip_coding, iszip_pixels_per_block)
                            if ierr != NC_NOERR:
                                if grp.data_model != 'NETCDF4': grp._enddef()
                                _ensure_nc_success(ierr)
                        ELSE:
                            msg = """
compression='szip' only works if linked version of hdf5 has szip functionality enabled"""
                            raise ValueError(msg)
                    if zstd:
                        IF HAS_ZSTANDARD_SUPPORT:
                            icomplevel = complevel
                            with nogil:
                                ierr = nc_def_var_zstandard(self._grpid, self._varid, icomplevel)
                            if ierr != NC_NOERR:
                                if grp.data_model != 'NETCDF4': grp._enddef()
                                _ensure_nc_success(ierr)
                        ELSE:
                            msg = """
compression='zstd' only works with netcdf-c >= 4.9.0.  To enable, install Cython, make sure you have
version 4.9.0 or higher netcdf-c with zstandard support, and rebuild netcdf4-python."""
                            raise ValueError(msg)
                    if bzip2:
                        IF HAS_BZIP2_SUPPORT:
                            icomplevel = complevel
                            with nogil:
                                ierr = nc_def_var_bzip2(self._grpid, self._varid, icomplevel)
                            if ierr != NC_NOERR:
                                if grp.data_model != 'NETCDF4': grp._enddef()
                                _ensure_nc_success(ierr)
                        ELSE:
                            msg = """
compression='bzip2' only works with netcdf-c >= 4.9.0.  To enable, install Cython, make sure you have
version 4.9.0 or higher netcdf-c with bzip2 support, and rebuild netcdf4-python."""
                            raise ValueError(msg)
                    if blosc_zstd or blosc_lz or blosc_lz4 or blosc_lz4hc or blosc_zlib:
                        IF HAS_BLOSC_SUPPORT:
                            iblosc_compressor = _blosc_dict[compression]
                            iblosc_shuffle = blosc_shuffle
                            iblosc_blocksize = 0 # not currently used by c lib
                            iblosc_complevel = complevel
                            with nogil:
                                ierr = nc_def_var_blosc(self._grpid, self._varid,\
                                    iblosc_compressor,\
                                    iblosc_complevel,iblosc_blocksize,\
                                    iblosc_shuffle)
                            if ierr != NC_NOERR:
                                if grp.data_model != 'NETCDF4': grp._enddef()
                                _ensure_nc_success(ierr)
                        ELSE:
                            msg = """
compression='blosc_*' only works with netcdf-c >= 4.9.0.  To enable, install Cython, make sure you have
version 4.9.0 or higher netcdf-c with blosc support, and rebuild netcdf4-python."""
                            raise ValueError(msg)
                # set checksum.
                if fletcher32 and ndims: # don't bother for scalar variable
                    with nogil:
                        ierr = nc_def_var_fletcher32(self._grpid, self._varid, 1)
                    if ierr != NC_NOERR:
                        if grp.data_model != 'NETCDF4': grp._enddef()
                        _ensure_nc_success(ierr)
                # set chunking stuff.
                if ndims: # don't bother for scalar variable.
                    if contiguous:
                        icontiguous = NC_CONTIGUOUS
                        if chunksizes is not None:
                            raise ValueError('cannot specify chunksizes for a contiguous dataset')
                    else:
                        icontiguous = NC_CHUNKED
                    if chunksizes is None:
                        chunksizesp = NULL
                    else:
                        if len(chunksizes) != len(dimensions):
                            if grp.data_model != 'NETCDF4': grp._enddef()
                            raise ValueError('chunksizes must be a sequence with the same length as dimensions')
                        chunksizesp = <size_t *>malloc(sizeof(size_t) * ndims)
                        for n from 0 <= n < ndims:
                            if not dimensions[n].isunlimited() and \
                               chunksizes[n] > dimensions[n].size:
                                msg = 'chunksize cannot exceed dimension size'
                                raise ValueError(msg)
                            chunksizesp[n] = chunksizes[n]
                    if chunksizes is not None or contiguous:
                        with nogil:
                            ierr = nc_def_var_chunking(self._grpid, self._varid, icontiguous, chunksizesp)
                        free(chunksizesp)
                        if ierr != NC_NOERR:
                            if grp.data_model != 'NETCDF4': grp._enddef()
                            _ensure_nc_success(ierr)
                # set endian-ness of variable
                if endian == 'little':
                    with nogil:
                        ierr = nc_def_var_endian(self._grpid, self._varid, NC_ENDIAN_LITTLE)
                elif endian == 'big':
                    with nogil:
                        ierr = nc_def_var_endian(self._grpid, self._varid, NC_ENDIAN_BIG)
                elif endian == 'native':
                    pass # this is the default format.
                else:
                    raise ValueError("'endian' keyword argument must be 'little','big' or 'native', got '%s'" % endian)
                # set quantization
                IF HAS_QUANTIZATION_SUPPORT:
                    if significant_digits is not None:
                        nsd = significant_digits
                        if quantize_mode == 'BitGroom':
                            with nogil:
                                ierr = nc_def_var_quantize(self._grpid,
                                       self._varid, NC_QUANTIZE_BITGROOM, nsd)
                        elif quantize_mode == 'GranularBitRound':
                            with nogil:
                                ierr = nc_def_var_quantize(self._grpid,
                                       self._varid, NC_QUANTIZE_GRANULARBR, nsd)
                        elif quantize_mode == 'BitRound':
                            ierr = nc_def_var_quantize(self._grpid,
                                       self._varid, NC_QUANTIZE_BITROUND, nsd)
                        else:
                            raise ValueError("'quantize_mode' keyword argument must be 'BitGroom','GranularBitRound' or 'BitRound', got '%s'" % quantize_mode)

                ELSE:
                    if significant_digits is not None:
                        msg = """
significant_digits kwarg only works with netcdf-c >= 4.9.0.  To enable, install Cython, make sure you have
version 4.9.0 or higher netcdf-c, and rebuild netcdf4-python. Otherwise, use least_significant_digit
kwarg for quantization."""
                        raise ValueError(msg)
                if ierr != NC_NOERR:
                    if grp.data_model != 'NETCDF4': grp._enddef()
                    _ensure_nc_success(ierr)
            else:
                if endian != 'native':
                    msg="only endian='native' allowed for NETCDF3 files"
                    raise RuntimeError(msg)
            # set a fill value for this variable if fill_value keyword
            # given.  This avoids the HDF5 overhead of deleting and
            # recreating the dataset if it is set later (after the enddef).
            if fill_value is not None:
                if not fill_value and isinstance(fill_value,bool):
                    # no filling for this variable if fill_value==False.
                    if not self._isprimitive:
                        # no fill values for VLEN and compound variables
                        # anyway.
                        ierr = 0
                    else:
                        with nogil:
                            ierr = nc_def_var_fill(self._grpid, self._varid, 1, NULL)
                    if ierr != NC_NOERR:
                        if grp.data_model != 'NETCDF4': grp._enddef()
                        _ensure_nc_success(ierr)
                else:
                    if self._isprimitive or self._isenum or \
                       (self._isvlen and self.dtype == str):
                        if self._isvlen and self.dtype == str:
                            _set_att(self._grp, self._varid, '_FillValue',\
                               _tostr(fill_value), xtype=xtype, force_ncstring=True)
                        else:
                            fillval = numpy.array(fill_value, self.dtype)
                            if not fillval.dtype.isnative: fillval.byteswap(True)
                            _set_att(self._grp, self._varid, '_FillValue',\
                                     fillval, xtype=xtype)
                    else:
                        raise AttributeError("cannot set _FillValue attribute for VLEN or compound variable")
            if least_significant_digit is not None:
                self.least_significant_digit = least_significant_digit
            # leave define mode if not a NETCDF4 format file.
            if grp.data_model != 'NETCDF4': grp._enddef()
        # count how many unlimited dimensions there are.
        self._nunlimdim = 0
        for dim in dimensions:
            if dim.isunlimited(): self._nunlimdim = self._nunlimdim + 1
        # set ndim attribute (number of dimensions).
        with nogil:
            ierr = nc_inq_varndims(self._grpid, self._varid, &numdims)
        _ensure_nc_success(ierr)
        self.ndim = numdims
        self._name = name
        # default for automatically applying scale_factor and
        # add_offset, and converting to/from masked arrays is True.
        self.scale = True
        self.mask = True
        # issue 809: default for converting arrays with no missing values to
        # regular numpy arrays
        self.always_mask = True
        # default is to automatically convert to/from character
        # to string arrays when _Encoding variable attribute is set.
        self.chartostring = True
        # propagate _ncstring_attrs__ setting from parent group.
        self._ncstring_attrs__ = grp._ncstring_attrs__
        if 'least_significant_digit' in self.ncattrs():
            self._has_lsd = True
        # avoid calling nc_get_vars for strided slices by default.
        # a fix for strided slice access using HDF5 was added
        # in 4.6.2.
        # always use nc_get_vars for strided access with OpenDAP (issue #838).
        if __netcdf4libversion__ >= "4.6.2" or\
           self._grp.filepath().startswith('http'):
            self._use_get_vars = True
        else:
            self._use_get_vars = False

    def __array__(self):
        # numpy special method that returns a numpy array.
        # allows numpy ufuncs to work faster on Variable objects
        # (issue 216).
        return self[...]

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        cdef int ierr, no_fill
        if not dir(self._grp):
            return 'Variable object no longer valid'
        ncdump = [repr(type(self))]
        show_more_dtype = True
        if self._iscompound:
            kind = 'compound'
        elif self._isvlen:
            kind = 'vlen'
        elif self._isenum:
            kind = 'enum'
        else:
            show_more_dtype = False
            kind = str(self.dtype)
        dimnames = tuple(_tostr(dimname) for dimname in self.dimensions)
        ncdump.append('%s %s(%s)' %\
            (kind, self._name, ', '.join(dimnames)))
        for name in self.ncattrs():
            ncdump.append('    %s: %s' % (name, self.getncattr(name)))
        if show_more_dtype:
            ncdump.append('%s data type: %s' % (kind, self.dtype))
        unlimdims = []
        for dimname in self.dimensions:
            dim = _find_dim(self._grp, dimname)
            if dim.isunlimited():
                unlimdims.append(dimname)
        if (self._grp.path != '/'): ncdump.append('path = %s' % self._grp.path)
        ncdump.append('unlimited dimensions: %s' % ', '.join(unlimdims))
        ncdump.append('current shape = %r' % (self.shape,))
        if __netcdf4libversion__ < '4.5.1' and\
            self._grp.file_format.startswith('NETCDF3'):
            # issue #908: no_fill not correct for NETCDF3 files before 4.5.1
            # before 4.5.1 there was no way to turn off filling on a
            # per-variable basis for classic files.
            no_fill=0
        else:
            with nogil:
                ierr = nc_inq_var_fill(self._grpid,self._varid,&no_fill,NULL)
            _ensure_nc_success(ierr)
        if self._isprimitive:
            if no_fill != 1:
                try:
                    fillval = self._FillValue
                    msg = 'filling on'
                except AttributeError:
                    fillval = default_fillvals[self.dtype.str[1:]]
                    if self.dtype.str[1:] in ['u1','i1']:
                        msg = 'filling on, default _FillValue of %s ignored' % fillval
                    else:
                        msg = 'filling on, default _FillValue of %s used' % fillval
                ncdump.append(msg)
            else:
                ncdump.append('filling off')


        return '\n'.join(ncdump)

    def _getdims(self):
        # Private method to get variables's dimension names
        cdef int ierr, numdims, n, nn
        cdef char namstring[NC_MAX_NAME+1]
        cdef int *dimids
        # get number of dimensions for this variable.
        with nogil:
            ierr = nc_inq_varndims(self._grpid, self._varid, &numdims)
        _ensure_nc_success(ierr)
        dimids = <int *>malloc(sizeof(int) * numdims)
        # get dimension ids.
        with nogil:
            ierr = nc_inq_vardimid(self._grpid, self._varid, dimids)
        _ensure_nc_success(ierr)
        # loop over dimensions, retrieve names.
        dimensions = ()
        for nn from 0 <= nn < numdims:
            with nogil:
                ierr = nc_inq_dimname(self._grpid, dimids[nn], namstring)
            _ensure_nc_success(ierr)
            name = namstring.decode('utf-8')
            dimensions = dimensions + (name,)
        free(dimids)
        return dimensions

    def _getname(self):
        # Private method to get name associated with instance
        cdef int err, _grpid
        cdef char namstring[NC_MAX_NAME+1]
        _grpid = self._grp._grpid
        with nogil:
            ierr = nc_inq_varname(_grpid, self._varid, namstring)
        _ensure_nc_success(ierr)
        return namstring.decode('utf-8')

    property name:
        """string name of Variable instance"""
        def __get__(self):
            return self._getname()
        def __set__(self,value):
            raise AttributeError("name cannot be altered")

    property datatype:
        """numpy data type (for primitive data types) or
        VLType/CompoundType/EnumType instance 
        (for compound, vlen  or enum data types)"""
        def __get__(self):
            if self._iscompound:
                return self._cmptype
            elif self._isvlen:
                return self._vltype
            elif self._isenum:
                return self._enumtype
            elif self._isprimitive:
                return self.dtype

    property shape:
        """find current sizes of all variable dimensions"""
        def __get__(self):
            shape = ()
            for dimname in self._getdims():
                # look in current group, and parents for dim.
                dim = _find_dim(self._grp,dimname)
                shape = shape + (len(dim),)
            return shape
        def __set__(self,value):
            raise AttributeError("shape cannot be altered")

    property size:
        """Return the number of stored elements."""
        def __get__(self):
            # issue #957: add int since prod(())=1.0
            return int(numpy.prod(self.shape))

    property dimensions:
        """get variables's dimension names"""
        def __get__(self):
            return self._getdims()
        def __set__(self,value):
            raise AttributeError("dimensions cannot be altered")


    def group(self):
        """
**`group(self)`**

return the group that this `Variable` is a member of."""
        return self._grp

    def ncattrs(self):
        """
**`ncattrs(self)`**

return netCDF attribute names for this `Variable` in a list."""
        return _get_att_names(self._grpid, self._varid)

    def setncattr(self,name,value):
        """
**`setncattr(self,name,value)`**

set a netCDF variable attribute using name,value pair.  Use if you need to set a
netCDF attribute with the same name as one of the reserved python
attributes."""
        cdef nc_type xtype
        xtype=-99
        # issue #959 - trying to set _FillValue results in mysterious
        # error when close method is called so catch it here. It is
        # already caught in __setattr__.
        if name == '_FillValue':
            msg='_FillValue attribute must be set when variable is '+\
            'created (using fill_value keyword to createVariable)'
            raise AttributeError(msg)
        if self._grp.data_model != 'NETCDF4': self._grp._redef()
        _set_att(self._grp, self._varid, name, value, xtype=xtype, force_ncstring=self._ncstring_attrs__)
        if self._grp.data_model != 'NETCDF4': self._grp._enddef()

    def setncattr_string(self,name,value):
        """
**`setncattr_string(self,name,value)`**

set a netCDF variable string attribute using name,value pair.
Use if you need to ensure that a netCDF attribute is created with type
`NC_STRING` if the file format is `NETCDF4`.
Use if you need to set an attribute to an array of variable-length strings."""
        cdef nc_type xtype
        xtype=-99
        if self._grp.data_model != 'NETCDF4':
            msg='file format does not support NC_STRING attributes'
            raise OSError(msg)
        _set_att(self._grp, self._varid, name, value, xtype=xtype, force_ncstring=True)

    def setncatts(self,attdict):
        """
**`setncatts(self,attdict)`**

set a bunch of netCDF variable attributes at once using a python dictionary.
This may be faster when setting a lot of attributes for a `NETCDF3`
formatted file, since nc_redef/nc_enddef is not called in between setting
each attribute"""
        if self._grp.data_model != 'NETCDF4': self._grp._redef()
        for name, value in attdict.items():
            _set_att(self._grp, self._varid, name, value)
        if self._grp.data_model != 'NETCDF4': self._grp._enddef()

    def getncattr(self,name,encoding='utf-8'):
        """
**`getncattr(self,name)`**

retrieve a netCDF variable attribute.  Use if you need to set a
netCDF attribute with the same name as one of the reserved python
attributes.

option kwarg `encoding` can be used to specify the
character encoding of a string attribute (default is `utf-8`)."""
        return _get_att(self._grp, self._varid, name, encoding=encoding)

    def delncattr(self, name):
        """
**`delncattr(self,name,value)`**

delete a netCDF variable attribute.  Use if you need to delete a
netCDF attribute with the same name as one of the reserved python
attributes."""
        cdef char *attname
        bytestr = _strencode(name)
        attname = bytestr
        if self._grp.data_model != 'NETCDF4': self._grp._redef()
        with nogil:
            ierr = nc_del_att(self._grpid, self._varid, attname)
        if self._grp.data_model != 'NETCDF4': self._grp._enddef()
        _ensure_nc_success(ierr)

    def filters(self):
        """
**`filters(self)`**

return dictionary containing HDF5 filter parameters."""
        cdef int ierr,ideflate,ishuffle,icomplevel,ifletcher32
        cdef int izstd=0
        cdef int ibzip2=0
        cdef int iblosc=0
        cdef int iszip=0
        cdef int iszip_coding=0
        cdef int iszip_pixels_per_block=0
        cdef int icomplevel_zstd=0 
        cdef int icomplevel_bzip2=0 
        cdef unsigned int iblosc_shuffle=0
        cdef unsigned int iblosc_compressor=0
        cdef unsigned int iblosc_blocksize=0
        cdef unsigned int iblosc_complevel=0
        filtdict = {'zlib':False,'szip':False,'zstd':False,'bzip2':False,'blosc':False,'shuffle':False,'complevel':0,'fletcher32':False}
        if self._grp.data_model not in ['NETCDF4_CLASSIC','NETCDF4']: return
        with nogil:
            ierr = nc_inq_var_deflate(self._grpid, self._varid, &ishuffle, &ideflate, &icomplevel)
        _ensure_nc_success(ierr)
        with nogil:
            ierr = nc_inq_var_fletcher32(self._grpid, self._varid, &ifletcher32)
        _ensure_nc_success(ierr)
        IF HAS_ZSTANDARD_SUPPORT:
            with nogil:
                ierr = nc_inq_var_zstandard(self._grpid, self._varid, &izstd,\
                       &icomplevel_zstd)
            if ierr != 0: izstd=0
            # _ensure_nc_success(ierr)
        IF HAS_BZIP2_SUPPORT:
            with nogil:
                ierr = nc_inq_var_bzip2(self._grpid, self._varid, &ibzip2,\
                       &icomplevel_bzip2)
            if ierr != 0: ibzip2=0
            #_ensure_nc_success(ierr)
        IF HAS_BLOSC_SUPPORT:
            with nogil:
                ierr = nc_inq_var_blosc(self._grpid, self._varid, &iblosc,\
                       &iblosc_compressor,&iblosc_complevel,&iblosc_blocksize,&iblosc_shuffle)
            if ierr != 0: iblosc=0
            #_ensure_nc_success(ierr)
        IF HAS_SZIP_SUPPORT:
            with nogil:
                ierr = nc_inq_var_szip(self._grpid, self._varid, &iszip_coding,\
                       &iszip_pixels_per_block)
            if ierr != 0:
                iszip=0
            else:
                if iszip_coding:
                    iszip=1
                else:
                    iszip=0
            #_ensure_nc_success(ierr)
        if ideflate:
            filtdict['zlib']=True
            filtdict['complevel']=icomplevel
        if izstd:
            filtdict['zstd']=True
            filtdict['complevel']=icomplevel_zstd
        if ibzip2:
            filtdict['bzip2']=True
            filtdict['complevel']=icomplevel_bzip2
        if iblosc:
            blosc_compressor = iblosc_compressor
            filtdict['blosc']={'compressor':_blosc_dict_inv[blosc_compressor],'shuffle':iblosc_shuffle}
            filtdict['complevel']=iblosc_complevel
        if iszip:
            szip_coding = iszip_coding
            filtdict['szip']={'coding':_szip_dict_inv[szip_coding],'pixels_per_block':iszip_pixels_per_block}
        if ishuffle:
            filtdict['shuffle']=True
        if ifletcher32:
            filtdict['fletcher32']=True
        return filtdict

    def quantization(self):
        """
**`quantization(self)`**

return number of significant digits and the algorithm used in quantization.
Returns None if quantization not active.
"""
        IF HAS_QUANTIZATION_SUPPORT:
            cdef int ierr, nsd, quantize_mode
            if self._grp.data_model not in ['NETCDF4_CLASSIC','NETCDF4']:
                return None
            else:
                with nogil:
                    ierr = nc_inq_var_quantize(self._grpid, self._varid, &quantize_mode, &nsd)
                _ensure_nc_success(ierr)
                if quantize_mode == NC_NOQUANTIZE:
                    return None
                else:
                    if quantize_mode == NC_QUANTIZE_GRANULARBR:
                        sig_digits = nsd
                        quant_mode = 'GranularBitRound'
                    elif quantize_mode == NC_QUANTIZE_BITROUND:
                        sig_digits = nsd # interpreted as bits, not decimal
                        quant_mode = 'BitRound'
                    else:
                        sig_digits = nsd
                        quant_mode = 'BitGroom'
                    return sig_digits, quant_mode
        ELSE:
            return None

    def endian(self):
        """
**`endian(self)`**

return endian-ness (`little,big,native`) of variable (as stored in HDF5 file)."""
        cdef int ierr, iendian
        if self._grp.data_model not in ['NETCDF4_CLASSIC','NETCDF4']:
            return 'native'
        with nogil:
            ierr = nc_inq_var_endian(self._grpid, self._varid, &iendian)
        _ensure_nc_success(ierr)
        if iendian == NC_ENDIAN_LITTLE:
            return 'little'
        elif iendian == NC_ENDIAN_BIG:
            return 'big'
        else:
            return 'native'

    def chunking(self):
        """
**`chunking(self)`**

return variable chunking information.  If the dataset is
defined to be contiguous (and hence there is no chunking) the word 'contiguous'
is returned.  Otherwise, a sequence with the chunksize for
each dimension is returned."""
        cdef int ierr, icontiguous, ndims
        cdef size_t *chunksizesp
        if self._grp.data_model not in ['NETCDF4_CLASSIC','NETCDF4']: return None
        ndims = self.ndim
        chunksizesp = <size_t *>malloc(sizeof(size_t) * ndims)
        with nogil:
            ierr = nc_inq_var_chunking(self._grpid, self._varid, &icontiguous, chunksizesp)
        _ensure_nc_success(ierr)
        chunksizes=[]
        for n from 0 <= n < ndims:
            chunksizes.append(chunksizesp[n])
        free(chunksizesp)
        if icontiguous:
            return 'contiguous'
        else:
            return chunksizes

    def get_var_chunk_cache(self):
        """
**`get_var_chunk_cache(self)`**

return variable chunk cache information in a tuple (size,nelems,preemption).
See netcdf C library documentation for `nc_get_var_chunk_cache` for
details."""
        cdef int ierr
        cdef size_t sizep, nelemsp
        cdef float preemptionp
        with nogil:
            ierr = nc_get_var_chunk_cache(self._grpid, self._varid, &sizep,
                   &nelemsp, &preemptionp)
        _ensure_nc_success(ierr)
        size = sizep; nelems = nelemsp; preemption = preemptionp
        return (size,nelems,preemption)

    def set_var_chunk_cache(self,size=None,nelems=None,preemption=None):
        """
**`set_var_chunk_cache(self,size=None,nelems=None,preemption=None)`**

change variable chunk cache settings.
See netcdf C library documentation for `nc_set_var_chunk_cache` for
details."""
        cdef int ierr
        cdef size_t sizep, nelemsp
        cdef float preemptionp
        # reset chunk cache size, leave other parameters unchanged.
        size_orig, nelems_orig, preemption_orig = self.get_var_chunk_cache()
        if size is not None:
            sizep = size
        else:
            sizep = size_orig
        if nelems is not None:
            nelemsp = nelems
        else:
            nelemsp = nelems_orig
        if preemption is not None:
            preemptionp = preemption
        else:
            preemptionp = preemption_orig
        with nogil:
            ierr = nc_set_var_chunk_cache(self._grpid, self._varid, sizep,
                   nelemsp, preemptionp)
        _ensure_nc_success(ierr)

    def __delattr__(self,name):
        # if it's a netCDF attribute, remove it
        if name not in _private_atts:
            self.delncattr(name)
        else:
            raise AttributeError(
            "'%s' is one of the reserved attributes %s, cannot delete. Use delncattr instead." % (name, tuple(_private_atts)))

    def __setattr__(self,name,value):
        # if name in _private_atts, it is stored at the python
        # level and not in the netCDF file.
        if name not in _private_atts:
            # if setting _FillValue or missing_value, make sure value
            # has same type and byte order as variable.
            if name == '_FillValue':
                msg='_FillValue attribute must be set when variable is '+\
                'created (using fill_value keyword to createVariable)'
                raise AttributeError(msg)
                #if self._isprimitive:
                #    value = numpy.array(value, self.dtype)
                #else:
                #    msg="cannot set _FillValue attribute for "+\
                #    "VLEN or compound variable"
                #    raise AttributeError(msg)
            elif name in ['valid_min','valid_max','valid_range','missing_value'] and self._isprimitive:
                # make sure these attributes written in same data type as variable.
                # also make sure it is written in native byte order
                # (the same as the data)
                valuea = numpy.array(value, self.dtype)
                # check to see if array cast is safe
                if _safecast(numpy.array(value),valuea):
                    value = valuea
                    if not value.dtype.isnative: value.byteswap(True)
                else: # otherwise don't do it, but issue a warning
                    msg="WARNING: %s cannot be safely cast to variable dtype" \
                    % name
                    warnings.warn(msg)
            self.setncattr(name, value)
        elif not name.endswith('__'):
            if hasattr(self,name):
                raise AttributeError(
                "'%s' is one of the reserved attributes %s, cannot rebind. Use setncattr instead." % (name, tuple(_private_atts)))
            else:
                self.__dict__[name]=value

    def __getattr__(self,name):
        # if name in _private_atts, it is stored at the python
        # level and not in the netCDF file.
        if name.startswith('__') and name.endswith('__'):
            # if __dict__ requested, return a dict with netCDF attributes.
            if name == '__dict__':
                names = self.ncattrs()
                values = []
                for name in names:
                    values.append(_get_att(self._grp, self._varid, name))
                return dict(zip(names, values))

            else:
                raise AttributeError
        elif name in _private_atts:
            return self.__dict__[name]
        else:
            return self.getncattr(name)

    def renameAttribute(self, oldname, newname):
        """
**`renameAttribute(self, oldname, newname)`**

rename a `Variable` attribute named `oldname` to `newname`."""
        cdef int ierr
        cdef char *oldnamec
        cdef char *newnamec
        bytestr = _strencode(oldname)
        oldnamec = bytestr
        bytestr = _strencode(newname)
        newnamec = bytestr
        with nogil:
            ierr = nc_rename_att(self._grpid, self._varid, oldnamec, newnamec)
        _ensure_nc_success(ierr)

    def __getitem__(self, elem):
        # This special method is used to index the netCDF variable
        # using the "extended slice syntax". The extended slice syntax
        # is a perfect match for the "start", "count" and "stride"
        # arguments to the nc_get_var() function, and is much more easy
        # to use.
        start, count, stride, put_ind =\
        _StartCountStride(elem,self.shape,dimensions=self.dimensions,grp=self._grp,use_get_vars=self._use_get_vars)
        datashape = _out_array_shape(count)
        if self._isvlen:
            data = numpy.empty(datashape, dtype='O')
        else:
            data = numpy.empty(datashape, dtype=self.dtype)

        # Determine which dimensions need to be
        # squeezed (those for which elem is an integer scalar).
        # The convention used is that for those cases,
        # put_ind for this dimension is set to -1 by _StartCountStride.
        squeeze = data.ndim * [slice(None),]
        for i,n in enumerate(put_ind.shape[:-1]):
            if n == 1 and put_ind.size > 0 and put_ind[...,i].ravel()[0] == -1:
                squeeze[i] = 0

        # Reshape the arrays so we can iterate over them.
        start = start.reshape((-1, self.ndim or 1))
        count = count.reshape((-1, self.ndim or 1))
        stride = stride.reshape((-1, self.ndim or 1))
        put_ind = put_ind.reshape((-1, self.ndim or 1))

        # Fill output array with data chunks.
        for (a,b,c,i) in zip(start, count, stride, put_ind):
            datout = self._get(a,b,c)
            if not hasattr(datout,'shape') or data.shape == datout.shape:
                data = datout
            else:
                shape = getattr(data[tuple(i)], 'shape', ())
                if self._isvlen and not len(self.dimensions):
                    # special case of scalar VLEN
                    data[0] = datout
                else:
                    data[tuple(i)] = datout.reshape(shape)

        # Remove extra singleton dimensions.
        if hasattr(data,'shape'):
            data = data[tuple(squeeze)]
        if hasattr(data,'ndim') and self.ndim == 0:
            # Make sure a numpy scalar array is returned instead of a 1-d array of
            # length 1.
            if data.ndim != 0: data = numpy.asarray(data[0])

        # if auto_scale mode set to True, (through
        # a call to set_auto_scale or set_auto_maskandscale),
        # perform automatic unpacking using scale_factor/add_offset.
        # if auto_mask mode is set to True (through a call to
        # set_auto_mask or set_auto_maskandscale), perform
        # automatic conversion to masked array using
        # missing_value/_Fill_Value.
        # applied for primitive and (non-string) vlen,
        # ignored for compound and enum datatypes.
        try: # check to see if scale_factor and add_offset is valid (issue 176).
            if hasattr(self,'scale_factor'): float(self.scale_factor)
            if hasattr(self,'add_offset'): float(self.add_offset)
            valid_scaleoffset = True
        except:
            valid_scaleoffset = False
            if self.scale:
                msg = 'invalid scale_factor or add_offset attribute, no unpacking done...'
                warnings.warn(msg)

        if self.mask and (self._isprimitive or self._isenum):\
            data = self._toma(data)
        else:
            # if attribute _Unsigned is "true", and variable has signed integer
            # dtype, return view with corresponding unsigned dtype (issue #656)
            if self.scale:  # only do this if autoscale option is on.
                is_unsigned = getattr(self, '_Unsigned', False) in ["true","True"]
                if is_unsigned and data.dtype.kind == 'i':
                    data=data.view('%su%s'%(data.dtype.byteorder,data.dtype.itemsize))

        if self.scale and\
           (self._isprimitive or (self._isvlen and self.dtype != str)) and\
           valid_scaleoffset:
            # if variable has scale_factor and add_offset attributes, apply
            # them.
            if hasattr(self, 'scale_factor') and hasattr(self, 'add_offset'):
                if self.add_offset != 0.0 or self.scale_factor != 1.0:
                    data = data*self.scale_factor + self.add_offset
                else:
                    data = data.astype(self.scale_factor.dtype) # issue 913
            # else if variable has only scale_factor attribute, rescale.
            elif hasattr(self, 'scale_factor') and self.scale_factor != 1.0:
                data = data*self.scale_factor
            # else if variable has only add_offset attribute, add offset.
            elif hasattr(self, 'add_offset') and self.add_offset != 0.0:
                data = data + self.add_offset

        # if _Encoding is specified for a character variable, return
        # a numpy array of strings with one less dimension.
        if self.chartostring and getattr(self.dtype,'kind',None) == 'S' and\
           getattr(self.dtype,'itemsize',None) == 1:
            encoding = getattr(self,'_Encoding',None)
            # should this only be done if self.scale = True?
            # should there be some other way to disable this?
            if encoding is not None:
                # only try to return a string array if rightmost dimension of
                # sliced data matches rightmost dimension of char variable
                if len(data.shape) > 0 and data.shape[-1] == self.shape[-1]:
                    # also make sure slice is along last dimension
                    matchdim = True
                    for cnt in count:
                        if cnt[-1] != self.shape[-1]:
                            matchdim = False
                            break
                    if matchdim:
                        data = chartostring(data, encoding=encoding)

        # if structure array contains char arrays, return view as strings
        # if _Encoding att set (issue #773)
        if self._iscompound and \
           self._cmptype.dtype != self._cmptype.dtype_view and \
           self.chartostring:
#          self.chartostring and getattr(self,'_Encoding',None) is not None:
                data = data.view(self._cmptype.dtype_view)
        return data

    def _toma(self,data):
        cdef int ierr, no_fill
        # if attribute _Unsigned is "true", and variable has signed integer
        # dtype, return view with corresponding unsigned dtype (issues #656,
        # #794)
        # _Unsigned attribute must be "true" or "True" (string). Issue #1232.
        is_unsigned = getattr(self, '_Unsigned', False) in ["True","true"]
        is_unsigned_int = is_unsigned and data.dtype.kind == 'i'
        if self.scale and is_unsigned_int:  # only do this if autoscale option is on.
            dtype_unsigned_int='%su%s' % (data.dtype.byteorder,data.dtype.itemsize)
            data = data.view(dtype_unsigned_int)
        # private function for creating a masked array, masking missing_values
        # and/or _FillValues.
        totalmask = numpy.zeros(data.shape, numpy.bool_)
        fill_value = None
        safe_missval = self._check_safecast('missing_value')
        if safe_missval:
            mval = numpy.array(self.missing_value, self.dtype)
            if self.scale and is_unsigned_int:
                mval = mval.view(dtype_unsigned_int)
            # create mask from missing values.
            mvalmask = numpy.zeros(data.shape, numpy.bool_)
            if mval.shape == (): # mval a scalar.
                mval = [mval] # make into iterable.
            for m in mval:
                # is scalar missing value a NaN?
                try:
                    mvalisnan = numpy.isnan(m)
                except TypeError: # isnan fails on some dtypes (issue 206)
                    mvalisnan = False
                if mvalisnan:
                    mvalmask += numpy.isnan(data)
                else:
                    mvalmask += data==m
            if mvalmask.any():
                # set fill_value for masked array
                # to missing_value (or 1st element
                # if missing_value is a vector).
                fill_value = mval[0]
                totalmask += mvalmask
        # set mask=True for data == fill value
        safe_fillval = self._check_safecast('_FillValue')
        if safe_fillval:
            fval = numpy.array(self._FillValue, self.dtype)
            if self.scale and is_unsigned_int:
                fval = fval.view(dtype_unsigned_int)
            # is _FillValue a NaN?
            try:
                fvalisnan = numpy.isnan(fval)
            except: # isnan fails on some dtypes (issue 202)
                fvalisnan = False
            if fvalisnan:
                mask = numpy.isnan(data)
            elif (data == fval).any():
                mask = data==fval
            else:
                mask = None
            if mask is not None:
                if fill_value is None:
                    fill_value = fval
                totalmask += mask
        # issue 209: don't return masked array if variable filling
        # is disabled.
        else:
            if __netcdf4libversion__ < '4.5.1' and\
                self._grp.file_format.startswith('NETCDF3'):
                # issue #908: no_fill not correct for NETCDF3 files before 4.5.1
                # before 4.5.1 there was no way to turn off filling on a
                # per-variable basis for classic files.
                no_fill=0
            else:
                with nogil:
                    ierr = nc_inq_var_fill(self._grpid,self._varid,&no_fill,NULL)
                _ensure_nc_success(ierr)
            # if no_fill is not 1, and not a byte variable, then use default fill value.
            # from http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-c/Fill-Values.html#Fill-Values
            # "If you need a fill value for a byte variable, it is recommended
            # that you explicitly define an appropriate _FillValue attribute, as
            # generic utilities such as ncdump will not assume a default fill
            # value for byte variables."
            # Explained here too:
            # http://www.unidata.ucar.edu/software/netcdf/docs/known_problems.html#ncdump_ubyte_fill
            # "There should be no default fill values when reading any byte
            # type, signed or unsigned, because the byte ranges are too
            # small to assume one of the values should appear as a missing
            # value unless a _FillValue attribute is set explicitly."
            # (do this only for non-vlens, since vlens don't have a default _FillValue)
            if not self._isvlen and (no_fill != 1 or self.dtype.str[1:] not in ['u1','i1']):
                fillval = numpy.array(default_fillvals[self.dtype.str[1:]],self.dtype)
                has_fillval = data == fillval
                # if data is an array scalar, has_fillval will be a boolean.
                # in that case convert to an array.
                if type(has_fillval) == bool: has_fillval=numpy.asarray(has_fillval)
                if has_fillval.any():
                    if fill_value is None:
                        fill_value = fillval
                    mask=data==fillval
                    totalmask += mask
        # set mask=True for data outside valid_min,valid_max.
        # (issue #576)
        validmin = None; validmax = None
        # if valid_range exists use that, otherwise
        # look for valid_min, valid_max.  No special
        # treatment of byte data as described at
        # http://www.unidata.ucar.edu/software/netcdf/docs/attribute_conventions.html).
        safe_validrange = self._check_safecast('valid_range')
        safe_validmin = self._check_safecast('valid_min')
        safe_validmax = self._check_safecast('valid_max')
        if safe_validrange and self.valid_range.size == 2:
            validmin = numpy.array(self.valid_range[0], self.dtype)
            validmax = numpy.array(self.valid_range[1], self.dtype)
        else:
            if safe_validmin:
                validmin = numpy.array(self.valid_min, self.dtype)
            if safe_validmax:
                validmax = numpy.array(self.valid_max, self.dtype)
        if validmin is not None and self.scale and is_unsigned_int:
            validmin = validmin.view(dtype_unsigned_int)
        if validmax is not None and self.scale and is_unsigned_int:
            validmax = validmax.view(dtype_unsigned_int)
        # http://www.unidata.ucar.edu/software/netcdf/docs/attribute_conventions.html).
        # "If the data type is byte and _FillValue
        # is not explicitly defined,
        # then the valid range should include all possible values.
        # Otherwise, the valid range should exclude the _FillValue
        # (whether defined explicitly or by default) as follows.
        # If the _FillValue is positive then it defines a valid maximum,
        #  otherwise it defines a valid minimum."
        byte_type = self.dtype.str[1:] in ['u1','i1']
        if safe_fillval:
            fval = numpy.array(self._FillValue, self.dtype)
        else:
            fval = numpy.array(default_fillvals[self.dtype.str[1:]],self.dtype)
            if byte_type: fval = None
        if self.dtype.kind != 'S': # don't set mask for character data
            # issues #761 and #748:  setting valid_min/valid_max to the
            # _FillVaue is too surprising for many users (despite the
            # netcdf docs attribute best practices suggesting clients
            # should do this).
            #if validmin is None and (fval is not None and fval <= 0):
            #    validmin = fval
            #if validmax is None and (fval is not None and fval > 0):
            #    validmax = fval
            if validmin is not None:
                totalmask += data < validmin
            if validmax is not None:
                totalmask += data > validmax
        if fill_value is None and fval is not None:
            fill_value = fval
        # if all else fails, use default _FillValue as fill_value
        # for masked array.
        if fill_value is None:
            fill_value = default_fillvals[self.dtype.str[1:]]
        # create masked array with computed mask
        masked_values = bool(totalmask.any())
        if masked_values:
            data = ma.masked_array(data,mask=totalmask,fill_value=fill_value)
        else:
            # issue #785: always return masked array, if no values masked
            data = ma.masked_array(data)
        # issue 515 scalar array with mask=True should be converted
        # to numpy.ma.MaskedConstant to be consistent with slicing
        # behavior of masked arrays.
        if data.shape == () and data.mask.all():
            # return a scalar numpy masked constant not a 0-d masked array,
            # so that data == numpy.ma.masked.
            data = data[()] # changed from [...] (issue #662)
        elif not self.always_mask and not masked_values:
            # issue #809: return a regular numpy array if requested
            # and there are no missing values
            data = numpy.array(data, copy=False)

        return data

    def _pack(self,data):
        # pack non-masked values using scale_factor and add_offset
        if hasattr(self, 'scale_factor') and hasattr(self, 'add_offset'):
            data = (data - self.add_offset)/self.scale_factor
            if self.dtype.kind in 'iu': data = numpy.around(data)
        elif hasattr(self, 'scale_factor'):
            data = data/self.scale_factor
            if self.dtype.kind in 'iu': data = numpy.around(data)
        elif hasattr(self, 'add_offset'):
            data = data - self.add_offset
            if self.dtype.kind in 'iu': data = numpy.around(data)
        if ma.isMA(data):
            # if underlying data in masked regions of masked array
            # corresponds to missing values, don't fill masked array -
            # just use underlying data instead
            if hasattr(self, 'missing_value') and \
               numpy.all(numpy.in1d(data.data[data.mask],self.missing_value)):
                data = data.data
            else:
                if hasattr(self, 'missing_value'):
                    # if missing value is a scalar, use it as fill_value.
                    # if missing value is a vector, raise an exception
                    # since we then don't know how to fill in masked values.
                    if numpy.array(self.missing_value).shape == ():
                        fillval = self.missing_value
                    else:
                        msg="cannot assign fill_value for masked array when missing_value attribute is not a scalar"
                        raise RuntimeError(msg)
                    if numpy.array(fillval).shape != ():
                        fillval = fillval[0]
                elif hasattr(self, '_FillValue'):
                    fillval = self._FillValue
                else:
                    fillval = default_fillvals[self.dtype.str[1:]]
                # some versions of numpy have trouble handling
                # MaskedConstants when filling - this is is
                # a workaround (issue #850)
                if data.shape == (1,) and data.mask.all():
                    data = numpy.array([fillval],self.dtype)
                else:
                    data = data.filled(fill_value=fillval)
        if self.dtype != data.dtype:
            data = data.astype(self.dtype) # cast data to var type, if necessary.
        return data

    def _assign_vlen(self, elem, data):
        """private method to assign data to a single item in a VLEN variable"""
        cdef size_t *startp
        cdef size_t *countp
        cdef int ndims, n
        cdef nc_vlen_t *vldata
        cdef char **strdata
        cdef ndarray data2
        if not self._isvlen:
            raise TypeError('_assign_vlen method only for use with VLEN variables')
        ndims = self.ndim
        msg="data can only be assigned to VLEN variables using integer indices"
        # check to see that elem is a tuple of integers.
        # handle negative integers.
        if _is_int(elem):
            if ndims > 1:
                raise IndexError(msg)
            if elem < 0:
                if self.shape[0]+elem >= 0:
                    elem = self.shape[0]+elem
                else:
                    raise IndexError("Illegal index")
        elif isinstance(elem, tuple):
            if len(elem) != ndims:
                raise IndexError("Illegal index")
            elemnew = []
            for n,e in enumerate(elem):
                if not _is_int(e):
                    raise IndexError(msg)
                elif e < 0:
                    enew = self.shape[n]+e
                    if enew < 0:
                        raise IndexError("Illegal index")
                    else:
                        elemnew.append(self.shape[n]+e)
                else:
                    elemnew.append(e)
            elem = tuple(elemnew)
        else:
            raise IndexError(msg)
        # set start, count
        if isinstance(elem, tuple):
            start = list(elem)
        else:
            start = [elem]
        count = [1]*ndims
        startp = <size_t *>malloc(sizeof(size_t) * ndims)
        countp = <size_t *>malloc(sizeof(size_t) * ndims)
        for n from 0 <= n < ndims:
            startp[n] = start[n]
            countp[n] = count[n]
        if self.dtype == str: # VLEN string
            strdata = <char **>malloc(sizeof(char *))
            # use _Encoding attribute to specify string encoding - if
            # not given, use 'utf-8'.
            encoding = getattr(self,'_Encoding','utf-8')
            bytestr = _strencode(data,encoding=encoding)
            strdata[0] = bytestr
            with nogil:
                ierr = nc_put_vara(self._grpid, self._varid,
                                   startp, countp, strdata)
            _ensure_nc_success(ierr)
            free(strdata)
        else: # regular VLEN
            if data.dtype != self.dtype:
                raise TypeError("wrong data type: should be %s, got %s" % (self.dtype,data.dtype))
            data2 = data
            vldata = <nc_vlen_t *>malloc(sizeof(nc_vlen_t))
            vldata[0].len = PyArray_SIZE(data2)
            vldata[0].p = PyArray_DATA(data2)
            with nogil:
                ierr = nc_put_vara(self._grpid, self._varid,
                                   startp, countp, vldata)
            _ensure_nc_success(ierr)
            free(vldata)
        free(startp)
        free(countp)

    def _check_safecast(self, attname):
        # check to see that variable attribute exists
        # can can be safely cast to variable data type.
        msg="""WARNING: %s not used since it
cannot be safely cast to variable data type""" % attname
        if hasattr(self, attname):
            att = numpy.array(self.getncattr(attname))
        else:
            return False
        try:
            atta = numpy.array(att, self.dtype)
        except ValueError:
            is_safe = False
            warnings.warn(msg)
            return is_safe
        is_safe = _safecast(att,atta)
        if not is_safe:
            warnings.warn(msg)
        return is_safe

    def __setitem__(self, elem, data):
        # This special method is used to assign to the netCDF variable
        # using "extended slice syntax". The extended slice syntax
        # is a perfect match for the "start", "count" and "stride"
        # arguments to the nc_put_var() function, and is much more easy
        # to use.

        # if _Encoding is specified for a character variable, convert
        # numpy array of strings to a numpy array of characters with one more
        # dimension.
        if self.chartostring and getattr(self.dtype,'kind',None) == 'S' and\
           getattr(self.dtype,'itemsize',None) == 1:
            # NC_CHAR variable
            encoding = getattr(self,'_Encoding',None)
            if encoding is not None:
                # _Encoding attribute is set
                # if data is a string or a bytes object, convert to a numpy string array
                # whose length is equal to the rightmost dimension of the
                # variable.
                if type(data) in [str,bytes]: data = numpy.asarray(data,dtype='S'+repr(self.shape[-1]))
                if data.dtype.kind in ['S','U'] and data.dtype.itemsize > 1:
                    # if data is a numpy string array, convert it to an array
                    # of characters with one more dimension.
                    data = stringtochar(data, encoding=encoding)

        # if structured data has strings (and _Encoding att set), create view as char arrays
        # (issue #773)
        if self._iscompound and \
           self._cmptype.dtype != self._cmptype.dtype_view and \
           _set_viewdtype(data.dtype) == self._cmptype.dtype_view and \
           self.chartostring:
#          self.chartostring and getattr(self,'_Encoding',None) is not None:
                # may need to cast input data to aligned type
                data = data.astype(self._cmptype.dtype_view).view(self._cmptype.dtype)

        if self._isvlen: # if vlen, should be object array (don't try casting)
            if self.dtype == str:
                # for string vars, if data is not an array
                # assume it is a python string and raise an error
                # if it is an array, but not an object array.
                if not isinstance(data, numpy.ndarray):
                    # issue 458, allow Ellipsis to be used for scalar var
                    if type(elem) == type(Ellipsis) and not\
                       len(self.dimensions): elem = 0
                    self._assign_vlen(elem, data)
                    return
                elif data.dtype.kind in ['S', 'U']:
                    if ma.isMA(data):
                        msg='masked arrays cannot be assigned by VLEN str slices'
                        raise TypeError(msg)
                    data = data.astype(object)
                elif data.dtype.kind != 'O':
                    msg = ('only numpy string, unicode or object arrays can '
                           'be assigned to VLEN str var slices')
                    raise TypeError(msg)
            else:
                # for non-string vlen arrays, if data is not multi-dim, or
                # not an object array, assume it represents a single element
                # of the vlen var.
                if not isinstance(data, numpy.ndarray) or data.dtype.kind != 'O':
                    # issue 458, allow Ellipsis to be used for scalar var
                    if type(elem) == type(Ellipsis) and not\
                       len(self.dimensions): elem = 0
                    # pack as integers if desired.
                    if self.scale:
                        data = self._pack(data)
                    self._assign_vlen(elem, data)
                    return

        # A numpy or masked array (or an object supporting the buffer interface) is needed.
        # Convert if necessary.
        if not ma.isMA(data) and not (hasattr(data,'data') and isinstance(data.data,memoryview)):
            # if auto scaling is to be done, don't cast to an integer yet.
            if self.scale and self.dtype.kind in 'iu' and \
               hasattr(self, 'scale_factor') or hasattr(self, 'add_offset'):
                data = numpy.array(data,numpy.float64)
            else:
                data = numpy.array(data,self.dtype)

        # for Enum variable, make sure data is valid.
        if self._isenum:
            test = numpy.zeros(data.shape,numpy.bool_)
            if ma.isMA(data):
                # fix for new behaviour in numpy.ma in 1.13 (issue #662)
                for val in self.datatype.enum_dict.values():
                    test += data.filled() == val
            else:
                for val in self.datatype.enum_dict.values():
                    test += data == val
            if not numpy.all(test):
                msg="trying to assign illegal value to Enum variable"
                raise ValueError(msg)

        start, count, stride, put_ind =\
        _StartCountStride(elem,self.shape,self.dimensions,self._grp,datashape=data.shape,put=True,use_get_vars=self._use_get_vars)
        datashape = _out_array_shape(count)

        # if a numpy scalar, create an array of the right size
        # and fill with scalar values.
        if data.shape == ():
            data = numpy.tile(data,datashape)
        # reshape data array if needed to conform with start,count,stride.
        if data.ndim != len(datashape) or\
           (data.shape != datashape and data.ndim > 1): # issue #1083
            # create a view so shape in caller is not modified (issue 90)
            try: # if extra singleton dims, just reshape
                data = data.view()
                data.shape = tuple(datashape)
            except ValueError: # otherwise broadcast
                data = numpy.broadcast_to(data, datashape)

        # Reshape these arrays so we can iterate over them.
        start = start.reshape((-1, self.ndim or 1))
        count = count.reshape((-1, self.ndim or 1))
        stride = stride.reshape((-1, self.ndim or 1))
        put_ind = put_ind.reshape((-1, self.ndim or 1))

        # quantize data if least_significant_digit attribute
        # exists (improves compression).
        if self._has_lsd:
            data = _quantize(data,self.least_significant_digit)

        if self.scale and self._isprimitive:
            # pack non-masked values using scale_factor and add_offset
            data = self._pack(data)

        # Fill output array with data chunks.
        for (a,b,c,i) in zip(start, count, stride, put_ind):
            dataput = data[tuple(i)]
            if dataput.size == 0: continue # nothing to write
            # convert array scalar to regular array with one element.
            if dataput.shape == ():
                if self._isvlen:
                    dataput=numpy.array(dataput,'O')
                else:
                    dataput=numpy.array(dataput,dataput.dtype)
            self._put(dataput,a,b,c)


    def __len__(self):
        if not self.shape:
            raise TypeError('len() of unsized object')
        else:
            return self.shape[0]


    def assignValue(self,val):
        """
**`assignValue(self, val)`**

assign a value to a scalar variable.  Provided for compatibility with
Scientific.IO.NetCDF, can also be done by assigning to an Ellipsis slice ([...])."""
        if len(self.dimensions):
            raise IndexError('to assign values to a non-scalar variable, use a slice')
        self[:]=val

    def getValue(self):
        """
**`getValue(self)`**

get the value of a scalar variable.  Provided for compatibility with
Scientific.IO.NetCDF, can also be done by slicing with an Ellipsis ([...])."""
        if len(self.dimensions):
            raise IndexError('to retrieve values from a non-scalar variable, use slicing')
        return self[slice(None)]

    def set_auto_chartostring(self,chartostring):
        """
**`set_auto_chartostring(self,chartostring)`**

turn on or off automatic conversion of character variable data to and
from numpy fixed length string arrays when the `_Encoding` variable attribute
is set.

If `chartostring` is set to `True`, when data is read from a character variable
(dtype = `S1`) that has an `_Encoding` attribute, it is converted to a numpy
fixed length unicode string array (dtype = `UN`, where `N` is the length
of the the rightmost dimension of the variable).  The value of `_Encoding`
is the unicode encoding that is used to decode the bytes into strings.

When numpy string data is written to a variable it is converted back to
indiviual bytes, with the number of bytes in each string equalling the
rightmost dimension of the variable.

The default value of `chartostring` is `True`
(automatic conversions are performed).
        """
        self.chartostring = bool(chartostring)

    def use_nc_get_vars(self,use_nc_get_vars):
        """
**`use_nc_get_vars(self,_use_get_vars)`**

enable the use of netcdf library routine `nc_get_vars`
to retrieve strided variable slices.  By default,
`nc_get_vars` may not used by default (depending on the
version of the netcdf-c library being used) since it may be
slower than multiple calls to the unstrided read routine `nc_get_vara`.
        """
        self._use_get_vars = bool(use_nc_get_vars)

    def set_auto_maskandscale(self,maskandscale):
        """
**`set_auto_maskandscale(self,maskandscale)`**

turn on or off automatic conversion of variable data to and
from masked arrays, automatic packing/unpacking of variable
data using `scale_factor` and `add_offset` attributes and
automatic conversion of signed integer data to unsigned integer
data if the `_Unsigned` attribute exists and is set to "true" (or "True").

If `maskandscale` is set to `True`, when data is read from a variable
it is converted to a masked array if any of the values are exactly
equal to the either the netCDF _FillValue or the value specified by the
missing_value variable attribute. The fill_value of the masked array
is set to the missing_value attribute (if it exists), otherwise
the netCDF _FillValue attribute (which has a default value
for each data type). If the variable has no missing_value attribute, the
_FillValue is used instead. If the variable has valid_min/valid_max and
missing_value attributes, data outside the specified range will be masked.
When data is written to a variable, the masked
array is converted back to a regular numpy array by replacing all the
masked values by the missing_value attribute of the variable (if it
exists).  If the variable has no missing_value attribute, the _FillValue
is used instead. 

If `maskandscale` is set to `True`, and the variable has a
`scale_factor` or an `add_offset` attribute, then data read
from that variable is unpacked using::

    data = self.scale_factor*data + self.add_offset

When data is written to a variable it is packed using::

    data = (data - self.add_offset)/self.scale_factor

If either scale_factor is present, but add_offset is missing, add_offset
is assumed zero.  If add_offset is present, but scale_factor is missing,
scale_factor is assumed to be one.
For more information on how `scale_factor` and `add_offset` can be
used to provide simple compression, see the
[PSL metadata conventions](http://www.esrl.noaa.gov/psl/data/gridded/conventions/cdc_netcdf_standard.shtml).

In addition, if `maskandscale` is set to `True`, and if the variable has an
attribute `_Unsigned` set to "true", and the variable has a signed integer data type,
a view to the data is returned with the corresponding unsigned integer data type.
This convention is used by the netcdf-java library to save unsigned integer
data in `NETCDF3` or `NETCDF4_CLASSIC` files (since the `NETCDF3`
data model does not have unsigned integer data types).

The default value of `maskandscale` is `True`
(automatic conversions are performed).
        """
        self.scale = self.mask = bool(maskandscale)

    def set_auto_scale(self,scale):
        """
**`set_auto_scale(self,scale)`**

turn on or off automatic packing/unpacking of variable
data using `scale_factor` and `add_offset` attributes.
Also turns on and off automatic conversion of signed integer data
to unsigned integer data if the variable has an `_Unsigned`
attribute set to "true" or "True".

If `scale` is set to `True`, and the variable has a
`scale_factor` or an `add_offset` attribute, then data read
from that variable is unpacked using::

    data = self.scale_factor*data + self.add_offset

When data is written to a variable it is packed using::

    data = (data - self.add_offset)/self.scale_factor

If either scale_factor is present, but add_offset is missing, add_offset
is assumed zero.  If add_offset is present, but scale_factor is missing,
scale_factor is assumed to be one.
For more information on how `scale_factor` and `add_offset` can be
used to provide simple compression, see the
[PSL metadata conventions](http://www.esrl.noaa.gov/psl/data/gridded/conventions/cdc_netcdf_standard.shtml).

In addition, if `scale` is set to `True`, and if the variable has an
attribute `_Unsigned` set to "true", and the variable has a signed integer data type,
a view to the data is returned with the corresponding unsigned integer datatype.
This convention is used by the netcdf-java library to save unsigned integer
data in `NETCDF3` or `NETCDF4_CLASSIC` files (since the `NETCDF3`
data model does not have unsigned integer data types).

The default value of `scale` is `True`
(automatic conversions are performed).
        """
        self.scale = bool(scale)

    def set_auto_mask(self,mask):
        """
**`set_auto_mask(self,mask)`**

turn on or off automatic conversion of variable data to and
from masked arrays .

If `mask` is set to `True`, when data is read from a variable
it is converted to a masked array if any of the values are exactly
equal to the either the netCDF _FillValue or the value specified by the
missing_value variable attribute. The fill_value of the masked array
is set to the missing_value attribute (if it exists), otherwise
the netCDF _FillValue attribute (which has a default value
for each data type). If the variable has no missing_value attribute, the
_FillValue is used instead. If the variable has valid_min/valid_max and
missing_value attributes, data outside the specified range will be masked.
When data is written to a variable, the masked
array is converted back to a regular numpy array by replacing all the
masked values by the missing_value attribute of the variable (if it
exists).  If the variable has no missing_value attribute, the _FillValue
is used instead. 

The default value of `mask` is `True`
(automatic conversions are performed).
        """
        self.mask = bool(mask)

    def set_always_mask(self,always_mask):
        """
**`set_always_mask(self,always_mask)`**

turn on or off conversion of data without missing values to regular
numpy arrays.

`always_mask` is a Boolean determining if automatic conversion of
masked arrays with no missing values to regular numpy arrays shall be
applied. Default is True. Set to False to restore the default behaviour
in versions prior to 1.4.1 (numpy array returned unless missing values are present,
otherwise masked array returned).
        """
        self.always_mask = bool(always_mask)

    def set_ncstring_attrs(self,ncstring_attrs):
        """
**`set_always_mask(self,ncstring_attrs)`**

turn on or off creating NC_STRING string attributes.

If `ncstring_attrs` is set to `True` then text attributes will be variable-length
NC_STRINGs.

The default value of `ncstring_attrs` is `False` (writing ascii text attributes as
NC_CHAR).

        """
        self._ncstring_attrs__ = bool(ncstring_attrs)

    def _put(self,ndarray data,start,count,stride):
        """Private method to put data into a netCDF variable"""
        cdef int ierr, ndims
        cdef npy_intp totelem
        cdef size_t *startp
        cdef size_t *countp
        cdef ptrdiff_t *stridep
        cdef char **strdata
        cdef void* elptr
        cdef char* databuff
        cdef ndarray dataarr
        cdef nc_vlen_t *vldata
        # rank of variable.
        ndims = len(self.dimensions)
        # make sure data is contiguous.
        # if not, make a local copy.
        if not PyArray_ISCONTIGUOUS(data):
            data = data.copy()
        # fill up startp,countp,stridep.
        totelem = 1
        negstride = 0
        sl = []
        startp = <size_t *>malloc(sizeof(size_t) * ndims)
        countp = <size_t *>malloc(sizeof(size_t) * ndims)
        stridep = <ptrdiff_t *>malloc(sizeof(ptrdiff_t) * ndims)
        for n from 0 <= n < ndims:
            count[n] = abs(count[n]) # make -1 into +1
            countp[n] = count[n]
            # for neg strides, reverse order (then flip that axis after data read in)
            if stride[n] < 0:
                negstride = 1
                stridep[n] = -stride[n]
                startp[n] = start[n]+stride[n]*(count[n]-1)
                stride[n] = -stride[n]
                sl.append(slice(None, None, -1)) # this slice will reverse the data
            else:
                startp[n] = start[n]
                stridep[n] = stride[n]
                sl.append(slice(None,None, 1))
            totelem = totelem*countp[n]
        # check to see that size of data array is what is expected
        # for slice given.
        dataelem = PyArray_SIZE(data)
        if totelem != dataelem:
            raise IndexError('size of data array does not conform to slice')
        if negstride:
            # reverse data along axes with negative strides.
            data = data[tuple(sl)].copy() # make sure a copy is made.
        if self._isprimitive or self._iscompound or self._isenum:
            # primitive, enum or compound data type.
            # if data type of array doesn't match variable,
            # try to cast the data.
            if self.dtype != data.dtype:
                data = data.astype(self.dtype) # cast data, if necessary.
            # byte-swap data in numpy array so that is has native
            # endian byte order (this is what netcdf-c expects -
            # issue #554, pull request #555)
            if not data.dtype.isnative:
                data = data.byteswap()
            # strides all 1 or scalar variable, use put_vara (faster)
            if sum(stride) == ndims or ndims == 0:
                with nogil:
                    ierr = nc_put_vara(self._grpid, self._varid,
                                       startp, countp, PyArray_DATA(data))
            else:
                with nogil:
                    ierr = nc_put_vars(self._grpid, self._varid,
                                       startp, countp, stridep, PyArray_DATA(data))
            _ensure_nc_success(ierr)
        elif self._isvlen:
            if data.dtype.char !='O':
                raise TypeError('data to put in string variable must be an object array containing Python strings')
            # flatten data array.
            data = data.flatten()
            if self.dtype == str:
                # convert all elements from strings to bytes
                # use _Encoding attribute to specify string encoding - if
                # not given, use 'utf-8'.
                encoding = getattr(self,'_Encoding','utf-8')
                for n in range(data.shape[0]):
                    data[n] = _strencode(data[n],encoding=encoding)
                # vlen string (NC_STRING)
                # loop over elements of object array, put data buffer for
                # each element in struct.
                # allocate struct array to hold vlen data.
                strdata = <char **>malloc(sizeof(char *)*totelem)
                for i from 0<=i<totelem:
                    strdata[i] = data[i]
                # strides all 1 or scalar variable, use put_vara (faster)
                if sum(stride) == ndims or ndims == 0:
                    with nogil:
                        ierr = nc_put_vara(self._grpid, self._varid,
                                           startp, countp, strdata)
                else:
                    raise IndexError('strides must all be 1 for string variables')
                    #with nogil:
                    #    ierr = nc_put_vars(self._grpid, self._varid,
                    #                       startp, countp, stridep, strdata)
                _ensure_nc_success(ierr)
                free(strdata)
            else:
                # regular vlen.
                # loop over elements of object array, put data buffer for
                # each element in struct.
                databuff = PyArray_BYTES(<ndarray>data)
                # allocate struct array to hold vlen data.
                vldata = <nc_vlen_t *>malloc(<size_t>totelem*sizeof(nc_vlen_t))
                for i from 0<=i<totelem:
                    elptr = (<void**>databuff)[0]
                    dataarr = <ndarray>elptr
                    if self.dtype != dataarr.dtype.str[1:]:
                        #dataarr = dataarr.astype(self.dtype) # cast data, if necessary.
                        # casting doesn't work ?? just raise TypeError
                        raise TypeError("wrong data type in object array: should be %s, got %s" % (self.dtype,dataarr.dtype))
                    vldata[i].len = PyArray_SIZE(dataarr)
                    vldata[i].p = PyArray_DATA(dataarr)
                    databuff = databuff + PyArray_STRIDES(data)[0]
                # strides all 1 or scalar variable, use put_vara (faster)
                if sum(stride) == ndims or ndims == 0:
                    with nogil:
                        ierr = nc_put_vara(self._grpid, self._varid,
                                           startp, countp, vldata)
                else:
                    raise IndexError('strides must all be 1 for vlen variables')
                    #with nogil:
                    #    ierr = nc_put_vars(self._grpid, self._varid,
                    #                       startp, countp, stridep, vldata)
                _ensure_nc_success(ierr)
                # free the pointer array.
                free(vldata)
        free(startp)
        free(countp)
        free(stridep)

    def _get(self,start,count,stride):
        """Private method to retrieve data from a netCDF variable"""
        cdef int ierr, ndims, totelem
        cdef size_t *startp
        cdef size_t *countp
        cdef ptrdiff_t *stridep
        cdef ndarray data, dataarr
        cdef void *elptr
        cdef char **strdata
        cdef nc_vlen_t *vldata
        # if one of the counts is negative, then it is an index
        # and not a slice so the resulting array
        # should be 'squeezed' to remove the singleton dimension.
        shapeout = ()
        squeeze_out = False
        for lendim in count:
            if lendim == -1:
                shapeout = shapeout + (1,)
                squeeze_out = True
            else:
                shapeout = shapeout + (lendim,)
        # rank of variable.
        ndims = len(self.dimensions)
        # fill up startp,countp,stridep.
        negstride = 0
        sl = []
        startp = <size_t *>malloc(sizeof(size_t) * ndims)
        countp = <size_t *>malloc(sizeof(size_t) * ndims)
        stridep = <ptrdiff_t *>malloc(sizeof(ptrdiff_t) * ndims)
        for n from 0 <= n < ndims:
            count[n] = abs(count[n]) # make -1 into +1
            countp[n] = count[n]
            # for neg strides, reverse order (then flip that axis after data read in)
            if stride[n] < 0:
                negstride = 1
                stridep[n] = -stride[n]
                startp[n] = start[n]+stride[n]*(count[n]-1)
                stride[n] = -stride[n]
                sl.append(slice(None, None, -1)) # this slice will reverse the data
            else:
                startp[n] = start[n]
                stridep[n] = stride[n]
                sl.append(slice(None,None, 1))
        if self._isprimitive or self._iscompound or self._isenum:
            data = numpy.empty(shapeout, self.dtype)
            # strides all 1 or scalar variable, use get_vara (faster)
            # if count contains a zero element, no data is being read
            if 0 not in count:
                if sum(stride) == ndims or ndims == 0:
                    with nogil:
                        ierr = nc_get_vara(self._grpid, self._varid,
                                           startp, countp, PyArray_DATA(data))
                else:
                    with nogil:
                        ierr = nc_get_vars(self._grpid, self._varid,
                                           startp, countp, stridep,
                                           PyArray_DATA(data))
            else:
                ierr = 0
            if ierr == NC_EINVALCOORDS:
                raise IndexError('index exceeds dimension bounds')
            elif ierr != NC_NOERR:
                _ensure_nc_success(ierr)
        elif self._isvlen:
            # allocate array of correct primitive type.
            data = numpy.empty(shapeout, 'O')
            # flatten data array.
            data = data.flatten()
            totelem = PyArray_SIZE(data)
            if self.dtype == str:
                # vlen string (NC_STRING)
                # allocate pointer array to hold string data.
                strdata = <char **>malloc(sizeof(char *) * totelem)
                # strides all 1 or scalar variable, use get_vara (faster)
                if sum(stride) == ndims or ndims == 0:
                    with nogil:
                        ierr = nc_get_vara(self._grpid, self._varid,
                                           startp, countp, strdata)
                else:
                    # FIXME: is this a bug in netCDF4?
                    raise IndexError('strides must all be 1 for string variables')
                    #with nogil:
                    #    ierr = nc_get_vars(self._grpid, self._varid,
                    #                       startp, countp, stridep, strdata)
                if ierr == NC_EINVALCOORDS:
                    raise IndexError
                elif ierr != NC_NOERR:
                   _ensure_nc_success(ierr)
                # loop over elements of object array, fill array with
                # contents of strdata.
                # use _Encoding attribute to decode string to bytes - if
                # not given, use 'utf-8'.
                encoding = getattr(self,'_Encoding','utf-8')
                for i from 0<=i<totelem:
                    if strdata[i]:
                        data[i] = strdata[i].decode(encoding)
                    else:
                        data[i] = "" # issue 915
                # reshape the output array
                data = numpy.reshape(data, shapeout)
                # free string data internally allocated in netcdf C lib
                with nogil:
                    ierr = nc_free_string(totelem, strdata)
                # free the pointer array
                free(strdata)
            else:
                # regular vlen
                # allocate struct array to hold vlen data.
                vldata = <nc_vlen_t *>malloc(totelem*sizeof(nc_vlen_t))
                for i in range(totelem):
                    vldata[i].len = 0
                    vldata[i].p = <void*>0
                # strides all 1 or scalar variable, use get_vara (faster)
                if sum(stride) == ndims or ndims == 0:
                    with nogil:
                        ierr = nc_get_vara(self._grpid, self._varid,
                                           startp, countp, vldata)
                else:
                    raise IndexError('strides must all be 1 for vlen variables')
                    #with nogil:
                    #    ierr = nc_get_vars(self._grpid, self._varid,
                    #                       startp, countp, stridep, vldata)
                if ierr == NC_EINVALCOORDS:
                    raise IndexError
                elif ierr != NC_NOERR:
                    _ensure_nc_success(ierr)
                # loop over elements of object array, fill array with
                # contents of vlarray struct, put array in object array.
                for i from 0<=i<totelem:
                    arrlen  = vldata[i].len
                    dataarr = numpy.empty(arrlen, self.dtype)
                    #dataarr.data = <char *>vldata[i].p
                    memcpy(PyArray_DATA(dataarr), vldata[i].p, dataarr.nbytes)
                    data[i] = dataarr
                # reshape the output array
                data = numpy.reshape(data, shapeout)
                # free vlen data internally allocated in netcdf C lib
                with nogil:
                    ierr = nc_free_vlens(totelem, vldata)
                # free the pointer array
                free(vldata)
        free(startp)
        free(countp)
        free(stridep)
        if negstride:
            # reverse data along axes with negative strides.
            data = data[tuple(sl)].copy() # make a copy so data is contiguous.
        # netcdf-c always returns data in native byte order,
        # regardless of variable endian-ness. Here we swap the
        # bytes if the variable dtype is not native endian, so the
        # dtype of the returned numpy array matches the variable dtype.
        # (pull request #555, issue #554).
        if not data.dtype.isnative:
            data.byteswap(True) # in-place byteswap
        if not self.dimensions:
            return data[0] # a scalar
        elif squeeze_out:
            return numpy.squeeze(data)
        else:
            return data

    def set_collective(self, value):
        """
**`set_collective(self,True_or_False)`**

turn on or off collective parallel IO access. Ignored if file is not
open for parallel access.
        """
        IF HAS_PARALLEL4_SUPPORT or HAS_PNETCDF_SUPPORT:
            # set collective MPI IO mode on or off
            if value:
                with nogil:
                    ierr = nc_var_par_access(self._grpid, self._varid,
                           NC_COLLECTIVE)
            else:
                with nogil:
                    ierr = nc_var_par_access(self._grpid, self._varid,
                           NC_INDEPENDENT)
            _ensure_nc_success(ierr)
        ELSE:
            pass # does nothing

    def get_dims(self):
        """
**`get_dims(self)`**

return a tuple of `Dimension` instances associated with this
`Variable`.
        """
        return tuple(_find_dim(self._grp, dim) for dim in self.dimensions)

    def __reduce__(self):
        # raise error is user tries to pickle a Variable object.
        raise NotImplementedError('Variable is not picklable')

# Compound datatype support.

cdef class CompoundType:
    """
A `CompoundType` instance is used to describe a compound data
type, and can be passed to the the `Dataset.createVariable` method of
a `Dataset` or `Group` instance.
Compound data types map to numpy structured arrays.
See `CompoundType.__init__` for more details.

The instance variables `dtype` and `name` should not be modified by
the user.
    """
    cdef public nc_type _nc_type
    cdef public dtype, dtype_view, name
    def __init__(self, grp, object dt, object dtype_name, **kwargs):
        """
        ***`__init__(group, datatype, datatype_name)`***

        CompoundType constructor.

        **`group`**: `Group` instance to associate with the compound datatype.

        **`datatype`**: A numpy dtype object describing a structured (a.k.a record)
        array.  Can be composed of homogeneous numeric or character data types, or
        other structured array data types.

        **`datatype_name`**: a Python string containing a description of the
        compound data type.

        ***Note 1***: When creating nested compound data types,
        the inner compound data types must already be associated with CompoundType
        instances (so create CompoundType instances for the innermost structures
        first).

        ***Note 2***: `CompoundType` instances should be created using the
        `Dataset.createCompoundType` method of a `Dataset` or
        `Group` instance, not using this class directly.
        """
        cdef nc_type xtype
        # convert dt to a numpy datatype object
        # and make sure the isalignedstruct flag is set to True
        # (so padding is added to the fields to match what a
        # C compiler would output for a similar C-struct).
        # This is needed because nc_get_vara is
        # apparently expecting the data buffer to include
        # padding to match what a C struct would have.
        # (this may or may not be still true, but empirical
        # evidence suggests that segfaults occur if this
        # alignment step is skipped - see issue #705).
        # numpy string subdtypes (i.e. 'S80') are
        # automatically converted to character array
        # subtypes (i.e. ('S1',80)).  If '_Encoding'
        # variable attribute is set, data will be converted
        # to and from the string array representation with views.
        dt = _set_alignment(numpy.dtype(dt))
        # create a view datatype for converting char arrays to/from strings
        dtview = _set_viewdtype(numpy.dtype(dt))
        if 'typeid' in kwargs:
            xtype = kwargs['typeid']
        else:
            xtype = _def_compound(grp, dt, dtype_name)
        self._nc_type = xtype
        self.dtype = dt
        self.dtype_view = dtview
        self.name = dtype_name

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return "%r: name = '%s', numpy dtype = %s" %\
            (type(self), self.name, self.dtype)

    def __reduce__(self):
        # raise error is user tries to pickle a CompoundType object.
        raise NotImplementedError('CompoundType is not picklable')

def _set_alignment(dt):
    # recursively set alignment flag in nested structured data type
    names = dt.names; formats = []
    for name in names:
        fmt = dt.fields[name][0]
        if fmt.kind == 'V':
            if fmt.shape == ():
                dtx = _set_alignment(dt.fields[name][0])
            else:
                if fmt.subdtype[0].kind == 'V': # structured dtype
                    raise TypeError('nested structured dtype arrays not supported')
                else:
                    dtx = dt.fields[name][0]
        else:
            # convert character string elements to char arrays
            if fmt.kind == 'S' and fmt.itemsize != 1:
                dtx = numpy.dtype('(%s,)S1' % fmt.itemsize)
            else:
                # primitive data type
                dtx = dt.fields[name][0]
        formats.append(dtx)
    # leave out offsets, they will be re-computed to preserve alignment.
    dtype_dict = {'names':names,'formats':formats}
    return numpy.dtype(dtype_dict, align=True)

def _set_viewdtype(dt):
    # recursively change character array dtypes to string dtypes
    names = dt.names; formats = []
    for name in names:
        fmt = dt.fields[name][0]
        if fmt.kind == 'V':
            if fmt.shape == ():
                dtx = _set_viewdtype(dt.fields[name][0])
            else:
                if fmt.subdtype[0].kind == 'V': # structured dtype
                    raise TypeError('nested structured dtype arrays not supported')
                elif fmt.subdtype[0].kind == 'S' and len(dt.fields[name][0].shape) == 1:
                    lenchar = dt.fields[name][0].shape[0]
                    dtx = numpy.dtype('S%s' % lenchar)
                else:
                    dtx = dt.fields[name][0]
        else:
            # primitive data type
            dtx = dt.fields[name][0]
        formats.append(dtx)
    dtype_dict = {'names':names,'formats':formats}
    return numpy.dtype(dtype_dict, align=True)

cdef _def_compound(grp, object dt, object dtype_name):
    # private function used to construct a netcdf compound data type
    # from a numpy dtype object by CompoundType.__init__.
    cdef nc_type xtype, xtype_tmp
    cdef int ierr, ndims, grpid
    cdef size_t offset, size
    cdef char *namstring
    cdef char *nested_namstring
    cdef int *dim_sizes
    bytestr = _strencode(dtype_name)
    namstring = bytestr
    size = dt.itemsize
    grpid = grp._grpid
    with nogil:
        ierr = nc_def_compound(grpid, size, namstring, &xtype)
    _ensure_nc_success(ierr)
    names = list(dt.fields.keys())
    formats = [v[0] for v in dt.fields.values()]
    offsets = [v[1] for v in dt.fields.values()]
    # make sure entries in lists sorted by offset.
    # (don't know why this is necessary, but it is for version 4.0.1)
    names = _sortbylist(names, offsets)
    formats = _sortbylist(formats, offsets)
    offsets.sort()
    for name, format, offset in zip(names, formats, offsets):
        bytestr = _strencode(name)
        namstring = bytestr
        if format.kind != 'V': # scalar primitive type
            try:
                xtype_tmp = _nptonctype[format.str[1:]]
            except KeyError:
                raise ValueError('Unsupported compound type element')
            with nogil:
                ierr = nc_insert_compound(grpid, xtype, namstring,
                                          offset, xtype_tmp)
            _ensure_nc_success(ierr)
        else:
            if format.shape ==  (): # nested scalar compound type
                # find this compound type in this group or it's parents.
                xtype_tmp = _find_cmptype(grp, format)
                bytestr = _strencode(name)
                nested_namstring = bytestr
                with nogil:
                    ierr = nc_insert_compound(grpid, xtype,\
                                              nested_namstring,\
                                              offset, xtype_tmp)
                _ensure_nc_success(ierr)
            else: # nested array compound element
                ndims = len(format.shape)
                dim_sizes = <int *>malloc(sizeof(int) * ndims)
                for n from 0 <= n < ndims:
                    dim_sizes[n] = format.shape[n]
                if format.subdtype[0].kind != 'V': # primitive type.
                    try:
                        xtype_tmp = _nptonctype[format.subdtype[0].str[1:]]
                    except KeyError:
                        raise ValueError('Unsupported compound type element')
                    with nogil:
                        ierr = nc_insert_array_compound(grpid,xtype,namstring,
                           offset,xtype_tmp,ndims,dim_sizes)
                    _ensure_nc_success(ierr)
                else: # nested array compound type.
                    raise TypeError('nested structured dtype arrays not supported')
                    # this code is untested and probably does not work, disable
                    # for now...
                #   # find this compound type in this group or it's parents.
                #   xtype_tmp = _find_cmptype(grp, format.subdtype[0])
                #   bytestr = _strencode(name)
                #   nested_namstring = bytestr
                #   with nogil:
                #       ierr = nc_insert_array_compound(grpid,xtype,\
                #                                       nested_namstring,\
                #                                       offset,xtype_tmp,\
                #                                       ndims,dim_sizes)
                #   _ensure_nc_success(ierr)
                free(dim_sizes)
    return xtype

cdef _find_cmptype(grp, dtype):
    # look for data type in this group and it's parents.
    # return datatype id when found, if not found, raise exception.
    cdef nc_type xtype
    match = False
    for cmpname, cmpdt in grp.cmptypes.items():
        xtype = cmpdt._nc_type
        names1 = dtype.fields.keys()
        names2 = cmpdt.dtype.fields.keys()
        formats1 = [v[0] for v in dtype.fields.values()]
        formats2 = [v[0] for v in cmpdt.dtype.fields.values()]
        formats2v = [v[0] for v in cmpdt.dtype_view.fields.values()]
        # match names, formats, but not offsets (they may be changed
        # by netcdf lib).
        if names1==names2 and formats1==formats2 or (formats1 == formats2v):
            match = True
            break
    if not match:
        try:
            parent_grp = grp.parent
        except AttributeError:
            raise ValueError("cannot find compound type in this group or parent groups")
        if parent_grp is None:
            raise ValueError("cannot find compound type in this group or parent groups")
        else:
            xtype = _find_cmptype(parent_grp,dtype)
    return xtype

cdef _read_compound(group, nc_type xtype, endian=None):
    # read a compound data type id from an existing file,
    # construct a corresponding numpy dtype instance,
    # then use that to create a CompoundType instance.
    # called by _get_vars, _get_types and _get_att.
    # Calls itself recursively for nested compound types.
    cdef int ierr, nf, numdims, ndim, classp, _grpid
    cdef size_t nfields, offset
    cdef nc_type field_typeid
    cdef int *dim_sizes
    cdef char field_namstring[NC_MAX_NAME+1]
    cdef char cmp_namstring[NC_MAX_NAME+1]
    # get name and number of fields.
    _grpid = group._grpid
    with nogil:
        ierr = nc_inq_compound(_grpid, xtype, cmp_namstring, NULL, &nfields)
    _ensure_nc_success(ierr)
    name = cmp_namstring.decode('utf-8')
    # loop over fields.
    names = []
    formats = []
    offsets = []
    for nf from 0 <= nf < nfields:
        with nogil:
            ierr = nc_inq_compound_field(_grpid,
                                         xtype,
                                         nf,
                                         field_namstring,
                                         &offset,
                                         &field_typeid,
                                         &numdims,
                                         NULL)
        _ensure_nc_success(ierr)
        dim_sizes = <int *>malloc(sizeof(int) * numdims)
        with nogil:
            ierr = nc_inq_compound_field(_grpid,
                                         xtype,
                                         nf,
                                         field_namstring,
                                         &offset,
                                         &field_typeid,
                                         &numdims,
                                         dim_sizes)
        _ensure_nc_success(ierr)
        field_name = field_namstring.decode('utf-8')
        names.append(field_name)
        offsets.append(offset)
        # if numdims=0, not an array.
        field_shape = ()
        if numdims != 0:
            for ndim from 0 <= ndim < numdims:
                field_shape = field_shape + (dim_sizes[ndim],)
        free(dim_sizes)
        # check to see if this field is a nested compound type.
        try:
            field_type =  _nctonptype[field_typeid]
            if endian is not None:
                format = endian + format
        except KeyError:
            with nogil:
                ierr = nc_inq_user_type(_grpid,
                       field_typeid,NULL,NULL,NULL,NULL,&classp)
            if classp == NC_COMPOUND: # a compound type
                # recursively call this function?
                field_type = _read_compound(group, field_typeid, endian=endian)
            else:
                raise KeyError('compound field of an unsupported data type')
        if field_shape != ():
            formats.append((field_type,field_shape))
        else:
            formats.append(field_type)
    # make sure entries in lists sorted by offset.
    names = _sortbylist(names, offsets)
    formats = _sortbylist(formats, offsets)
    offsets.sort()
    # create a dict that can be converted into a numpy dtype.
    dtype_dict = {'names':names,'formats':formats,'offsets':offsets}
    return CompoundType(group, dtype_dict, name, typeid=xtype)

# VLEN datatype support.

cdef class VLType:
    """
A `VLType` instance is used to describe a variable length (VLEN) data
type, and can be passed to the the `Dataset.createVariable` method of
a `Dataset` or `Group` instance. See
`VLType.__init__` for more details.

The instance variables `dtype` and `name` should not be modified by
the user.
    """
    cdef public nc_type _nc_type
    cdef public dtype, name
    def __init__(self, grp, object dt, object dtype_name, **kwargs):
        """
        **`__init__(group, datatype, datatype_name)`**

        VLType constructor.

        **`group`**: `Group` instance to associate with the VLEN datatype.

        **`datatype`**: An numpy dtype object describing the component type for the
        variable length array.

        **`datatype_name`**: a Python string containing a description of the
        VLEN data type.

        ***`Note`***: `VLType` instances should be created using the
        `Dataset.createVLType` method of a `Dataset` or
        `Group` instance, not using this class directly.
        """
        cdef nc_type xtype
        if 'typeid' in kwargs:
            xtype = kwargs['typeid']
        else:
            xtype, dt = _def_vlen(grp, dt, dtype_name)
        self._nc_type = xtype
        self.dtype = dt
        if dt == str:
            self.name = None
        else:
            self.name = dtype_name

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        if self.dtype == str:
            return '%r: string type' % (type(self),)
        else:
            return "%r: name = '%s', numpy dtype = %s" %\
                (type(self), self.name, self.dtype)

    def __reduce__(self):
        # raise error is user tries to pickle a VLType object.
        raise NotImplementedError('VLType is not picklable')

cdef _def_vlen(grp, object dt, object dtype_name):
    # private function used to construct a netcdf VLEN data type
    # from a numpy dtype object or python str object by VLType.__init__.
    cdef nc_type xtype, xtype_tmp
    cdef int ierr, ndims, grpid
    cdef size_t offset, size
    cdef char *namstring
    cdef char *nested_namstring
    grpid = grp._grpid
    if dt == str: # python string, use NC_STRING
        xtype = NC_STRING
        # dtype_name ignored
    else: # numpy datatype
        bytestr = _strencode(dtype_name)
        namstring = bytestr
        dt = numpy.dtype(dt) # convert to numpy datatype.
        if dt.str[1:] in _supportedtypes:
            # find netCDF primitive data type corresponding to
            # specified numpy data type.
            xtype_tmp = _nptonctype[dt.str[1:]]
            with nogil:
                ierr = nc_def_vlen(grpid, namstring, xtype_tmp, &xtype);
            _ensure_nc_success(ierr)
        else:
            raise KeyError("unsupported datatype specified for VLEN")
    return xtype, dt

cdef _read_vlen(group, nc_type xtype, endian=None):
    # read a VLEN data type id from an existing file,
    # construct a corresponding numpy dtype instance,
    # then use that to create a VLType instance.
    # called by _get_types, _get_vars.
    cdef int ierr, grpid
    cdef size_t vlsize
    cdef nc_type base_xtype
    cdef char vl_namstring[NC_MAX_NAME+1]
    grpid = group._grpid
    if xtype == NC_STRING:
        dt = str
        name = None
    else:
        with nogil:
            ierr = nc_inq_vlen(grpid, xtype, vl_namstring, &vlsize, &base_xtype)
        _ensure_nc_success(ierr)
        name = vl_namstring.decode('utf-8')
        try:
            datatype = _nctonptype[base_xtype]
            if endian is not None: datatype = endian + datatype
            dt = numpy.dtype(datatype) # see if it is a primitive type
        except KeyError:
            raise KeyError("unsupported component type for VLEN")
    return VLType(group, dt, name, typeid=xtype)

# Enum datatype support.

cdef class EnumType:
    """
A `EnumType` instance is used to describe an Enum data
type, and can be passed to the the `Dataset.createVariable` method of
a `Dataset` or `Group` instance. See
`EnumType.__init__` for more details.

The instance variables `dtype`, `name` and `enum_dict` should not be modified by
the user.
    """
    cdef public nc_type _nc_type
    cdef public dtype, name, enum_dict
    def __init__(self, grp, object dt, object dtype_name, object enum_dict, **kwargs):
        """
        **`__init__(group, datatype, datatype_name, enum_dict)`**

        EnumType constructor.

        **`group`**: `Group` instance to associate with the VLEN datatype.

        **`datatype`**: An numpy integer dtype object describing the base type
        for the Enum.

        **`datatype_name`**: a Python string containing a description of the
        Enum data type.

        **`enum_dict`**: a Python dictionary containing the Enum field/value
        pairs.

        ***`Note`***: `EnumType` instances should be created using the
        `Dataset.createEnumType` method of a `Dataset` or
        `Group` instance, not using this class directly.
        """
        cdef nc_type xtype
        if 'typeid' in kwargs:
            xtype = kwargs['typeid']
        else:
            xtype, dt = _def_enum(grp, dt, dtype_name, enum_dict)
        self._nc_type = xtype
        self.dtype = dt
        self.name = dtype_name
        self.enum_dict = enum_dict

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return "%r: name = '%s', numpy dtype = %s, fields/values =%s" %\
            (type(self), self.name, self.dtype, self.enum_dict)

    def __reduce__(self):
        # raise error is user tries to pickle a EnumType object.
        raise NotImplementedError('EnumType is not picklable')

cdef _def_enum(grp, object dt, object dtype_name, object enum_dict):
    # private function used to construct a netCDF Enum data type
    # from a numpy dtype object or python str object by EnumType.__init__.
    cdef nc_type xtype, xtype_tmp
    cdef int ierr, grpid
    cdef char *namstring
    cdef ndarray value_arr
    bytestr = _strencode(dtype_name)
    namstring = bytestr
    grpid = grp._grpid
    dt = numpy.dtype(dt) # convert to numpy datatype.
    if dt.str[1:] in _intnptonctype.keys():
        # find netCDF primitive data type corresponding to
        # specified numpy data type.
        xtype_tmp = _intnptonctype[dt.str[1:]]
        with nogil:
            ierr = nc_def_enum(grpid, xtype_tmp, namstring, &xtype)
        _ensure_nc_success(ierr)
    else:
        msg="unsupported datatype specified for ENUM (must be integer)"
        raise KeyError(msg)
    # insert named members into enum type.
    for field in enum_dict:
        value_arr = numpy.array(enum_dict[field],dt)
        bytestr = _strencode(field)
        namstring = bytestr
        with nogil:
            ierr = nc_insert_enum(grpid, xtype, namstring,
                   PyArray_DATA(value_arr))
        _ensure_nc_success(ierr)
    return xtype, dt

cdef _read_enum(group, nc_type xtype, endian=None):
    # read a Enum data type id from an existing file,
    # construct a corresponding numpy dtype instance,
    # then use that to create a EnumType instance.
    # called by _get_types, _get_vars.
    cdef int ierr, grpid, nmem
    cdef ndarray enum_val
    cdef nc_type base_xtype
    cdef char enum_namstring[NC_MAX_NAME+1]
    cdef size_t nmembers
    grpid = group._grpid
    # get name, datatype, and number of members.
    with nogil:
        ierr = nc_inq_enum(grpid, xtype, enum_namstring, &base_xtype, NULL,\
                &nmembers)
    _ensure_nc_success(ierr)
    enum_name = enum_namstring.decode('utf-8')
    try:
        datatype = _nctonptype[base_xtype]
        if endian is not None: datatype = endian + datatype
        dt = numpy.dtype(datatype) # see if it is a primitive type
    except KeyError:
        raise KeyError("unsupported component type for ENUM")
    # loop over members, build dict.
    enum_dict = {}
    enum_val = numpy.empty(1,dt)
    for nmem from 0 <= nmem < nmembers:
        with nogil:
            ierr = nc_inq_enum_member(grpid, xtype, nmem, \
                                      enum_namstring,PyArray_DATA(enum_val))
        _ensure_nc_success(ierr)
        name = enum_namstring.decode('utf-8')
        enum_dict[name] = enum_val.item()
    return EnumType(group, dt, enum_name, enum_dict, typeid=xtype)

cdef _strencode(pystr,encoding=None):
    # encode a string into bytes.  If already bytes, do nothing.
    # uses 'utf-8' for default encoding.
    if encoding is None:
        encoding = 'utf-8'
    try:
        return pystr.encode(encoding)
    except (AttributeError, UnicodeDecodeError):
        return pystr # already bytes or unicode?

def _to_ascii(bytestr):
    # encode a byte string to an ascii encoded string.
    return str(bytestr,encoding='ascii')

def stringtoarr(string,NUMCHARS,dtype='S'):
    """
**`stringtoarr(a, NUMCHARS,dtype='S')`**

convert a string to a character array of length `NUMCHARS`

**`a`**:  Input python string.

**`NUMCHARS`**:  number of characters used to represent string
(if len(a) < `NUMCHARS`, it will be padded on the right with blanks).

**`dtype`**:  type of numpy array to return.  Default is `'S'`, which
means an array of dtype `'S1'` will be returned.  If dtype=`'U'`, a
unicode array (dtype = `'U1'`) will be returned.

returns a rank 1 numpy character array of length NUMCHARS with datatype `'S1'`
(default) or `'U1'` (if dtype=`'U'`)"""
    if dtype not in ["S","U"]:
        raise ValueError("dtype must string or unicode ('S' or 'U')")
    arr = numpy.zeros(NUMCHARS,dtype+'1')
    arr[0:len(string)] = tuple(string)
    return arr

def stringtochar(a,encoding='utf-8'):
    """
**`stringtochar(a,encoding='utf-8')`**

convert a string array to a character array with one extra dimension

**`a`**:  Input numpy string array with numpy datatype `'SN'` or `'UN'`, where N
is the number of characters in each string.  Will be converted to
an array of characters (datatype `'S1'` or `'U1'`) of shape `a.shape + (N,)`.

optional kwarg `encoding` can be used to specify character encoding (default
`utf-8`). If `encoding` is 'none' or 'bytes', a `numpy.string_` the input array
is treated a raw byte strings (`numpy.string_`).

returns a numpy character array with datatype `'S1'` or `'U1'`
and shape `a.shape + (N,)`, where N is the length of each string in a."""
    dtype = a.dtype.kind
    if dtype not in ["S","U"]:
        raise ValueError("type must string or unicode ('S' or 'U')")
    if encoding in ['none','None','bytes']:
        b = numpy.array(tuple(a.tobytes()),'S1')
    else:
        b = numpy.array(tuple(a.tobytes().decode(encoding)),dtype+'1')
    b.shape = a.shape + (a.itemsize,)
    return b

def chartostring(b,encoding='utf-8'):
    """
**`chartostring(b,encoding='utf-8')`**

convert a character array to a string array with one less dimension.

**`b`**:  Input character array (numpy datatype `'S1'` or `'U1'`).
Will be converted to a array of strings, where each string has a fixed
length of `b.shape[-1]` characters.

optional kwarg `encoding` can be used to specify character encoding (default
`utf-8`). If `encoding` is 'none' or 'bytes', a `numpy.string_` btye array is
returned.

returns a numpy string array with datatype `'UN'` (or `'SN'`) and shape
`b.shape[:-1]` where where `N=b.shape[-1]`."""
    dtype = b.dtype.kind
    if dtype not in ["S","U"]:
        raise ValueError("type must be string or unicode ('S' or 'U')")
    if encoding in ['none','None','bytes']:
        bs = b.tobytes()
    else:
        bs = b.tobytes().decode(encoding)
    slen = int(b.shape[-1])
    if encoding in ['none','None','bytes']:
        a = numpy.array([bs[n1:n1+slen] for n1 in range(0,len(bs),slen)],'S'+repr(slen))
    else:
        a = numpy.array([bs[n1:n1+slen] for n1 in range(0,len(bs),slen)],'U'+repr(slen))
    a.shape = b.shape[:-1]
    return a

class MFDataset(Dataset):
    """
Class for reading multi-file netCDF Datasets, making variables
spanning multiple files appear as if they were in one file.
Datasets must be in `NETCDF4_CLASSIC, NETCDF3_CLASSIC, NETCDF3_64BIT_OFFSET
or NETCDF3_64BIT_DATA` format (`NETCDF4` Datasets won't work).

Adapted from [pycdf](http://pysclint.sourceforge.net/pycdf) by Andre Gosselin.

Example usage (See `MFDataset.__init__` for more details):

```python
>>> import numpy as np
>>> # create a series of netCDF files with a variable sharing
>>> # the same unlimited dimension.
>>> for nf in range(10):
...     with Dataset("mftest%s.nc" % nf, "w", format='NETCDF4_CLASSIC') as f:
...         f.createDimension("x",None)
...         x = f.createVariable("x","i",("x",))
...         x[0:10] = np.arange(nf*10,10*(nf+1))
>>> # now read all those files in at once, in one Dataset.
>>> f = MFDataset("mftest*nc")
>>> print(f.variables["x"][:])
[ 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47
 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71
 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95
 96 97 98 99]
```
    """

    def __init__(self, files, check=False, aggdim=None, exclude=[],
            master_file=None):
        """
        **`__init__(self, files, check=False, aggdim=None, exclude=[],
        master_file=None)`**

        Open a Dataset spanning multiple files, making it look as if it was a
        single file. Variables in the list of files that share the same
        dimension (specified with the keyword `aggdim`) are aggregated. If
        `aggdim` is not specified, the unlimited is aggregated.  Currently,
        `aggdim` must be the leftmost (slowest varying) dimension of each
        of the variables to be aggregated.

        **`files`**: either a sequence of netCDF files or a string with a
        wildcard (converted to a sorted list of files using glob)  If
        the `master_file` kwarg is not specified, the first file
        in the list will become the "master" file, defining all the
        variables with an aggregation dimension which may span
        subsequent files. Attribute access returns attributes only from "master"
        file. The files are always opened in read-only mode.

        **`check`**: True if you want to do consistency checking to ensure the
        correct variables structure for all of the netcdf files.  Checking makes
        the initialization of the MFDataset instance much slower. Default is
        False.

        **`aggdim`**: The name of the dimension to aggregate over (must
        be the leftmost dimension of each of the variables to be aggregated).
        If None (default), aggregate over the unlimited dimension.

        **`exclude`**: A list of variable names to exclude from aggregation.
        Default is an empty list.

        **`master_file`**: file to use as "master file", defining all the
        variables with an aggregation dimension and all global attributes.
       """

        # Open the master file in the base class, so that the CDFMF instance
        # can be used like a CDF instance.
        if isinstance(files, str):
            if files.startswith('http'):
                msg='cannot using file globbing for remote (OPeNDAP) datasets'
                raise ValueError(msg)
            else:
                files = sorted(glob(files))

        if not files:
           msg='no files specified (file list is empty)'
           raise OSError(msg)

        if master_file is not None:
            if master_file not in files:
                raise ValueError('master_file not in files list')
            else:
                master = master_file
        else:
            master = files[0]

        # Open the master again, this time as a classic CDF instance. This will avoid
        # calling methods of the CDFMF subclass when querying the master file.
        cdfm = Dataset(master)
        # copy attributes from master.
        for name, value in cdfm.__dict__.items():
            self.__dict__[name] = value

        # Make sure the master defines a dim with name aggdim,
        # or an unlimited dimension.
        aggDimId = None
        for dimname,dim in cdfm.dimensions.items():
            if aggdim is None:
                if dim.isunlimited():
                    aggDimId = dim
                    aggDimName = dimname
            else:
                if dimname == aggdim:
                    aggDimId = dim
                    aggDimName = dimname
        if aggDimId is None:
            raise OSError("master dataset %s does not have a aggregation dimension" % master)

        # Get info on all aggregation variables defined in the master.
        # Make sure the master defines at least one aggregation variable.
        masterRecVar = {}
        for vName,v in cdfm.variables.items():
            # skip variables specified in exclude list.
            if vName in exclude: continue
            dims = v.dimensions
            shape = v.shape
            dtype = v.dtype
            # Be careful: we may deal with a scalar (dimensionless) variable.
            # Unlimited dimension always occupies index 0.
            if (len(dims) > 0 and aggDimName == dims[0]):
                masterRecVar[vName] = (dims, shape, dtype)
        if len(masterRecVar) == 0:
            raise OSError("master dataset %s does not have any variables to aggregate" % master)

        # Create the following:
        #   cdf       list of Dataset instances
        #   cdfVLen   list unlimited dimension lengths in each CDF instance
        #   cdfRecVar dictionary indexed by the aggregation var names; each key holds
        #             a list of the corresponding Variable instance, one for each
        #             cdf file of the file set
        cdf = []
        self._cdf = cdf        # Store this now, because dim() method needs it
        cdfVLen = []
        cdfRecVar = {}

        # Open each remaining file in read-only mode.
        # Make sure each file defines the same aggregation variables as the master
        # and that the variables are defined in the same way (name, shape and type)
        for f in files:
            if f == master:
                part = cdfm
            else:
                part = Dataset(f)
            if cdfRecVar == {}:
                empty_cdfRecVar = True
            else:
                empty_cdfRecVar = False
            varInfo = part.variables
            for v in masterRecVar.keys():
                if check:
                    # Make sure master rec var is also defined here.
                    if v not in varInfo.keys():
                        raise OSError("aggregation variable %s not defined in %s" % (v, f))

                    #if not vInst.dimensions[0] != aggDimName:

                    masterDims, masterShape, masterType = masterRecVar[v][:3]
                    extDims = varInfo[v].dimensions
                    extShape = varInfo[v].shape
                    extType = varInfo[v].dtype
                    # Check that dimension names are identical.
                    if masterDims != extDims:
                        raise OSError("variable %s : dimensions mismatch between "
                                       "master %s (%s) and extension %s (%s)" %
                                       (v, master, masterDims, f, extDims))

                    # Check that the ranks are identical, and the dimension lengths are
                    # identical (except for that of the unlimited dimension, which of
                    # course may vary.
                    if len(masterShape) != len(extShape):
                        raise OSError("variable %s : rank mismatch between "
                                       "master %s (%s) and extension %s (%s)" %
                                       (v, master, len(masterShape), f, len(extShape)))
                    if masterShape[1:] != extShape[1:]:
                        raise OSError("variable %s : shape mismatch between "
                                       "master %s (%s) and extension %s (%s)" %
                                       (v, master, masterShape, f, extShape))

                    # Check that the data types are identical.
                    if masterType != extType:
                        raise OSError("variable %s : data type mismatch between "
                                       "master %s (%s) and extension %s (%s)" %
                                       (v, master, masterType, f, extType))

                    # Everything ok.
                    if empty_cdfRecVar:
                        cdfRecVar[v] = [part.variables[v]]
                    else:
                        cdfRecVar[v].append(part.variables[v])
                else:
                    # No making sure of anything -- assume this is ok..
                    if empty_cdfRecVar:
                        cdfRecVar[v] = [part.variables[v]]
                    else:
                        cdfRecVar[v].append(part.variables[v])

            cdf.append(part)
            cdfVLen.append(len(part.dimensions[aggDimName]))

        # Attach attributes to the MFDataset instance.
        # A local __setattr__() method is required for them.
        self._files = files            # list of cdf file names in the set
        self._cdfVLen = cdfVLen              # list of unlimited lengths
        self._cdfTLen = sum(cdfVLen) # total length
        self._cdfRecVar = cdfRecVar          # dictionary of Variable instances for all
                                             # the aggregation variables
        self._dims = cdfm.dimensions
        self._grps = cdfm.groups
        for dimname, dim in self._dims.items():
            if dimname == aggDimName:
                self._dims[dimname] = _Dimension(dimname, dim, self._cdfVLen, self._cdfTLen)
        self._vars = cdfm.variables
        for varname,var in self._vars.items():
            if varname in self._cdfRecVar.keys():
                self._vars[varname] = _Variable(self, varname, var, aggDimName)
        self._file_format = []
        self._data_model = []
        self._disk_format = []
        for dset in self._cdf:
            if dset.file_format == 'NETCDF4' or dset.data_model == 'NETCDF4':
                raise ValueError('MFNetCDF4 only works with NETCDF3_* and NETCDF4_CLASSIC formatted files, not NETCDF4')
            self._file_format.append(dset.file_format)
            self._data_model.append(dset.data_model)
            self._disk_format.append(dset.disk_format)
        self._path = '/'

    def __setattr__(self, name, value):
        """override base class attribute creation"""
        self.__dict__[name] = value

    def __getattribute__(self, name):
        if name in ['variables','dimensions','file_format','groups',\
                    'data_model','disk_format','path']:
            if name == 'dimensions': return self._dims
            if name == 'variables': return self._vars
            if name == 'file_format': return self._file_format
            if name == 'data_model': return self._data_model
            if name == 'disk_format': return self._disk_format
            if name == 'path': return self._path
            if name == 'groups': return self._grps
        else:
            return Dataset.__getattribute__(self, name)

    def ncattrs(self):
        """
        **`ncattrs(self)`**

        return the netcdf attribute names from the master file.
        """
        return self._cdf[0].__dict__.keys()

    def close(self):
        """
        **`close(self)`**

        close all the open files.
        """
        for dset in self._cdf:
            dset.close()

    def __repr__(self):
        ncdump = [repr(type(self))]
        dimnames = tuple(str(dimname) for dimname in self.dimensions.keys())
        varnames = tuple(str(varname) for varname in self.variables.keys())
        grpnames = ()
        if self.path == '/':
            ncdump.append('root group (%s data model, file format %s):' %
                    (self.data_model[0], self.disk_format[0]))
        else:
            ncdump.append('group %s:' % self.path)
        for name in self.ncattrs():
            ncdump.append('    %s: %s' % (name, self.__dict__[name]))
        ncdump.append('    dimensions = %s' % str(dimnames))
        ncdump.append('    variables = %s' % str(varnames))
        ncdump.append('    groups = %s' % str(grpnames))
        return '\n'.join(ncdump)

    def __reduce__(self):
        # raise error is user tries to pickle a MFDataset object.
        raise NotImplementedError('MFDataset is not picklable')

class _Dimension:
    def __init__(self, dimname, dim, dimlens, dimtotlen):
        self.dimlens = dimlens
        self.dimtotlen = dimtotlen
        self._name = dimname
    def __len__(self):
        return self.dimtotlen
    def isunlimited(self):
        return True
    def __repr__(self):
        if self.isunlimited():
            return "%r (unlimited): name = '%s', size = %s" %\
                (type(self), self._name, len(self))
        else:
            return "%r: name = '%s', size = %s" %\
                (type(self), self._name, len(self))

class _Variable:
    def __init__(self, dset, varname, var, recdimname):
        self.dimensions = var.dimensions
        self._dset = dset
        self._grp = dset
        self._mastervar = var
        self._recVar = dset._cdfRecVar[varname]
        self._recdimname = recdimname
        self._recLen = dset._cdfVLen
        self.dtype = var.dtype
        self._name = var._name
        # copy attributes from master.
        for name, value in var.__dict__.items():
            self.__dict__[name] = value
    def typecode(self):
        return self.dtype
    def ncattrs(self):
        return self._mastervar.__dict__.keys()
    def __getattr__(self,name):
        if name == 'shape': return self._shape()
        if name == 'ndim': return len(self._shape())
        if name == 'name':  return self._name
        try:
            return self.__dict__[name]
        except:
            raise AttributeError(name)
    def __repr__(self):
        ncdump = [repr(type(self))]
        dimnames = tuple(str(dimname) for dimname in self.dimensions)
        ncdump.append('%s %s%s' % (self.dtype, self._name, dimnames))
        for name in self.ncattrs():
            ncdump.append('    %s: %s' % (name, self.__dict__[name]))
        unlimdims = []
        for dimname in self.dimensions:
            dim = _find_dim(self._grp, dimname)
            if dim.isunlimited():
                unlimdims.append(str(dimname))
        ncdump.append('unlimited dimensions = %r' % (tuple(unlimdims),))
        ncdump.append('current size = %r' % (self.shape,))
        return '\n'.join(ncdump)
    def __len__(self):
        if not self._shape:
            raise TypeError('len() of unsized object')
        else:
            return self._shape()[0]
    def _shape(self):
        recdimlen = len(self._dset.dimensions[self._recdimname])
        return (recdimlen,) + self._mastervar.shape[1:]
    def set_auto_chartostring(self,val):
        for v in self._recVar:
            v.set_auto_chartostring(val)
    def set_auto_maskandscale(self,val):
        for v in self._recVar:
            v.set_auto_maskandscale(val)
    def set_auto_mask(self,val):
        for v in self._recVar:
            v.set_auto_mask(val)
    def set_auto_scale(self,val):
        for v in self._recVar:
            v.set_auto_scale(val)
    def set_always_mask(self,val):
        for v in self._recVar:
            v.set_always_mask(val)
    def __getitem__(self, elem):
        """Get records from a concatenated set of variables."""

        # This special method is used to index the netCDF variable
        # using the "extended slice syntax". The extended slice syntax
        # is a perfect match for the "start", "count" and "stride"
        # arguments to the nc_get_var() function, and is much more easy
        # to use.
        start, count, stride, put_ind =\
        _StartCountStride(elem, self.shape)
        datashape = _out_array_shape(count)
        data = ma.empty(datashape, dtype=self.dtype)

        # Determine which dimensions need to be squeezed
        # (those for which elem is an integer scalar).
        # The convention used is that for those cases,
        # put_ind for this dimension is set to -1 by _StartCountStride.
        squeeze = data.ndim * [slice(None),]
        for i,n in enumerate(put_ind.shape[:-1]):
            if n == 1 and put_ind[...,i].ravel()[0] == -1:
                squeeze[i] = 0

        # Reshape the arrays so we can iterate over them.
        strt = start.reshape((-1, self.ndim or 1))
        cnt = count.reshape((-1, self.ndim or 1))
        strd = stride.reshape((-1, self.ndim or 1))
        put_ind = put_ind.reshape((-1, self.ndim or 1))

        # Fill output array with data chunks.
        # Number of variables making up the MFVariable.Variable.
        nv = len(self._recLen)
        for (start,count,stride,ind) in zip(strt, cnt, strd, put_ind):
            # make sure count=-1 becomes count=1
            count = [abs(cnt) for cnt in count]
            if (numpy.array(stride) < 0).any():
                raise IndexError('negative strides not allowed when slicing MFVariable Variable instance')
            # Start, stop and step along 1st dimension, eg the unlimited
            # dimension.
            sta = start[0]
            step = stride[0]
            stop = sta + count[0] * step

            # Build a list representing the concatenated list of all records in
            # the MFVariable variable set. The list is composed of 2-elem lists
            # each holding:
            #  the record index inside the variables, from 0 to n
            #  the index of the Variable instance to which each record belongs
            idx = []    # list of record indices
            vid = []    # list of Variable indices
            for n in range(nv):
                k = self._recLen[n]     # number of records in this variable
                idx.extend(range(k))
                vid.extend([n] * k)

            # Merge the two lists to get a list of 2-elem lists.
            # Slice this list along the first dimension.
            lst = list(zip(idx, vid)).__getitem__(slice(sta, stop, step))

            # Rebuild the slicing expression for dimensions 1 and ssq.
            newSlice = [slice(None, None, None)]
            for n in range(1, len(start)):   # skip dimension 0
                s = slice(start[n],start[n] + count[n] * stride[n], stride[n])
                newSlice.append(s)

            # Apply the slicing expression to each var in turn, extracting records
            # in a list of arrays.
            lstArr = []
            ismasked = False
            for n in range(nv):
                # Get the list of indices for variable 'n'.
                idx = [i for i,numv in lst if numv == n]
                if idx:
                    # Rebuild slicing expression for dimension 0.
                    newSlice[0] = slice(idx[0], idx[-1] + 1, step)
                    # Extract records from the var, and append them to a list
                    # of arrays.
                    dat = Variable.__getitem__(self._recVar[n],tuple(newSlice))
                    if ma.isMA(dat) and not ismasked:
                        ismasked=True
                        fill_value = dat.fill_value
                    lstArr.append(dat)
            if ismasked:
                lstArr = ma.concatenate(lstArr)
            else:
                lstArr = numpy.concatenate(lstArr)
            if lstArr.dtype != data.dtype: data = data.astype(lstArr.dtype)
            # sometimes there are legitimate singleton dimensions, in which
            # case the array shapes won't conform. If so, a ValueError will
            # result, and no squeeze will be done.
            try:
                data[tuple(ind)] = lstArr.squeeze()
            except ValueError:
                data[tuple(ind)] = lstArr

        # Remove extra singleton dimensions.
        data = data[tuple(squeeze)]

        # if no masked elements, return numpy array.
        if ma.isMA(data) and not data.mask.any():
            data = data.filled()

        return data


class MFTime(_Variable):
    """
Class providing an interface to a MFDataset time Variable by imposing a unique common
time unit and/or calendar to all files.

Example usage (See `MFTime.__init__` for more details):

```python
>>> import numpy as np
>>> f1 = Dataset("mftest_1.nc","w", format="NETCDF4_CLASSIC")
>>> f2 = Dataset("mftest_2.nc","w", format="NETCDF4_CLASSIC")
>>> f1.createDimension("time",None)
>>> f2.createDimension("time",None)
>>> t1 = f1.createVariable("time","i",("time",))
>>> t2 = f2.createVariable("time","i",("time",))
>>> t1.units = "days since 2000-01-01"
>>> t2.units = "days since 2000-02-01"
>>> t1.calendar = "standard"
>>> t2.calendar = "standard"
>>> t1[:] = np.arange(31)
>>> t2[:] = np.arange(30)
>>> f1.close()
>>> f2.close()
>>> # Read the two files in at once, in one Dataset.
>>> f = MFDataset("mftest_*nc")
>>> t = f.variables["time"]
>>> print(t.units)
days since 2000-01-01
>>> print(t[32])  # The value written in the file, inconsistent with the MF time units.
1
>>> T = MFTime(t)
>>> print(T[32])
32
```
    """

    def __init__(self, time, units=None, calendar=None):
        """
        **`__init__(self, time, units=None, calendar=None)`**

        Create a time Variable with units consistent across a multifile
        dataset.

        **`time`**: Time variable from a `MFDataset`.

        **`units`**: Time units, for example, `'days since 1979-01-01'`. If `None`,
        use the units from the master variable.

        **`calendar`**: Calendar overload to use across all files, for example,
        `'standard'` or `'gregorian'`. If `None`, check that the calendar attribute
        is present on each variable and values are unique across files raising a
        `ValueError` otherwise.
        """
        import datetime
        self.__time = time

        # copy attributes from master time variable.
        for name, value in time.__dict__.items():
            self.__dict__[name] = value

        # Make sure calendar attribute present in all files if no default calendar
        # is provided. Also assert this value is the same across files.
        if calendar is None:
            calendars = [None] * len(self._recVar)
            for idx, t in enumerate(self._recVar):
                if not hasattr(t, 'calendar'):
                    msg = 'MFTime requires that the time variable in all files ' \
                          'have a calendar attribute if no default calendar is provided.'
                    raise ValueError(msg)
                else:
                    calendars[idx] = t.calendar
            calendars = set(calendars)
            if len(calendars) > 1:
                msg = 'MFTime requires that the same time calendar is ' \
                      'used by all files if no default calendar is provided.'
                raise ValueError(msg)
            else:
                calendar = list(calendars)[0]

        # Set calendar using the default or the unique calendar value across all files.
        self.calendar = calendar

        # Override units if units is specified.
        self.units = units or time.units

        # Reference date to compute the difference between different time units.
        ref_date = datetime.datetime(1900,1,1)
        ref_num = date2num(ref_date, self.units, self.calendar)

        # Create delta vector: delta = ref_num(ref_date) - num(ref_date)
        # So that ref_num(date) = num(date) + delta
        self.__delta = numpy.empty(len(self), time.dtype)

        i0 = 0; i1 = 0
        for i,v in enumerate(self._recVar):
            n = self._recLen[i] # Length of time vector.
            num = date2num(ref_date, v.units, self.calendar)
            i1 += n
            self.__delta[i0:i1] = ref_num - num
            i0 += n


    def __getitem__(self, elem):
        return self.__time[elem] + self.__delta[elem]
