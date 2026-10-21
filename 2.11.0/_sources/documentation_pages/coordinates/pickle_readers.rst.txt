.. Contains the formatted docstrings for the serialization of universe located
.. mainly in 'MDAnalysis/libs/pickle_file_io.py'
.. _serialization:

*********************************************************
Serialization of Coordinate Readers
*********************************************************

To achieve a working implementation of parallelism, this document illustrates
the basic idea of how different coordinate readers are being serialized in MDAnalysis,
and what developers should do to serialize a new reader.

To make sure every Trajectory reader can be successfully
serialized, we implement picklable I/O classes (see :ref:`implemented-fileio`).
When the file is pickled, filename and other necessary attributes of the open
file handle are saved. On unpickling, the file is opened by filename.
This means that for a successful unpickle, the original file still has to
be accessible with its filename. To retain the current frame of the trajectory,
:func:`_read_frame(previous frame)` will be called during unpickling.

.. _how_to_serialize_a_new_reader:

How to serialize a new reader
-----------------------------

File Access
^^^^^^^^^^^
If the new reader uses :func:`util.anyopen()`
(e.g. :class:`MDAnalysis.coordinates.PDB.PDBReader`),
the reading handler can be pickled without modification.
If the new reader uses I/O classes from other package
(e.g. :class:`MDAnalysis.coordinates.GSD.GSDReader`),
and cannot be pickled natively, create a new picklable class inherited from
the file class in that package
(e.g. :class:`MDAnalysis.coordinates.GSD.GSDPicklable`),
adding :func:`__getstate__`,
:func:`__setstate__` functions (or :func:`__reduce__` if needed. Consult the
pickle `documentation <https://docs.python.org/3/library/pickle.html>`_ of python)
to allow file handler serialization.

To seek or not to seek
^^^^^^^^^^^^^^^^^^^^^^
Some I/O classes support :func:`seek` and :func:`tell` functions to allow the file
to be pickled with an offset. It is normally not needed for MDAnalysis with
random access. But if error occurs during testing, find a way to make the offset work.
Maybe this I/O class supports frame indexing? Maybe the file handler inside this I/O
class supports offset?

For example, in :class:`MDAnalysis.coordinates.TRZ.TRZReader`,
:func:`_read_frame` is implemented by :func:`_seek` ing the file into
its previous frame and :func:`_read_next_timestep`, so the offset of the file is crucial
for such machinery to work.

Miscellaneous
^^^^^^^^^^^^^
If pickle still fails due to some unpicklable attributes, try to find a way
to pickle those, or write custom :func:`__getstate__` and :func:`__setstate__`
methods for the reader.

If the new reader is written in Cython, read :class:`lib.formats.libmdaxdr` and
:class:`lib.formats.libdcd` files as references.

.. _test_pickle:

Tests
-----
_SingleFrameReader Test
^^^^^^^^^^^^^^^^^^^^^^^
If the new reader is a single-frame reader, the basic test should normally
inherited from :class:`_SingleFrameReader`, where the pickliablity is tested.

BaseReaderTest and MultiframeReaderTest
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If the test for the new reader uses :class:`BaseReaderTest` or
:class:`MultiframeReaderTest`, whether the current timestep information is
saved (the former), whether its relative position is maintained,
i.e. next() reads the right next timestep, and whether its last timestep
can be pickled, are already tested.

File handler Test
^^^^^^^^^^^^^^^^^
If the new reader accesses the file with :func:`util.anyopen`, add necessary
tests inside ``parallelism/test_multiprocessing.py`` for the reader.

If the new reader accessed the file with a new picklable I/O class,
add necessary tests inside ``utils/test_pickleio.py`` for the I/O class,
``parallelism/test_multiprocessing.py`` for the reader.

.. _implemented-fileio:

Currently implemented picklable IO Formats
------------------------------------------

* :class:`MDAnalysis.lib.picklable_file_io.FileIOPicklable`
* :class:`MDAnalysis.lib.picklable_file_io.BufferIOPicklable`
* :class:`MDAnalysis.lib.picklable_file_io.TextIOPicklable`
* :class:`MDAnalysis.lib.picklable_file_io.BZ2Picklable`
* :class:`MDAnalysis.lib.picklable_file_io.GzipPicklable`
* :class:`MDAnalysis.coordinates.GSD.GSDPicklable`
* :class:`MDAnalysis.coordinates.TRJ.NCDFPicklable`
* :class:`MDAnalysis.coordinates.chemfiles.ChemfilesPicklable`
* :class:`MDAnalysis.coordinates.H5MD.H5PYPicklable`
