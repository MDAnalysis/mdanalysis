.. Contains the formatted docstrings for the serialization of universe located 
.. mainly in 'MDAnalysis/libs/pickle_file_io.py'
.. _serialization:

*********************************************************
Serialization of Universes
*********************************************************

As we approach the exascale barrier, researchers are handling increasingly 
large volumes of molecular dynamics (MD) data. Whilst MDAnalysis is a flexible
and relatively fast framework for complex analysis tasks in MD simulations, 
implementing a parallel computing framework would play a pivotal role in 
accelerating the time to solution for such large datasets.

To achieve a flawless implementation of parallelism, this document illustrates
the basic idea of how current :class:`Universe` is being serialized. And what people
should do to serialize a new reader.

:class:`Universe` is serialized by serializing each necessary compartment---
:attr:`_topology`, :attr:`_trajectory`, :attr:`anchor_name`; a new :class:`Universe` is
reconstructed with these elements. 

To make sure every Trajectory reader can be successfully
serialized, we implement picklable I/O classes (see :ref:`implemented-fileio`).
When the file is pickled, filename and other necessary attributes of the open 
file handle are saved. On unpickling, the file is opened by filename.
This means that for a successful unpickle, the original file still has to
be accessible with its filename. To retain the current frame of the trajectory,
:func:`_read_frame(previous frame)` will be called during unpickling.

Here is an example of how a :class:`Universe` is (un)pickled:

.. code-block:: python

    u = MDAnalysis.Universe(topology, trajectory)
    u.trajectory[5]
    u_pickled = pickle.loads(pickle.dumps(u)
    u_pickled.ts.frame == u.trajectory.ts.frame

.. _simple_parallel_analysis:

Simple analysis case
--------------------

.. code-block:: python

    def cog(u, ag, frame_id):
        u.trajectory[frame_id]
        return ag.center_of_geometry()

    u = MDAnalysis.Universe(GRO, XTC)
    ag = u.atoms[2:5]

    p = multiprocessing.Pool(2)
    result = np.array([p.apply(cog, args=(u, ag, i))
                       for i in range(u.n_frames)])
    p.close()

Note to make sure :class:`Atomgroup` finds its anchored :class:`Universe` after pickling,
the :class:`Universe` has to be pickled first.

.. _how_to_serialize_a_new_reader:

How to serialize a new reader
-----------------------------

File Access
^^^^^^^^^^^
If the new reader uses :func:`util.anyopen()` (e.g. PDB), the reading handler
can be pickled without modification.
If the new reader uses I/O classes from other package (e.g. GSD), and cannot
be pickled natively, create a new picklable class inherited from 
the file class in that package (e.g. GSDPicklable), adding :func:`__getstate__`,
:func:`__setstate__` functions to allow file handler serialization.

To seek or not to seek
^^^^^^^^^^^^^^^^^^^^^^
Some I/O class supports :func:`seek` and :func:`tell` functions to allow the file 
to be pickled with an offset. It is normally not needed for MDAnalysis with
random access. But if error occurs, find a way to make the offset work.

Miscellaneous
^^^^^^^^^^^^^
If pickle still fails due to some unpicklable attributes, try to find a way
to pickle those, or write custom :func:`__getstate__` and :func:`__setstate__`
methods for the reader.

If the new reader is written in Cython, read :class:`lib.formats.libmdaxdr` and
:class:`lib.formats.libdcd` files as references.

Tests
^^^^^
If the test for the new reader uses :class:`BaseReaderTest`, whether
the current timestep information is saved, and whether its relative
position is maintained, i.e. next() read the right next timestep,
are already tested.

If the new reader accesses the file with :func:`util.anyopen`, add necessary
testes inside ``parallelism/test_multiprocessing.py`` for the reader.

If the new reader accessed the file with new picklable I/O class,
add necessary tests inside ``utils/test_pickleio.py`` for I/O class,
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
