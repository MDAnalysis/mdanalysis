.. Contains the formatted docstrings for the serialization of universe located 
.. mainly in 'MDAnalysis/libs/pickle_file_io.py'
.. _serialization:

*********************************************************
Serialization of Universes
*********************************************************

.. module:: MDAnalysis.libs.picklable_file_io


.. code-block:: python

    u = MDAnalysis.Universe(topology, trajectory)
    u.trajectory[5]
    u_pickled = pickle.loads(pickle.dumps(u)
    u_pickled.ts.frame == u.trajectory.ts.frame

.. _how_to_serialize_a_new_reader:

How to serialize a new reader
-----------------------------

TODO


.. _implemented-fileio:

Currently implemented picklable FileIO Formats
----------------------------------------------

    :class:`MDAnalysis.lib.picklable_file_io.FileIOPicklable`
    :class:`MDAnalysis.lib.picklable_file_io.BufferIOPicklable`
    :class:`MDAnalysis.lib.picklable_file_io.TextIOPicklable`
    :class:`MDAnalysis.lib.picklable_file_io.BZ2Picklable`
    :class:`MDAnalysis.lib.picklable_file_io.GzipPicklable`
    :class:`MDAnalysis.coordinates.GSD.GSDPicklable`
    :class:`MDAnalysis.coordinates.TRJ.NCDFPicklable`
    :class:`MDAnalysis.coordinates.chemfiles.ChemfilesPicklable`
