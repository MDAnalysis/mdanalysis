.. _load_desmond:

########################################
Loading Desmond MD files with MDAnalysis
########################################

Classes to read a topology from a DESRES_ Molecular Structure file
format (DMS_) coordinate files (as used by the Desmond_ MD package).

.. _DESRES: http://www.deshawresearch.com
.. _Desmond: http://www.deshawresearch.com/resources_desmond.html
.. _DMS: http://www.deshawresearch.com/Desmond_Users_Guide-0.7.pdf


.. _load_dms:

Loading DMS files
-----------------

Using a DMS file for the topology will populate the ``id``, ``atomnum``, ``mass``,
``charge``, ``name``, ``chainID``, ``resid``, ``resname``, ``segid`` and ``type``
attributes.  Atom types are guessed from the atom name.  Additionally bond
information will be read.
For implementation details see :mod:`MDAnalysis.topology.DMSParser`.

When used as a trajectory file, DMS files provide box information, positions
and if available, velocities.  Only a single frame is read.
For implementation details see :mod:`MDAnalysis.coordinates.DMS`.
