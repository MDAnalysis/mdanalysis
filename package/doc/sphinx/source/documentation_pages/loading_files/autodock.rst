.. _load_pdbqt:

######################################
Loading Autodock files with MDAnalysis
######################################

MDAnalysis can read PDBQT_ files from AutoDock_ to provide both
topology and coordinate data. It is also possible to write PDBQT files.

When used as a topology file MDAnalysis will populate the ``record_type``, ``id``,
``name``, ``altLoc``, ``resname``, ``segid``,
``resid``, ``icode``, ``occupancy``, ``tempfactor``, ``charge`` and ``type`` attributes.
Atom elements are guessed and masses are guessed and set to 0.0 if not known.
For implementation details see :mod:`MDAnalysis.topology.PDBQTParser`.

When used as a coordinate file the coordiantes are read as a single frame;
multi model PDBQT files are not supported.
For implementation details of reading and writing PDBQT files see
:mod:`MDAnalysis.coordinates.PDBQT`.



.. _PDBQT:
   http://autodock.scripps.edu/faqs-help/faq/what-is-the-format-of-a-pdbqt-file
.. _AutoDock:
   http://autodock.scripps.edu/
