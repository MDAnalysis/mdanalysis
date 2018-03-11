.. _loading_amber:

###################################
Loading Amber files with MDAnalysis
###################################

MDAnalysis can read PRMTOP files as topologies and
NCDF and ascii coordinate and restart files as trajectories from
`Amber MD`_ simulations.
Typically using NCDF as the trajectory format will give the best performance
as all other trajectory formats are ascii based.

.. _Amber MD: http://ambermd.org

.. _load_amber_top:

Loading Amber PRMTOP files
--------------------------

MDAnalysis reads `TOP format`_ files with an extension of ``TOP``,
``PRMTOP`` or ``PARM7``.
The below table shows the relationship between Amber flags and the attribute
created in MDAnalysis:
The ``type_indices`` attribute is unique to Amber formats and is
an *integer* representation of the atom type, rather than the
typical string representation found throughout MDAnalysis.

.. table:: Mapping of Amber flags to MDAnalysis attributes

  +-----------------+----------------------+
  | AMBER flag      | MDAnalysis attribute |
  +=================+======================+
  | ATOM_NAME       | names                |
  +-----------------+----------------------+
  | CHARGE          | charges              |
  +-----------------+----------------------+
  | ATOMIC_NUMBER   | elements             |
  +-----------------+----------------------+
  | MASS            | masses               |
  +-----------------+----------------------+
  | ATOM_TYPE_INDEX | type_indices         |
  +-----------------+----------------------+
  | AMBER_ATOM_TYPE | types                |
  +-----------------+----------------------+
  | RESIDUE_LABEL   | resnames             |
  +-----------------+----------------------+
  | RESIDUE_POINTER | residues             |
  +-----------------+----------------------+

.. note::

   The Amber charge is converted to electron charges as used in
   MDAnalysis and other packages. To get back Amber charges, multiply
   by 18.2223.

For implementation details, see
:mod:`MDAnalysis.topology.TOPParser`.

.. _`TOP format`: http://ambermd.org/formats.html#topo.cntrl

.. _load_amber_ncdf:

Loading Amber NCDF files
------------------------

The `AMBER netcdf`_ format make use of NetCDF_ (Network Common Data
Form) format. Such binary trajectories are recognized in MDAnalysis by
the '.ncdf' or '.nc' suffix.

Binary trajectories can also contain velocities and forces, and can record the
exact time
step. In principle, the trajectories can be in different units than the AMBER
defaults of ångström and picoseconds but at the moment MDAnalysis only supports
those and will raise a :exc:`NotImplementedError` if anything else is detected.

For implementation details see :class:`MDAnalysis.coordinates.TRJ.NCDFReader`.

.. _AMBER netcdf: http://ambermd.org/netcdf/nctraj.xhtml
.. _NetCDF: http://www.unidata.ucar.edu/software/netcdf


.. _load_amber_trj:

Loading TRJ files
-----------------

ASCII `AMBER TRJ coordinate files`_
are handled by the :class:`TRJReader`. It is also possible to directly
read *bzip2* or *gzip* compressed files.

AMBER ASCII trajectories are recognised by the suffix '.trj',
'.mdcrd' or '.crdbox (possibly with an additional '.gz' or '.bz2').

.. Note::

   In the AMBER community, these trajectories are often saved with the
   suffix '.crd' but this extension conflicts with the CHARMM CRD
   format and MDAnalysis *will not correctly autodetect AMBER ".crd"
   trajectories*. Instead, explicitly provide the ``format="TRJ"``
   argument to :class:`~MDAnalysis.core.universe.Universe`::

     u = MDAnalysis.Universe("top.prmtop", "traj.crd", format="TRJ")

   In this way, the format is correct set.

.. rubric:: Limitations

* Periodic boxes are only stored as box lengths A, B, C in an AMBER
  trajectory; the reader always assumes that these are orthorhombic
  boxes.

* The trajectory does not contain time information so we simply set
  the time step to 1 ps (or the user could provide it as kwarg *dt*)

* Trajectories with fewer than 4 atoms probably fail to be read (BUG).

* If the trajectory contains exactly *one* atom then it is always
  assumed to be non-periodic (for technical reasons).

* Velocities are currently *not supported* as ASCII trajectories.

For implementation details see
:class:`MDAnalysis.coordinates.TRJ.TRJReader`.

.. _AMBER TRJ coordinate files: http://ambermd.org/formats.html#trajectory

.. _load_amber_restart:

Loading Amber restart files
---------------------------

A single frame of positions can be read from `Amber restart files`_,
which require a file suffix of ``.inpcrd`` or ``.restrt``.
See :mod:`MDAnalysis.coordinates.INPCRD` for implementation details.

.. _Amber restart files: http://ambermd.org/formats.html#restart

