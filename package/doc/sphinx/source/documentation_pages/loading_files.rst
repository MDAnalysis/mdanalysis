##################
Loading your data
##################

   
Basic Usage
===========


MDAnalysis aims to read any and all molecular simulation data,
from whatever files you provide to it.
This is done via the `Universe` object,
which accepts paths to files as its arguments.
Usually two files are required to create a `Universe` object:
A topology file, which provides information about the names and types of atoms,
and a trajectory file, which gives the positions of atoms over time.


This generally corresponds to the file required to set up the simulation acting as the topology,
while the results file of the simulation provides the trajectory.
For example to load results from a CHARMM simulation,
we provide a path to the PSF file to act as a topology,
and a path to the DCD results to act as the trajectory
::

   import MDAnalysis as mda

   u = mda.Universe('adk.psf', 'adk_dims.dcd')


.. note:: If a file which also provides coordinates is used as a topology, no trajectory
	  information is read from this file.  Ie the first frame will come from the trajectory
	  unless the `all_coordinates` keyword is set to ``True``. 

   
Single file Universes
---------------------

Occasionally a file may contain both topology and trajectory information,
in such cases it is sufficient to provide only a single filename to Universe
::

   import MDAnalysis as mda

   u = mda.Universe('myfile.pdb')


Concatenating trajectories
--------------------------
   
It is also possible to read multiple consecutive trajectories,
(for example if a simulation was restarted),
by passing a list of trajectory filenames.
In this example, frames from `traj1.trr` and `traj2.trr` will be concatenated when iterating through the trajectory.
::

   import MDAnalysis as mda

   u = mda.Universe('topol.tpr', ['traj1.trr', 'traj2.trr'])


Supported formats and further details
=====================================




This table lists all currently supported file formats in MDAnalysis,
whether they can act as either Topology or Trajectory files,
as well as links to the relevant documentation pages.
In general MDAnalysis will try and extract all available data from a
given file, for full details of what is extracted from each file format
consult the relevant documentation page.


Generally the format of a file is automatically detected from the extension,
for example a file called `system.xyz` is recognised as an XYZ file.
This can be overriden by supplying the `topology_format` and `format` keyword
arguments to Universe.
A full list of valid values for these keywords are given in the below table.


.. note:: It is possible to pass tarballed/zipped versions of files.  The
	  format detection will work around this.


.. _Supported formats:

.. csv-table:: Table of supported formats
   :header: "Source", "Name", "Format", "Topology", "Trajectory", "I/O"

   ":ref:`Amber <loading_amber>`", ":ref:`PARM parameter/topology <load_amber_top>`", "TOP, PRMTOP, PARM7", "Yes", "Yes", "r"
   "", ":ref:`Ascii trajectory <load_amber_trj>`", "TRJ, MDCRD", "No", "Yes", "r"
   "", ":ref:`Ascii restart <load_amber_restart>`", "INPCRD, RESTRT", "No", "Yes", "r"
   "", ":ref:`NetCFD trajectory <load_amber_ncdf>`", "NCDF, NC", "Minimal", "Yes", "r/w"
   ":ref:`Poisson Boltzmann <load_apbs>`", ":ref:`PQR files <load_pqr>`", "PQR", "Yes", Yes", "r"
   ":ref:`Autodock <load_pdbqt>`", ":ref:`Autodock PDBQT files <load_pdbqt>`", "PDBQT", "Yes", "Yes", "r"
   ":ref:`Charmm <load_charmm>`", ":ref:`PSF files <load_psf>`", "PSF", "Yes", "No", "r"
   "", ":ref:`Binary DCD files <load_dcd>`", "DCD", "Minimal", "Yes", "r/w"
   "", ":ref:`Ascii trajectory <load_crd>`", "CRD", "Minimal", "Yes", "r"
   ":ref:`Desmond MD <load_desmond>`", ":ref:`DMS trajectory <load_dms>`", "DMS", "Yes", "Yes", "r"
   ":ref:`DL Poly <load_dlpoly>`", ":ref:`Ascii History <load_history>`", "HISTORY", "Yes", "Yes", "r"
   "", ":ref:`Ascii config <load_config>`", "CONFIG", "Yes", "Yes", "r"
   ":ref:`GAMESS <load_gamess>`", ":ref:`GAMESS <load_gms>`", "GMS, LOG, OUT", "Yes", "Yes", "r"
   ":ref:`Gromacs <loading_gromacs>`", ":ref:`Gromos <load_gro>`", "GRO", "Yes", "Yes", "r/w"
   "", ":ref:`TPR file <load_tpr>`", "TPR", "Yes", "No", "r"
   "", ":ref:`TRR trajectory <load_trr>`", "TRR", "Minimal", "Yes", "r/w"
   "", ":ref:`XTC trajectory <load_trr>`", "XTC", "Minimal", "Yes", "r/w"
   ":ref:`Hoomd <load_hoomd>`", ":ref:`XML Topology <load_xml>`", "XML", "Yes", "Yes", "r"
   "", ":ref:`Global simulation data? <load_gsd>`", "GSD", "No", Yes", "r"
   ":ref:`IBIsCO and YASP trajectories <load_ibisco>`", ":ref:`Binary trajectory <load_ibisco>`", "TRZ", "Minimal", "Yes", "r/w"
   ":ref:`Lammps <load_lammps>`", ":ref:`Data file <load_data>`", "DATA", "Yes", "Yes", "r"
   "", ":ref:`Binary DCD <load_lammps_dcd>`", "DCD", "Minimal", "Yes", "r/w"
   ":ref:`Protein Databank <load_databank>`", ":ref:`PDB <load_pdb>`", "PDB, ENT, XPDB", "Yes", Yes", "r/w"
   "", ":ref:`Macromolecular transmission format <load_mmtf>`", "MMTF", "Yes", "Yes", "r"
   ":ref:`Tinker <load_tinker>`", ":ref:`Extended XYZ <load_txyz>`", "TXYZ", "Yes", "Yes", "r"
   ":ref:`Tripos <load_tripos>`", ":ref:`MOL2 <load_mol2>`", "MOL2", "Yes", "Yes", "r/w"
   ":ref:`XYZ files <load_xyz>`", ":ref:`Ascii XYZ files <load_xyz>`", "XYZ", "Yes", "Yes", "r/w"

.. toctree::
   :maxdepth: 2
   :hidden:

   ./loading_files/amber
   ./loading_files/apbs
   ./loading_files/autodock
   ./loading_files/charmm
   ./loading_files/desmond
   ./loading_files/dlpoly
   ./loading_files/gamess
   ./loading_files/gromacs
   ./loading_files/hoomd
   ./loading_files/ibisco
   ./loading_files/lammps
   ./loading_files/protein
   ./loading_files/tinker
   ./loading_files/tripos
   ./loading_files/xyz


