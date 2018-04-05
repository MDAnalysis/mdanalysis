.. _load_xyz:

#################################
Loading XYZ files with MDAnalysis
#################################

MDAnalysis is able to read and write the ubiquitous XYZ format to provide
a basic topology and trajectory.  The format follow matches the
`VMD xyzplugin`_ and assumes that columns are whitespace delimited.

.. _`VMD xyzplugin`: http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/xyzplugin.html

When used as a topology the atom names will be read from the file with the
atom types and masses guessed.  Default resids and segids of ``1`` and ``SYSTEM``
will be assigned to all atoms.
For implementation details see :mod:`MDAnalysis.topology.XYZParser`.

When used as a trajectory MDAnalysis will read multiple frames of coordinates.
It is assumed that all frames have (n_atoms + 2) lines.
For implementation details see :mod:`MDAnalysis.coordinates.XYZ`.

.. Note
   The XYZ format provides no box information, any Universe created from
   XYZ files will have this information set to 0.

