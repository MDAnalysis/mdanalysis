.. _load_apbs:

#############################################
Loading Poisson Boltzmann files in MDAnalysis
#############################################

MDAnalysis is able to read and write PQR files produced by
Poisson Boltzmann.

.. _load_pqr:

Loading PQR files
-----------------

Read atoms with charges from a PQR_ file (as written by PDB2PQR_). The
following is adopted from the description of the PQR_ format as used by APBS_:

*MDAnalysis* reads very loosely-formatted PQR files: all fields are
**whitespace-delimited** rather than the strict column formatting mandated
by the PDB_ format. This more liberal formatting allows coordinates
which are larger/smaller than ±999 Å.

MDAnalysis reads data on a per-line basis from PQR files using the following format::

   recordName serial atomName residueName chainID residueNumber X Y Z charge radius

If this fails it is assumed that the *chainID* was omitted and the shorter
format is read::

   recordName serial atomName residueName residueNumber X Y Z charge radius

Anything else will raise a :exc:`ValueError`.

The whitespace is the most important feature of this format: fields
*must* be separated by at least one space or tab character. The fields
are:

*recordName*
    A string which specifies the type of PQR entry and should either be ATOM or
    HETATM.
*serial*
    An integer which provides the atom index (but note that MDAnalysis renumbers
    atoms so one cannot rely on the *serial*)
*atomName*
    A string which provides the atom name.
*residueName*
    A string which provides the residue name.
*chainID*
    An optional string which provides the chain ID of the atom.
*residueNumber*
    An integer which provides the residue index.
*X Y Z*
    Three floats which provide the atomic coordiantes.
*charge*
    A float which provides the atomic charge (in electrons).
*radius*
    A float which provides the atomic radius (in Å).

Clearly, this format can deviate wildly from PDB_ due to the use of whitespaces
rather than specific column widths and alignments. This deviation can be
particularly significant when large coordinate values are used.

Output should look like this (although the only real requirement is
*whitespace* separation between *all* entries). The chainID is optional
and can be omitted::

  ATOM      1  N    MET     1     -11.921   26.307   10.410 -0.3000 1.8500
  ATOM     36  NH1  ARG     2      -6.545   25.499    3.854 -0.8000 1.8500
  ATOM     37 HH11  ARG     2      -6.042   25.480    4.723  0.4600 0.2245


.. Warning:: Fields *must be white-space separated* or data are read
             incorrectly. PDB formatted files are *not* guaranteed to be
             white-space separated so extra care should be taken when quickly
             converting PDB files into PQR files using simple scripts.

For example, PQR files created with PDB2PQR_ and the `--whitespace`
option are guaranteed to conform to the above format::

   pdb2pqr --ff=charmm --whitespace 4ake.pdb 4ake.pqr

.. seealso:: For implementation details, see :mod:`MDAnalysis.coordinates.PQR`
   
.. _PQR:     https://apbs-pdb2pqr.readthedocs.io/en/latest/formats/pqr.html
.. _APBS:    https://apbs-pdb2pqr.readthedocs.io/en/latest/apbs/index.html
.. _PDB2PQR: https://apbs-pdb2pqr.readthedocs.io/en/latest/pdb2pqr/index.html
.. _PDB:     http://www.wwpdb.org/documentation/file-format
