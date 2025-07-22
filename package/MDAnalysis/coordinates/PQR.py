# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#


"""
PQR file format --- :mod:`MDAnalysis.coordinates.PQR`
=====================================================

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


Notes
-----
The PQR format does not provide a means by which to provide box information.
In all cases the `dimensions` attribute will be set to `None`.


.. _PQR:     https://apbs.readthedocs.io/en/latest/formats/pqr.html
.. _APBS:    https://apbs.readthedocs.io/en/latest/
.. _PDB2PQR: https://pdb2pqr.readthedocs.io/en/latest/
.. _PDB:     http://www.wwpdb.org/documentation/file-format
"""
import itertools
import numpy as np
import warnings

from ..lib import util
from . import base


class PQRReader(base.SingleFrameReaderBase):
    """Read a PQR_ file into MDAnalysis.

    .. _PQR:
        https://apbs.readthedocs.io/en/latest/formats/pqr.html

    .. versionchanged:: 0.11.0
       Frames now 0-based instead of 1-based
    """

    format = "PQR"
    units = {"time": None, "length": "Angstrom"}

    # how to slice fields[x:y] to grab coordinates
    _SLICE_INDICES = {
        "ORIGINAL": (-5, -2),
        "NO_CHAINID": (-5, -2),
        "GROMACS": (-6, -3),
    }

    def _read_first_frame(self):
        from ..topology.PQRParser import PQRParser

        flavour = None

        coords = []
        with util.openany(self.filename) as pqrfile:
            for line in pqrfile:
                if line.startswith(("ATOM", "HETATM")):
                    if flavour is None:
                        flavour = PQRParser.guess_flavour(line)
                        x, y = self._SLICE_INDICES[flavour]

                    fields = line.split()
                    # convert all entries at the end once for optimal speed
                    coords.append(fields[x:y])

        self.n_atoms = len(coords)
        self.ts = self._Timestep.from_coordinates(coords, **self._ts_kwargs)
        self.ts.frame = 0  # 0-based frame number
        if self.convert_units:
            # in-place !
            self.convert_pos_from_native(self.ts._pos)

    def Writer(self, filename, **kwargs):
        """Returns a PQRWriter for *filename*.

        Parameters
        ----------
        filename : str
            filename of the output PQR file

        Returns
        -------
        :class:`PQRWriter`
        """
        return PQRWriter(filename, **kwargs)


class PQRWriter(base.WriterBase):
    """Write a single coordinate frame in whitespace-separated PQR format.

    Charges ("Q") are taken from the
    :attr:`MDAnalysis.core.groups.Atom.charge` attribute while
    radii are obtaine from the
    :attr:`MDAnalysis.core.groups.Atom.radius` attribute.

    * If the segid is 'SYSTEM' then it will be set to the empty
      string. Otherwise the first letter will be used as the chain ID.
    * The serial number always starts at 1 and increments sequentially
      for the atoms.

    The output format is similar to ``pdb2pqr --whitespace``.


    .. versionadded:: 0.9.0
    .. versionchanged:: 2.6.0
       Files are now written in `wt` mode, and keep extensions, allowing
       for files to be written under compressed formats
    """

    format = "PQR"
    units = {"time": None, "length": "Angstrom"}

    # serial, atomName, residueName, chainID, residueNumber, XYZ, charge, radius
    fmt_ATOM = (
        "ATOM {serial:6d} {name:<4}  {resname:<3} {chainid:1.1}"
        " {resid:4d}   {pos[0]:-8.3f} {pos[1]:-8.3f}"
        " {pos[2]:-8.3f} {charge:-7.4f} {radius:6.4f}\n"
    )
    fmt_remark = "REMARK   {0} {1}\n"

    def __init__(self, filename, convert_units=True, **kwargs):
        """Set up a PQRWriter with full whitespace separation.

        Parameters
        ----------
        filename : str
             output filename
        convert_units: bool (optional)
            units are converted to the MDAnalysis base format; [``True``]
        remarks : str (optional)
             remark lines (list of strings) or single string to be added to the
             top of the PQR file
        """
        self.filename = util.filename(filename, ext="pqr", keep=True)
        self.convert_units = (
            convert_units  # convert length and time to base units
        )
        self.remarks = kwargs.pop("remarks", "PQR file written by MDAnalysis")

    def write(self, selection, frame=None):
        """Write selection at current trajectory frame to file.

        Parameters
        ----------
        selection : AtomGroup or Universe
            MDAnalysis AtomGroup or Universe
        frame : int (optional)
            optionally move to frame index `frame`; by default, write the
            current frame


        .. versionchanged:: 0.11.0
           Frames now 0-based instead of 1-based

        """
        # write() method that complies with the Trajectory API
        try:
            u = selection.universe
        except AttributeError:
            errmsg = "Input obj is neither an AtomGroup or Universe"
            raise TypeError(errmsg) from None

        if frame is not None:
            u.trajectory[frame]  # advance to frame
        else:
            try:
                frame = u.trajectory.ts.frame
            except AttributeError:
                frame = 0  # should catch cases when we are analyzing a single frame(?)

        atoms = selection.atoms  # make sure to use atoms (Issue 46)
        coordinates = (
            atoms.positions
        )  # can write from selection == Universe (Issue 49)
        if self.convert_units:
            self.convert_pos_to_native(
                coordinates
            )  # inplace because coordinates is already a copy

        # Check atom attributes
        # required:
        # - name
        # - resname
        # - chainid
        # - resid
        # - position
        # - charge
        # - radius
        attrs = {}
        missing_topology = []
        for attr, dflt in (
            ("names", itertools.cycle(("X",))),
            ("resnames", itertools.cycle(("UNK",))),
            ("resids", itertools.cycle((1,))),
            ("charges", itertools.cycle((0.0,))),
            ("radii", itertools.cycle((1.0,))),
        ):
            try:
                attrs[attr] = getattr(atoms, attr)
            except AttributeError:
                attrs[attr] = dflt
                missing_topology.append(attr)

        # chainids require special handling
        # try chainids, then segids
        # if neither, use ' '
        # if 'SYSTEM', use ' '
        try:
            attrs["chainids"] = atoms.chainids
        except AttributeError:
            try:
                attrs["chainids"] = atoms.segids
            except AttributeError:
                pass
        if not "chainids" in attrs or all(attrs["chainids"] == "SYSTEM"):
            attrs["chainids"] = itertools.cycle((" ",))

        if "charges" in missing_topology:
            total_charge = 0.0
        else:
            total_charge = atoms.total_charge()

        if missing_topology:
            warnings.warn(
                "Supplied AtomGroup was missing the following attributes: "
                "{miss}. These will be written with default values. "
                "".format(miss=", ".join(missing_topology))
            )

        with util.openany(self.filename, "wt") as pqrfile:
            # Header / Remarks
            # The *remarknumber* is typically 1 but :program:`pdb2pgr`
            # also uses 6 for the total charge and 5 for warnings.
            for rem in util.asiterable(self.remarks):
                pqrfile.write(self.fmt_remark.format(rem, 1))
            pqrfile.write(
                self.fmt_remark.format(
                    "Input: frame {0} of {1}".format(
                        frame, u.trajectory.filename
                    ),
                    5,
                )
            )
            pqrfile.write(
                self.fmt_remark.format(
                    "total charge: {0:+8.4f} e".format(total_charge), 6
                )
            )

            # Atom descriptions and coords
            for atom_index, (
                pos,
                name,
                resname,
                chainid,
                resid,
                charge,
                radius,
            ) in enumerate(
                zip(
                    coordinates,
                    attrs["names"],
                    attrs["resnames"],
                    attrs["chainids"],
                    attrs["resids"],
                    attrs["charges"],
                    attrs["radii"],
                ),
                start=1,
            ):
                # pad so that only 4-letter atoms are left-aligned
                name = " " + name if len(name) < 4 else name

                pqrfile.write(
                    self.fmt_ATOM.format(
                        serial=atom_index,
                        name=name,
                        resname=resname,
                        chainid=chainid,
                        resid=resid,
                        pos=pos,
                        charge=charge,
                        radius=radius,
                    )
                )
