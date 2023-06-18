# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
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


"""CRD structure files in MDAnalysis --- :mod:`MDAnalysis.coordinates.CRD`
===========================================================================

Read and write coordinates in CHARMM CARD coordinate format (suffix
"crd"). The CHARMM "extended format" is handled automatically.

"""
import itertools
import numpy as np
import warnings

from ..exceptions import NoDataError
from ..lib import util
from . import base


class CRDReader(base.SingleFrameReaderBase):
    """CRD reader that implements the standard and extended CRD coordinate formats

    .. versionchanged:: 0.11.0
       Now returns a ValueError instead of FormatError.
       Frames now 0-based instead of 1-based.
    """
    format = 'CRD'
    units = {'time': None, 'length': 'Angstrom'}

    def _read_first_frame(self):
        # EXT:
        #      (i10,2x,a)  natoms,'EXT'
        #      (2I10,2X,A8,2X,A8,3F20.10,2X,A8,2X,A8,F20.10)
        #      iatom,ires,resn,typr,x,y,z,segid,rid,wmain
        # standard:
        #      (i5) natoms
        #      (2I5,1X,A4,1X,A4,3F10.5,1X,A4,1X,A4,F10.5)
        #      iatom,ires,resn,typr,x,y,z,segid,orig_resid,wmain
        coords_list = []
        with util.openany(self.filename) as crdfile:
            extended = False
            natoms = 0
            for linenum, line in enumerate(crdfile):
                if line.strip().startswith('*') or line.strip() == "":
                    continue  # ignore TITLE and empty lines
                fields = line.split()
                if len(fields) <= 2:
                    # should be the natoms line
                    natoms = int(fields[0])
                    extended = (fields[-1] == 'EXT')
                    continue
                # process coordinates
                try:
                    if extended:
                        coords_list.append(np.array(line[45:100].split()[0:3], dtype=float))
                    else:
                        coords_list.append(np.array(line[20:50].split()[0:3], dtype=float))
                except Exception:
                    errmsg = (f"Check CRD format at line {linenum}: "
                              f"{line.rstrip()}")
                    raise ValueError(errmsg) from None

        self.n_atoms = len(coords_list)

        self.ts = self._Timestep.from_coordinates(np.array(coords_list),
                                                  **self._ts_kwargs)
        self.ts.frame = 0  # 0-based frame number
        # if self.convert_units:
        #    self.convert_pos_from_native(self.ts._pos)             # in-place !

        # sanity check
        if self.n_atoms != natoms:
            raise ValueError("Found %d coordinates in %r but the header claims that there "
                              "should be %d coordinates." % (self.n_atoms, self.filename, natoms))

    def Writer(self, filename, **kwargs):
        """Returns a CRDWriter for *filename*.

        Parameters
        ----------
        filename: str
            filename of the output CRD file

        Returns
        -------
        :class:`CRDWriter`

        """
        return CRDWriter(filename, **kwargs)


class CRDWriter(base.WriterBase):
    """CRD writer that implements the CHARMM CRD coordinate format.

    It automatically writes the CHARMM EXT extended format if there
    are more than 99,999 atoms.

    Requires the following attributes to be present:
    - resids
    - resnames
    - names
    - chainIDs
    - tempfactors

    .. versionchanged:: 0.11.0
       Frames now 0-based instead of 1-based
    .. versionchanged:: 2.2.0
       CRD extended format can now be explicitly requested with the
       `extended` keyword
    .. versionchanged:: 2.6.0
       Files are now written in `wt` mode, and keep extensions, allowing
       for files to be written under compressed formats
    """
    format = 'CRD'
    units = {'time': None, 'length': 'Angstrom'}

    fmt = {
        #crdtype = 'extended'
        #fortran_format = '(2I10,2X,A8,2X,A8,3F20.10,2X,A8,2X,A8,F20.10)'
        "ATOM_EXT": ("{serial:10d}{totRes:10d}  {resname:<8.8s}  {name:<8.8s}"
                     "{pos[0]:20.10f}{pos[1]:20.10f}{pos[2]:20.10f}  "
                     "{chainID:<8.8s}  {resSeq:<8d}{tempfactor:20.10f}\n"),
        "NUMATOMS_EXT": "{0:10d} EXT\n",
        #crdtype = 'standard'
        #fortran_format = '(2I5,1X,A4,1X,A4,3F10.5,1X,A4,1X,A4,F10.5)'
        "ATOM": ("{serial:5d}{totRes:5d} {resname:<4.4s} {name:<4.4s}"
                 "{pos[0]:10.5f}{pos[1]:10.5f}{pos[2]:10.5f} "
                 "{chainID:<4.4s} {resSeq:<4d}{tempfactor:10.5f}\n"),
        "TITLE": "* FRAME {frame} FROM {where}\n",
        "NUMATOMS": "{0:5d}\n",
    }

    def __init__(self, filename, **kwargs):
        """
        Parameters
        ----------
        filename : str or :class:`~MDAnalysis.lib.util.NamedStream`
             name of the output file or a stream

        extended : bool (optional)
             By default, noextended CRD format is used [``False``].
             However, extended CRD format can be forced by
             specifying `extended` ``=True``. Note that the extended format
             is *always* used if the number of atoms exceeds 99,999, regardless
             of the setting of `extended`.

             .. versionadded:: 2.2.0
        """

        self.filename = util.filename(filename, ext='crd', keep=True)
        self.crd = None

        # account for explicit crd format, if requested
        self.extended = kwargs.pop("extended", False)

    def write(self, selection, frame=None):
        """Write selection at current trajectory frame to file.

        Parameters
        ----------
        selection : AtomGroup
             group of atoms to be written
        frame : int (optional)
             Move the trajectory to frame `frame`; by default, write
             the current frame.

        """
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
                frame = 0  # should catch cases when we are analyzing a single PDB (?)


        atoms = selection.atoms  # make sure to use atoms (Issue 46)
        coor = atoms.positions  # can write from selection == Universe (Issue 49)

        n_atoms = len(atoms)
        # Detect which format string we're using to output (EXT or not)
        # *len refers to how to truncate various things,
        # depending on output format!
        if self.extended or n_atoms > 99999:
            at_fmt = self.fmt['ATOM_EXT']
            serial_len = 10
            resid_len = 8
            totres_len = 10
        else:
            at_fmt = self.fmt['ATOM']
            serial_len = 5
            resid_len = 4
            totres_len = 5

        # Check for attributes, use defaults for missing ones
        attrs = {}
        missing_topology = []
        for attr, default in (
                ('resnames', itertools.cycle(('UNK',))),
                # Resids *must* be an array because we index it later
                ('resids', np.ones(n_atoms, dtype=int)),
                ('names', itertools.cycle(('X',))),
                ('tempfactors', itertools.cycle((0.0,))),
        ):
            try:
                attrs[attr] = getattr(atoms, attr)
            except (NoDataError, AttributeError):
                attrs[attr] = default
                missing_topology.append(attr)
        # ChainIDs - Try ChainIDs first, fall back to Segids
        try:
            attrs['chainIDs'] = atoms.chainIDs
        except (NoDataError, AttributeError):
            # try looking for segids instead
            try:
                attrs['chainIDs'] = atoms.segids
            except (NoDataError, AttributeError):
                attrs['chainIDs'] = itertools.cycle(('',))
                missing_topology.append(attr)
        if missing_topology:
            warnings.warn(
                "Supplied AtomGroup was missing the following attributes: "
                "{miss}. These will be written with default values. "
                "".format(miss=', '.join(missing_topology)))

        with util.openany(self.filename, 'wt') as crd:
            # Write Title
            crd.write(self.fmt['TITLE'].format(
                frame=frame, where=u.trajectory.filename))
            crd.write("*\n")

            # Write NUMATOMS
            if self.extended or n_atoms > 99999:
                crd.write(self.fmt['NUMATOMS_EXT'].format(n_atoms))
            else:
                crd.write(self.fmt['NUMATOMS'].format(n_atoms))

            # Write all atoms

            current_resid = 1
            resids = attrs['resids']
            for i, pos, resname, name, chainID, resid, tempfactor in zip(
                    range(n_atoms), coor, attrs['resnames'], attrs['names'],
                    attrs['chainIDs'], attrs['resids'], attrs['tempfactors']):
                if not i == 0 and resids[i] != resids[i-1]:
                    current_resid += 1

                # Truncate numbers
                serial = util.ltruncate_int(i + 1, serial_len)
                resid = util.ltruncate_int(resid, resid_len)
                current_resid = util.ltruncate_int(current_resid, totres_len)

                crd.write(at_fmt.format(
                    serial=serial, totRes=current_resid, resname=resname,
                    name=name, pos=pos, chainID=chainID,
                    resSeq=resid, tempfactor=tempfactor))
