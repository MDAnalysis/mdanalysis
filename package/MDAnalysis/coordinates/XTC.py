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
"""
XTC trajectory files --- :mod:`MDAnalysis.coordinates.XTC`
==========================================================

Read and write GROMACS XTC trajectories.

See Also
--------
MDAnalysis.coordinates.TRR: Read and write GROMACS TRR trajectory files.
MDAnalysis.coordinates.XDR: BaseReader/Writer for XDR based formats
"""

import errno
from . import base
from .XDR import XDRBaseReader, XDRBaseWriter
from ..lib.formats.libmdaxdr import XTCFile
from ..lib.mdamath import triclinic_vectors, triclinic_box


class XTCWriter(XDRBaseWriter):
    """Writer for the Gromacs XTC trajectory format.

    XTC is a compressed trajectory format from Gromacs. The trajectory is saved
    with reduced precision (3 decimal places by default) compared to other
    lossless formarts like TRR and DCD. The main advantage of XTC files is that
    they require significantly less disk space and the loss of precision is
    usually not a problem.
    """
    format = 'XTC'
    multiframe = True
    units = {'time': 'ps', 'length': 'nm'}
    _file = XTCFile

    def __init__(self, filename, n_atoms, convert_units=True,
                 precision=3, **kwargs):
        """
        Parameters
        ----------
        filename : str
            filename of the trajectory
        n_atoms : int
            number of atoms to write
        convert_units : bool (optional)
            convert into MDAnalysis units
        precision : float (optional)
            set precision of saved trjactory to this number of decimal places.
        """
        super(XTCWriter, self).__init__(filename, n_atoms, convert_units,
                                        **kwargs)
        self.precision = precision

    def _write_next_frame(self, ag):
        """Write information associated with ``ag`` at current frame into trajectory

        Parameters
        ----------
        ag : AtomGroup or Universe

        See Also
        --------
        <FormatWriter>.write(AtomGroup/Universe)
        The normal write() method takes a more general input


        .. versionchanged:: 1.0.0
           Added ability to use either AtomGroup or Universe.
        .. versionchanged:: 2.0.0
           Deprecated support for Timestep argument has now been removed.
           Use AtomGroup or Universe as an input instead.
        """
        try:
            # Atomgroup?
            ts = ag.ts
        except AttributeError:
            try:
                # Universe?
                ts = ag.trajectory.ts
            except AttributeError:
                errmsg = "Input obj is neither an AtomGroup or Universe"
                raise TypeError(errmsg) from None

        xyz = ts.positions.copy()
        time = ts.time
        step = ts.frame
        dimensions = ts.dimensions

        if self._convert_units:
            xyz = self.convert_pos_to_native(xyz, inplace=False)
            dimensions = self.convert_dimensions_to_unitcell(ts, inplace=False)

        box = triclinic_vectors(dimensions)
        # libmdaxdr will multiply the coordinates by precision. This means for
        # a precision of 3 decimal places we need to pass 1000.0 to the xdr
        # library.
        precision = 10.0 ** self.precision
        self._xdr.write(xyz, box, step, time, precision)


class XTCReader(XDRBaseReader):
    """Reader for the Gromacs XTC trajectory format.

    XTC is a compressed trajectory format from Gromacs. The trajectory is saved
    with reduced precision (3 decimal places) compared to other lossless
    formarts like TRR and DCD. The main advantage of XTC files is that they
    require significantly less disk space and the loss of precision is usually
    not a problem.

    Notes
    -----
    See :ref:`Notes on offsets <offsets-label>` for more information about
    offsets.



    """
    format = 'XTC'
    units = {'time': 'ps', 'length': 'nm'}
    _writer = XTCWriter
    _file = XTCFile

    def _read_next_timestep(self, ts=None):
        """
        copy next frame into timestep
        
        versionadded:: 2.4.0
            XTCReader implements this method so that it can use
            read_direct_x method of XTCFile to read the data directly
            into the timestep rather than copying it from a temporary array.
        """
        if self._frame == self.n_frames - 1:
            raise IOError(errno.EIO, 'trying to go over trajectory limit')
        if ts is None:
            ts = self.ts
        if ts.has_positions:
            frame = self._xdr.read_direct_x(ts.positions)
        else:
            frame = self._xdr.read()
        self._frame += 1
        self._frame_to_ts(frame, ts)
        return ts

    def _frame_to_ts(self, frame, ts):
        """convert a xtc-frame to a mda TimeStep"""
        ts.frame = self._frame
        ts.time = frame.time
        ts.data['step'] = frame.step
        ts.dimensions = triclinic_box(*frame.box)

        if self._sub is not None:
            ts.positions = frame.x[self._sub]
        else:
            ts.positions = frame.x
        if self.convert_units:
            self.convert_pos_from_native(ts.positions)
            if ts.dimensions is not None:
                self.convert_pos_from_native(ts.dimensions[:3])

        return ts
