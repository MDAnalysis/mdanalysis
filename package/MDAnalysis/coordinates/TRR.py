# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.mdanalysis.org
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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
"""
TRR trajectory files --- :mod:`MDAnalysis.coordinates.TRR`
==========================================================

Read and write GROMACS TRR trajectories.

See Also
--------
MDAnalysis.coordinates.XTC: Read and write GROMACS XTC trajectory files.
MDAnalysis.coordinates.XDR: BaseReader/Writer for XDR based formats
"""
from __future__ import absolute_import

from .XDR import XDRBaseReader, XDRBaseWriter
from ..lib.formats.libmdaxdr import TRRFile
from ..lib.mdamath import triclinic_vectors, triclinic_box


class TRRWriter(XDRBaseWriter):
    """Writer for the Gromacs TRR format.

    The Gromacs TRR trajectory format is a lossless format. The TRR format can
    store *velocoties* and *forces* in addition to the coordinates. It is also
    used by other Gromacs tools to store and process other data such as modes
    from a principal component analysis.

    If the data dictionary of a :class:`Timestep` contains the key
    'lambda' the corresponding value will be used as the lambda value
    for written TRR file.  If ``None`` is found the lambda is set to 0.

    """

    format = 'TRR'
    multiframe = True
    units = {'time': 'ps', 'length': 'nm', 'velocity': 'nm/ps',
             'force': 'kJ/(mol*nm)'}
    _file = TRRFile

    def write_next_timestep(self, ts):
        """Write timestep object into trajectory.

        Parameters
        ----------
        ts : :class:`~base.Timestep`

        See Also
        --------
        <FormatWriter>.write(AtomGroup/Universe/TimeStep)
        The normal write() method takes a more general input
        """
        xyz = None
        if ts.has_positions:
            xyz = ts.positions.copy()
            if self._convert_units:
                self.convert_pos_to_native(xyz)

        velo = None
        if ts.has_velocities:
            velo = ts.velocities.copy()
            if self._convert_units:
                self.convert_velocities_to_native(velo)

        forces = None
        if ts.has_forces:
            forces = ts.forces.copy()
            if self._convert_units:
                self.convert_forces_to_native(forces)

        time = ts.time
        step = ts.frame

        if self._convert_units:
            dimensions = self.convert_dimensions_to_unitcell(ts, inplace=False)

        box = triclinic_vectors(dimensions)

        lmbda = 0
        if 'lambda' in ts.data:
            lmbda = ts.data['lambda']

        self._xdr.write(xyz, velo, forces, box, step, time, lmbda,
                        self.n_atoms)


class TRRReader(XDRBaseReader):
    """Reader for the Gromacs TRR format.

    The Gromacs TRR trajectory format is a lossless format. The TRR format can
    store *velocoties* and *forces* in addition to the coordinates. It is also
    used by other Gromacs tools to store and process other data such as modes
    from a principal component analysis.

    The lambda value is written in the data dictionary of the returned
    :class:`Timestep`

    Notes
    -----
    See :ref:`Notes on offsets <offsets-label>` for more information about
    offsets.

    """
    format = 'TRR'
    units = {'time': 'ps', 'length': 'nm', 'velocity': 'nm/ps',
             'force': 'kJ/(mol*nm)'}
    _writer = TRRWriter
    _file = TRRFile

    def _frame_to_ts(self, frame, ts):
        """convert a trr-frame to a mda TimeStep"""
        ts.time = frame.time
        ts.frame = self._frame
        ts.data['step'] = frame.step

        ts.has_positions = frame.hasx
        ts.has_velocities = frame.hasv
        ts.has_forces = frame.hasf
        ts.dimensions = triclinic_box(*frame.box)

        if self.convert_units:
            self.convert_pos_from_native(ts.dimensions[:3])

        if ts.has_positions:
            if self._sub is not None:
                ts.positions = frame.x[self._sub]
            else:
                ts.positions = frame.x
            if self.convert_units:
                self.convert_pos_from_native(ts.positions)

        if ts.has_velocities:
            if self._sub is not None:
                ts.velocities = frame.v[self._sub]
            else:
                ts.velocities = frame.v
            if self.convert_units:
                self.convert_velocities_from_native(ts.velocities)

        if ts.has_forces:
            if self._sub is not None:
                ts.forces = frame.f[self._sub]
            else:
                ts.forces = frame.f
            if self.convert_units:
                self.convert_forces_from_native(ts.forces)

        ts.data['lambda'] = frame.lmbda

        return ts
