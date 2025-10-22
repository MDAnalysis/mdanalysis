# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
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
TPR file format --- :mod:`MDAnalysis.coordinates.TPR`
======================================================

Class for reading positions and velocities from GROMACS TPR files.


Reading TPR files
-----------------
MDAnalysis can read positions and velocities from GROMACS TPR files,
and can also do so using different versions of the topology and
coordinate file for a given system.


For example, both of these are supported modes for reading
in positions and velocities from GROMACS TPR files::

   >>> u = mda.Universe(TPR2020, TPR2024_4)
   >>> u = mda.Universe(TPR2020)


Classes
-------

.. autoclass:: TPRReader
   :members:
   :inherited-members:

"""

from . import base
from ..lib import util
from .timestep import Timestep
import MDAnalysis.topology.tpr.utils as tpr_utils
import MDAnalysis.topology.tpr.setting as S

import logging

logger = logging.getLogger("MDAnalysis.coordinates.TPR")

import numpy as np


class TPRReader(base.SingleFrameReaderBase):
    # TODO: reduce duplication with `TPRparser`;
    # we could share some state for the position
    # in the binary file to avoid re-reading topology
    # or perhaps combine the topology and coordinate reading
    # with some inheritance shenanigans?
    """Class supporting read in of positions and velocities from GROMACS TPR files.

    .. versionadded:: 2.10.0
    """
    format = "TPR"
    units = {"length": "nm", "velocity": "nm/ps"}
    _Timestep = Timestep

    def _read_first_frame(self):
        # Read header/move over topology
        # TODO: reduce duplication with TPRparser perhaps...
        with util.openany(self.filename, mode="rb") as infile:
            tprf = infile.read()
        data = tpr_utils.TPXUnpacker(tprf)
        try:
            th = tpr_utils.read_tpxheader(data)  # tpxheader
        except (EOFError, ValueError):
            msg = f"{self.filename}: Invalid tpr coordinate file or cannot be recognized"
            logger.critical(msg)
            raise IOError(msg)

        self.ts = ts = self._Timestep(th.natoms, **self._ts_kwargs)
        self.n_atoms = th.natoms

        # Starting with gromacs 2020 (tpx version 119), the body of the file
        # is encoded differently. We change the unpacker accordingly.
        if th.fver >= S.tpxv_AddSizeField and th.fgen >= 27:
            actual_body_size = len(data.get_buffer()) - data.get_position()
            if actual_body_size == 4 * th.sizeOfTprBody:
                # See issue #2428.
                msg = (
                    "TPR files produced with beta versions of gromacs 2020 "
                    "are not supported."
                )
                logger.critical(msg)
                raise IOError(msg)
            data = tpr_utils.TPXUnpacker2020.from_unpacker(data)

        state_ngtc = th.ngtc  # done init_state() in src/gmxlib/tpxio.c
        if th.bBox:
            tpr_utils.extract_box_info(data, th.fver)

        if state_ngtc > 0:
            if th.fver < 69:  # redundancy due to  different versions
                tpr_utils.ndo_real(data, state_ngtc)
            tpr_utils.ndo_real(
                data, state_ngtc
            )  # relevant to Berendsen tcoupl_lambda

        tpr_top = tpr_utils.do_mtop(
            data, th.fver, tpr_resid_from_one=True, precision=th.precision
        )

        if th.bX:
            self.ts._pos = np.asarray(
                tpr_utils.ndo_rvec(data, th.natoms), dtype=np.float32
            )
        if th.bV:
            self.ts.velocities = np.asarray(
                tpr_utils.ndo_rvec(data, th.natoms), dtype=np.float32
            )
            self.ts.has_velocities = True
