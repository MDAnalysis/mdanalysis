# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

from . import base
from ..lib import util
from .timestep import Timestep
import MDAnalysis.topology.tpr.utils as tpr_utils
import MDAnalysis.topology.tpr.setting as S


import numpy as np


class TPRReader(base.SingleFrameReaderBase):
    # TODO: reduce duplication with `TPRparser`;
    # we could share some state for the position
    # in the binary file to avoid re-reading topology
    # or perhaps combine the topology and coordinate reading
    # with some inheritance shenanigans?
    format = "TPR"
    units = {"length": "nm", "velocity": "nm/ps"}
    _Timestep = Timestep

    def _read_first_frame(self):
        # Read header/move over topology
        # TODO: reduce duplication with TPRparser perhaps...
        with util.openany(self.filename, mode='rb') as infile:
            tprf = infile.read()
        data = tpr_utils.TPXUnpacker(tprf)
        try:
            th = tpr_utils.read_tpxheader(data)                    # tpxheader
        except (EOFError, ValueError):
            msg = f"{self.filename}: Invalid tpr file or cannot be recognized"
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

        state_ngtc = th.ngtc         # done init_state() in src/gmxlib/tpxio.c
        if th.bBox:
            tpr_utils.extract_box_info(data, th.fver)

        if state_ngtc > 0:
            if th.fver < 69:                      # redundancy due to  different versions
                tpr_utils.ndo_real(data, state_ngtc)
            tpr_utils.ndo_real(data, state_ngtc)        # relevant to Berendsen tcoupl_lambda

        if th.bTop:
            tpr_top = tpr_utils.do_mtop(data, th.fver,
                                        tpr_resid_from_one=True)
        else:
            msg = f"{self.filename}: No topology found in tpr file"
            logger.critical(msg)
            raise IOError(msg)

        if th.bX:
            self.ts._pos = np.asarray(tpr_utils.ndo_rvec(data, th.natoms),
                                      dtype=np.float32)
        if th.bV:
            self.ts._velocities = np.asarray(tpr_utils.ndo_rvec(data, th.natoms),
                                             dtype=np.float32)
            self.ts.has_velocities = True
