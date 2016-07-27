"""MMTF Reader

Insert docs here
"""

import mmtf

from . import base


class MMTFReader(base.SingleFrameReader):
    format = 'MMTF'

    def _read_first_frame(self):
        # TOOD: Check units?
        top = mmtf.parse(self.filename)
        self.n_atoms = top.num_atoms

        self.ts = ts = self._Timestep(self.n_atoms,
                                      **self._ts_kwargs)
        ts._pos[:, 0] = top.x_coord_list.copy()
        ts._pos[:, 1] = top.y_coord_list.copy()
        ts._pos[:, 2] = top.z_coord_list.copy()
        ts.dimensions = top.unit_cell

        ts.data['occupancy'] = top.occupancy_list.copy()

        return ts



