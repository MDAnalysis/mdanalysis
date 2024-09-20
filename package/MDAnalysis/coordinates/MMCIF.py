import numpy as np
import gemmi
import logging
from . import base

logger = logging.getLogger("MDAnalysis.coordinates.MMCIF")


class MMCIFReader(base.SingleFrameReaderBase):
    """Reads from an MMCIF file"""

    format = "MMCIF"
    units = {"time": None, "length": "Angstrom"}

    def _read_first_frame(self):
        structure = gemmi.read_structure(self.filename)
        coords = np.array(
            [
                [*at.pos.tolist()]
                for model in structure
                for chain in model
                for res in chain
                for at in res
            ]
        )
        self.n_atoms = len(coords)
        self.ts = self._Timestep.from_coordinates(coords, **self._ts_kwargs)
        self.ts.frame = 0

    def Writer(self, filename, n_atoms=None, **kwargs):
        raise NotImplementedError

    def close(self):
        pass


class MMCIFWriter(base.WriterBase):
    pass
