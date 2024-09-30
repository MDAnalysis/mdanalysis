import logging
import warnings

import numpy as np

from . import base

try:
    import gemmi

    HAS_GEMMI = True
except ImportError:
    HAS_GEMMI = False

logger = logging.getLogger("MDAnalysis.coordinates.MMCIF")


def get_coordinates(model: gemmi.Model) -> np.ndarray:
    """Get coordinates of all atoms in the `gemmi.Model` object.

    Parameters
    ----------
    model
        input `gemmi.Model`, e.g. `gemmi.read_structure('file.cif')[0]`

    Returns
    -------
        np.ndarray, shape [n, 3], where `n` is the number of atoms in the structure.
    """
    return np.array(
        [[*at.pos.tolist()] for chain in model for res in chain for at in res]
    )


class MMCIFReader(base.SingleFrameReaderBase):
    """Reads from an MMCIF file"""

    format = "MMCIF"
    units = {"time": None, "length": "Angstrom"}

    def _read_first_frame(self):
        structure = gemmi.read_structure(self.filename)
        if len(structure) > 1:
            warnings.warn(
                f"File {self.filename} has {len(structure)} models, but only the first one will be read"
            )
        if len(structure) > 1:
            warnings.warn(
                "MMCIF model {self.filename} contains {len(model)=} different models, "
                "but only the first one will be used to assign the topology"
            )  # TODO: if the structures represent timestamps, can parse them with :func:`get_coordinates`.

        model = structure[0]
        coords = get_coordinates(model)
        self.n_atoms = len(coords)
        self.ts = self._Timestep.from_coordinates(coords, **self._ts_kwargs)
        self.ts.frame = 0

    def close(self):
        pass
