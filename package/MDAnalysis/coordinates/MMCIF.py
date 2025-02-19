"""MMCIF structure files in MDAnalysis --- :mod:`MDAnalysis.coordinates.MMCIF`
==========================================================================

MDAnalysis reads coordinates from MMCIF (macromolecular Crystallographic Information File) files
using the gemmi library as a backend. MMCIF is a more modern and flexible alternative to the PDB format,
capable of storing detailed structural and experimental data about biological macromolecules.

MMCIF files use a structured, tabular format with key-value pairs to store both coordinate and atom information.
The format supports multiple models/frames, though this implementation currently only reads the first model
and provides warning messages for multi-model files.

Basic usage
-----------

Reading an MMCIF file is straightforward:

    .. code-block:: python

        import MDAnalysis as mda
        u = mda.Universe("structure.cif")


The reader will automatically detect if the structure contains placeholder unit cell information
(usually it's the case for cryoEM structures, and cell parameters are (1, 1, 1, 90, 90, 90))
and set dimensions to None in that case.

Capabilities
------------

The MMCIF reader implementation uses the gemmi library to parse files and extract coordinates
and unit cell information. Currently only reading capability is supported, with the following
features:

- Single frame/model reading
- Unit cell dimensions detection
- Support for compressed .cif.gz files
- Automatic handling of placeholder unit cells for cryoEM structures

Examples
--------

Basic structure loading::

    .. code-block:: python
        # Load structure from MMCIF
        u = mda.Universe("structure.cif")

        # or from cif.gz file
        u = mda.Universe("structure.cif.gz")


The coordinates are stored in Angstroms.

Classes
-------

.. autoclass:: MMCIFReader
   :members:
   :inherited-members:

See Also
--------
- wwPDB MMCIF Resources: http://mmcif.wwpdb.org
- Gemmi library documentation: https://gemmi.readthedocs.io

.. versionadded:: 2.9.0
"""

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


def get_coordinates(model: "gemmi.Model") -> np.ndarray:
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
    """Reads from an MMCIF file using ``gemmi`` library as a backend.

    Notes
    -----

    If the structure represents an ensemble, only the first structure in the ensemble
    is read here (and a warning is thrown). Also, if the structure has a placeholder "CRYST1"
    record (1, 1, 1, 90, 90, 90), it's set to ``None`` instead.

    .. versionadded:: 2.8.0
    """

    format = ["cif", "cif.gz", "mmcif"]
    units = {"time": None, "length": "Angstrom"}

    def _read_first_frame(self):
        structure = gemmi.read_structure(self.filename)
        cell_dims = np.array(
            [
                getattr(structure.cell, name)
                for name in ("a", "b", "c", "alpha", "beta", "gamma")
            ]
        )
        if len(structure) > 1:
            warnings.warn(  # FIXME: add tests for this
                f"File {self.filename} has {len(structure)} models, but only the first one will be read"
            )
        if len(structure) > 1:
            warnings.warn(  # FIXME: add tests for this
                "MMCIF model {self.filename} contains {len(model)=} different models, "
                "but only the first one will be used to assign the topology"
            )  # TODO: if the structures represent timestamps, can parse them with :func:`get_coordinates`.

        model = structure[0]
        coords = get_coordinates(model)
        self.n_atoms = len(coords)
        self.ts = self._Timestep.from_coordinates(coords, **self._ts_kwargs)
        if np.allclose(cell_dims, np.array([1.0, 1.0, 1.0, 90.0, 90.0, 90.0])):
            warnings.warn(  # FIXME: add tests for this
                "1 A^3 CRYST1 record,"
                " this is usually a placeholder."
                " Unit cell dimensions will be set to None."
            )
            self.ts.dimensions = None  # FIXME: add tests for this
        else:
            self.ts.dimensions = cell_dims
        self.ts.frame = 0
