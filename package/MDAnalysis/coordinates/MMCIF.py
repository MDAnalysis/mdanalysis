# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
"""
MMCIF structure files in MDAnalysis --- :mod:`MDAnalysis.coordinates.MMCIF`
===========================================================================

.. versionadded:: 2.11.0

MDAnalysis reads coordinates from MMCIF (macromolecular Crystallographic
Information File) files, also known as PDBx/mmCIF format, using the
`gemmi <https://gemmi.readthedocs.io>`_ library as a backend. MMCIF is a
more modern and flexible alternative to the PDB format, capable of storing
detailed structural and experimental data about biological macromolecules.

MMCIF files use a structured, tabular format with key-value pairs to store
both coordinate and atom information. The format supports multiple
models/frames, though this implementation currently only reads the first
model and provides warning messages for multi-model files.

The reader automatically detects if the structure contains placeholder unit
cell information (usually the case for cryoEM structures, where cell
parameters are (1, 1, 1, 90, 90, 90)) and sets dimensions to ``None``
in that case.

Basic usage
-----------

.. code-block:: python

    import MDAnalysis as mda

    u = mda.Universe("structure.cif")

    # or from a compressed file
    u = mda.Universe("structure.cif.gz")

See Also
--------
* `wwPDB MMCIF Resources <http://mmcif.wwpdb.org>`_
* `Gemmi library documentation <https://gemmi.readthedocs.io>`_

Classes and Functions
---------------------

.. autoclass:: MMCIFReader
   :members:
   :inherited-members:

"""

import logging
import warnings
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from ..lib import util
from . import base

if TYPE_CHECKING:
    from gemmi import Model, Structure

try:
    import gemmi

    HAS_GEMMI = True
except ImportError:
    HAS_GEMMI = False

logger = logging.getLogger("MDAnalysis.coordinates.MMCIF")


def _read_gemmi_structure(filename: str | Path) -> "Structure":
    # This function exists because of some lacking methods in the gemmi Python API.
    # Within gemmi in C++, one can call `read_structure` and in-memory, string, and filepath
    # arguments will all be accepted:
    # https://github.com/project-gemmi/gemmi/blob/4416e298f204b7b57bf5b3051d7efd4fe02957cf/include/gemmi/mmread.hpp#L86

    # However, for MDA to similarly accept common input types like streams (open File-like objs and StringIO objs)
    # as well as pathlib.Path() objects, we have to use the Python API methods available currently (as of 0.7.3)
    # with a string as a common target for all input types.
    # For this, we call gemmi.cif.read_string (https://gemmi.readthedocs.io/en/latest/cif.html#reading) to handle CIF
    # strings and gemmi.read_pdb_string to handle PDB strings (no one method can handle both formats currently Py-side)

    # openany() is called instead of passing file paths (when available) differently from streams;
    # even though reading the file into a string is less efficient, this is easier to maintain.

    # If the gemmi Python API is extended, this function can be simplified/removed and replaced with something like
    # gemmi.read_structure
    with util.openany(filename) as f:
        content_as_str = f.read()
    try:
        # String -> Doc -> Block -> Structure
        # making Structure from first Block in Document as is done internally in gemmi:
        # https://github.com/project-gemmi/gemmi/blob/4416e298f204b7b57bf5b3051d7efd4fe02957cf/include/gemmi/mmcif.hpp#L32
        return gemmi.make_structure_from_block(
            gemmi.cif.read_string(content_as_str)[0]
        )
    except ValueError as e:
        try:
            return gemmi.read_pdb_string(content_as_str)
        except ValueError:
            raise e


def _get_coordinates(model: "Model") -> np.ndarray:
    """Get coordinates of all atoms in the `gemmi.Model` object.

    Parameters
    ----------
    model
        input ``gemmi.Model``, e.g. ``gemmi.read_structure('file.cif')[0]``

    Returns
    -------
        np.ndarray, shape [n, 3], where ``n`` is the number of atoms in the structure.
    """
    return np.array(
        [[*at.pos.tolist()] for chain in model for res in chain for at in res]
    )


class MMCIFReader(base.SingleFrameReaderBase):
    """Reads from an MMCIF file using :mod:`gemmi` as a backend.

    Notes
    -----

    If the structure represents an ensemble, only the first structure in the ensemble
    is read here (and a warning is thrown). Also, if the structure has a placeholder "CRYST1"
    record (1, 1, 1, 90, 90, 90), it's set to ``None`` instead.

    .. versionadded:: 2.11.0
    """

    format = ["cif", "cif.gz", "mmcif", "mmcif.gz"]
    units = {"time": None, "length": "Angstrom"}

    def _read_first_frame(self):
        structure = self._get_structure()
        cell_dims = np.array(
            [
                getattr(structure.cell, name)
                for name in ("a", "b", "c", "alpha", "beta", "gamma")
            ]
        )
        if len(structure) > 1:
            wmsg = (
                f"File {self.filename} has {len(structure)} models, "
                "but only the first one will be read"
            )
            warnings.warn(wmsg)
            logger.warning(wmsg)

        model = structure[0]
        coords = _get_coordinates(model)
        self.n_atoms = len(coords)
        self.ts = self._Timestep.from_coordinates(coords, **self._ts_kwargs)
        if np.allclose(cell_dims, np.array([1.0, 1.0, 1.0, 90.0, 90.0, 90.0])):
            wmsg = (
                "1 A^3 CRYST1 record,"
                " this is usually a placeholder."
                " Unit cell dimensions will be set to None."
            )
            warnings.warn(wmsg)
            logger.warning(wmsg)
            self.ts.dimensions = None
        else:
            self.ts.dimensions = cell_dims
        self.ts.frame = 0

    def _get_structure(self):
        return _read_gemmi_structure(self.filename)
