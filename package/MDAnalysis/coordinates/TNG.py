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
TNG trajectory files --- :mod:`MDAnalysis.coordinates.TNG`
==========================================================


The TNG format :footcite:p:`Lundborg2014` is a format used in `GROMACS`_ for
storage of trajectory and topology information. The TNG format allows a wide
range of compression algorithms and unlike the compressed XTC format can also
store velocities and forces in addition to positions.

The classes in this module are based on the `pytng`_ package for reading TNG
files. The reader is directed to the `pytng documentation`_ for further reading
about how pytng works under the hood. Pytng is an optional dependency and must
be installed to use this Reader. Use of the reader without pytng installed will
raise an `ImportError`. The variable `HAS_PYTNG` indicates whether this module
has pytng availble.

In addition to particle-dependent trajectory information like positions,
forces and velocities, the TNG format can store trajectory metadata and
other arbitrary time dependent data. Additional information can range from
the virial and pressure components to the molecular topology of the system.
This is enabled by a block based system in which binary flags indicate the
presence or absence of various data blocks. The structure of a TNG file is
provided  in the `TNG specification`_. The TNG paper :footcite:p:`Lundborg2014` and
the `pytng documentation`_ are also good resources. The user is encouraged to
read the full list of `TNG blocks`_ to understand the full power of the TNG
format.



Current Limitations
-------------------
Currently there is only a Reader for TNG files and no Writer. This will depend
on upstream changes in pytng. Additionally, it is not currently possible to
read the molecular topology from a TNG file.

There are also some limitations to reading TNG files with pytng.
Currently we do not allow data to be read *off stride*. In essence this
means that all of the critical trajectory data (positions, box, velocities
(if present), forces (if present)) must share the same stride in trajectory
integrator steps. These critical blocks in the TNG file are henceforth called
*special blocks*. Optional blocks (all blocks that are not special blocks)
will not be read if they do not share an integrator step with, or are not
divisible by the shared integrator step of the special blocks.


References
----------

.. footbibliography::


.. Links

.. _GROMACS:
    https://www.gromacs.org/
.. _pytng:
    https://github.com/MDAnalysis/pytng
.. _pytng documentation:
    https://www.mdanalysis.org/pytng/
.. _TNG blocks:
    https://www.mdanalysis.org/pytng/documentation_pages/Blocks.html
.. _TNG specification:
    https://gitlab.com/hugomacdermott/tng/-/blob/master/Trajectoryformatspecification.mk

"""

import warnings
from typing import List, Optional

import numpy as np
from MDAnalysis.coordinates import base
from MDAnalysis.coordinates.timestep import Timestep
from MDAnalysis.lib.mdamath import triclinic_box
from MDAnalysis.lib.util import store_init_arguments

from ..due import Doi, due

try:
    import pytng
except ImportError:
    # Indicates whether pytng is found.
    HAS_PYTNG = False
else:
    # Indicates whether pytng is found.
    HAS_PYTNG = True


class TNGReader(base.ReaderBase):
    r"""Reader for the TNG format

    The TNG format :footcite:p:`Lundborg2014` is a format used in `GROMACS`_ for
    storage of trajectory and topology information. This reader is are based on
    the `pytng`_ package for reading TNG files. The reader is directed to the
    `pytng documentation`_ and `TNG specification`_ for further reading.

    The TNG format allows a wide range of compression
    algorithms and unlike the compressed XTC format can also store velocities
    and forces in addition to positions.

    The contents of the *special blocks* (positions, box, velocities, forces)
    are read into the timestep if present. Additional blocks are read into the
    `ts.data` dictionary if they are available at the current frame.

    .. versionadded:: 2.4.0
    """

    format = "TNG"
    # NOTE: Time units are in seconds unlike other GROMACS formats
    units = {
        "time": "second",
        "length": "nm",
        "velocity": "nm/ps",
        "force": "kJ/(mol*nm)",
    }

    _box_blockname = "TNG_TRAJ_BOX_SHAPE"
    _positions_blockname = "TNG_TRAJ_POSITIONS"
    _velocities_blockname = "TNG_TRAJ_VELOCITIES"
    _forces_blockname = "TNG_TRAJ_FORCES"
    _special_blocks = [
        _box_blockname,
        _positions_blockname,
        _velocities_blockname,
        _forces_blockname,
    ]

    @due.dcite(Doi("10.1002/jcc.23495"), description="The TNG paper", path=__name__)
    @store_init_arguments
    def __init__(self, filename: str, convert_units: bool = True, **kwargs):
        """Initialize a TNG trajectory

        Parameters
        ----------
        filename : str
            filename of the trajectory
        convert_units : bool (optional)
            convert into MDAnalysis units

        """
        if not HAS_PYTNG:
            raise ImportError("TNGReader: To read TNG files please install pytng")

        super(TNGReader, self).__init__(filename, **kwargs)

        self.filename = filename
        self.convert_units = convert_units

        self._file_iterator = pytng.TNGFileIterator(self.filename, "r")
        self.n_atoms = self._file_iterator.n_atoms
        self._n_steps = self._file_iterator.n_steps

        # names of the blocks
        self._block_names = list(self._file_iterator.block_ids.keys())
        # block ids, dict of C long long
        self._block_dictionary = self._file_iterator.block_ids
        self._block_strides = self._file_iterator.block_strides
        self._data_frames = self._file_iterator.n_data_frames
        # init all special blocks to not present
        self._special_block_present = {k: False for k in self._special_blocks}

        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)

        # check for all the special blocks
        self._has_box = self._box_blockname in self._block_names
        if self._has_box:
            self._special_block_present[self._box_blockname] = True
            self.ts.dimensions = np.zeros(6, dtype=np.float32)

        self._has_positions = self._positions_blockname in self._block_names
        if self._has_positions:
            self._special_block_present[self._positions_blockname] = True
            self.ts.positions = self._file_iterator.make_ndarray_for_block_from_name(
                self._positions_blockname
            )

        self._has_velocities = self._velocities_blockname in self._block_names
        if self._has_velocities:
            self._special_block_present[self._velocities_blockname] = True
            self.ts.velocities = self._file_iterator.make_ndarray_for_block_from_name(
                self._velocities_blockname
            )

        self._has_forces = self._forces_blockname in self._block_names
        if self._has_forces:
            self._special_block_present[self._forces_blockname] = True
            self.ts.forces = self._file_iterator.make_ndarray_for_block_from_name(
                self._forces_blockname
            )

        # check for any additional blocks that will be read into ts.data
        self._additional_blocks = [
            block for block in self._block_names if block not in self._special_blocks
        ]
        self._check_strides_and_frames()
        self._frame = 0
        # box needs a temporary place to be stored as ts.dimensions is
        # wrong shape initially
        self._box_temp = self._file_iterator.make_ndarray_for_block_from_name(
            self._box_blockname
        )
        self._frame = -1
        self._read_next_timestep()

    def _check_strides_and_frames(self):
        """
        Check that the strides and frame numbers of the blocks in the TNG
        file match up so that the file can be iterated over retrieving data at
        each integrator step
        """

        strides = []
        n_data_frames = []

        for block in self.special_blocks:
            stride = self._block_strides[block]
            strides.append(stride)
            n_data_frame = self._data_frames[block]
            n_data_frames.append(n_data_frame)

        strides_eq = all(v == strides[0] for v in strides)
        frames_eq = all(v == n_data_frames[0] for v in n_data_frames)

        if (not strides_eq) or (not frames_eq):
            raise IOError(
                "Strides of TNG special blocks not equal, file cannot be read"
            )

        self._global_stride = strides[0]
        # NOTE frame number is 0 indexed so increment
        self._n_frames = n_data_frames[0] + 1

        self._additional_blocks_to_read = []
        for block in self._additional_blocks:
            stride_add = self._block_strides[block]
            if stride_add != self._global_stride:
                if stride_add % self._global_stride:  # pragma: no cover
                    warnings.warn(
                        f"TNG additional block {block} does not match"
                        " strides of other blocks and is not"
                        " divisible by the global stride_length."
                        " It will not be read"
                    )
                else:
                    self._additional_blocks_to_read.append(block)  # pragma: no cover
            else:
                self._additional_blocks_to_read.append(block)

    def close(self):
        """close the reader"""
        self._file_iterator._close()

    @staticmethod
    def parse_n_atoms(filename: str) -> int:
        """parse the number of atoms in the TNG trajectory file

        Parameters
        ----------
        filename : str
            The name of the TNG file

        Returns
        -------
        n_atoms : int
            The number of atoms in the TNG file

        """
        if not HAS_PYTNG:
            raise ImportError("TNGReader: To read TNG files please install pytng")
        with pytng.TNGFileIterator(filename, "r") as tng:
            n_atoms = tng.n_atoms
        return n_atoms

    @property
    def n_frames(self) -> int:
        """number of frames in trajectory

        Returns
        -------
        n_frames : int
            The number of data containing frames in the trajectory
        """
        return self._n_frames

    @property
    def blocks(self) -> List[str]:
        """list of the blocks that are in the file

        Returns
        -------

        blocks : list
            The block names present in the TNG trajectory
        """
        return self._block_names

    @property
    def special_blocks(self) -> List[str]:
        """list of the special blocks that are in the file

        Returns
        -------

        special_blocks : list
            The special block names (positions, box, velocities, forces)
            present in the TNG trajectory
        """
        return [k for k, v in self._special_block_present.items() if v]

    @property
    def additional_blocks(self) -> List[str]:
        """list of the additional (non-special) blocks that are being read from
        the trajectory. This may be exclude some blocks present in the file if
        they do not fall on the same trajectory stride as the positions and
        velocities.

        Returns
        -------

        additional_blocks : list
            The additional block names in the TNG trajectory
        """
        return self._additional_blocks_to_read

    def _reopen(self):
        """reopen the trajectory"""
        self.ts.frame = 0
        self._frame = -1
        self._file_iterator._close()
        self._file_iterator._open(self.filename, "r")

    def _frame_to_step(self, frame: int) -> int:
        """Convert a frame number to an integrator step

        Parameters
        ----------
        frame : int
            The frame number

        Returns
        -------
        integrators_step :  int
            The integrator step of the frame in the TNG file

        """
        return frame * self._global_stride

    def _read_frame(self, i: int) -> Timestep:
        """read frame i

        Parameters
        ----------
        i : int
            The trajectory frame to be read

        Returns
        -------
        ts : Timestep
            Data from frame i encapsulated in an MDA `:class:Timestep`
        """
        self._frame = i - 1
        ts = self._read_next_timestep()
        return ts

    def _read_next_timestep(self, ts: Optional[Timestep] = None) -> Timestep:
        """Read next frame into a timestep

        Parameters
        ----------
        ts : Timestep
            The timestep to read the data into

        Returns
        -------
        ts :  Timestep
            The timestep filled with data from the next step
        """
        if self._frame == self.n_frames - 1:
            raise IOError("trying to go over trajectory limit")
        if ts is None:
            ts = self.ts
        # convert from frames to integrator steps
        step = self._frame_to_step(self._frame + 1)
        iterator_step = self._file_iterator.read_step(step)
        self._frame += 1
        ts = self._frame_to_ts(iterator_step, ts)
        return ts

    def _frame_to_ts(
        self, curr_step: "pytng.TNGCurrentIntegratorStep", ts: Timestep
    ) -> Timestep:
        """convert a TNGCurrentIteratorStep to an MDA Timestep

        Parameters
        ----------
        curr_step : pytng.TNGCurrentIntegratorStep
            The current timestep in the TNG trajectory

        ts : Timestep
            The current timestep in the MDA trajectory

        Returns
        -------
        ts: Timestep
            The current timestep in the MDA trajectory with data from TNG
            trajectory integrator step

        Raises
        ------
        IOError
            The reading of one of the attributes from the TNG file failed
        """

        ts.frame = self._frame
        time = curr_step.get_time()
        if self.convert_units:
            time = self.convert_time_from_native(time)
        ts.time = time
        ts.data["step"] = curr_step.step

        if self._has_box:
            curr_step.get_box(self._box_temp)
            ts.dimensions = triclinic_box(*self._box_temp.reshape(3, 3))
            if not curr_step.read_success:
                raise IOError("Failed to read box from TNG file")
            if self.convert_units:
                self.convert_pos_from_native(ts.dimensions[:3])
        if self._has_positions:
            curr_step.get_positions(ts.positions)
            if not curr_step.read_success:
                raise IOError("Failed to read positions from TNG file")
            if self.convert_units:
                self.convert_pos_from_native(ts.positions)
        if self._has_velocities:
            curr_step.get_velocities(ts.velocities)
            if not curr_step.read_success:
                raise IOError("Failed to read velocities from TNG file")
            if self.convert_units:
                self.convert_velocities_from_native(ts.velocities)
        if self._has_forces:
            curr_step.get_forces(ts.forces)
            if not curr_step.read_success:
                raise IOError("Failed to read forces from TNG file")
            if self.convert_units:
                self.convert_forces_from_native(self.ts.forces)

        for block in self._additional_blocks_to_read:
            add_block_stride = self._block_strides[block]
            # check we are on stride for our block
            if not (add_block_stride % self._global_stride):
                block_data = self._file_iterator.make_ndarray_for_block_from_name(block)
                # additional blocks read into ts.data dictionary
                ts.data[block] = curr_step.get_blockid(
                    self._block_dictionary[block], block_data
                )
                if not curr_step.read_success:
                    raise IOError(
                        f"Failed to read additional block {block} from TNG file"
                    )
        return ts

    def __getstate__(self):
        """Make a dictionary of the class state to pickle Reader instance.

        Must be done manually as pytng uses a non-trivial`__cinit__`.
        """
        state = self.__dict__.copy()
        # cant have PyTNG file iterator in as it is non-pickable
        del state["_file_iterator"]
        return state

    def __setstate__(self, state):
        """Restore class from `state` dictionary in unpickling of Reader
        instance
        """
        self.__dict__ = state
        # reconstruct file iterator
        self._file_iterator = pytng.TNGFileIterator(self.filename, "r")

        # unlike self._read_frame(self._frame),
        # the following lines update the state of the C-level file iterator
        # without updating the ts object.
        # This is necessary to preserve the modification,
        # e.g. changing coordinates, in the ts object.
        # see PR #3722 for more details.
        step = self._frame_to_step(self._frame)
        _ = self._file_iterator.read_step(step)

    def Writer(self):
        """Writer for TNG files

        Raises
        ------
        NotImplementedError
            Currently there is not writer for TNG files pending implementation
            upstream in pytng.
        """
        raise NotImplementedError("There is currently no writer for TNG files")
