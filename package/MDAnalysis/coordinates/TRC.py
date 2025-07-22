# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
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


"""GROMOS11 trajectory reader --- :mod:`MDAnalysis.coordinates.TRC`
====================================================================

Reads coordinates, timesteps and box-sizes from GROMOS11 TRC trajectories.

To load the trajectory into :class:`~MDAnalysis.core.universe.Universe`,
you need to provide topology information using a topology file such as
a PDB::

    import MDAnalysis as mda
    u = mda.Universe("topology.pdb", ["md_1.trc.gz","md_2.trc.gz"],
    continuous=True)

.. Note::
   The GROMOS trajectory format organizes its data in blocks. A block starts
   with a blockname in capital letters (example: POSITIONRED) and ends with
   a line containing only ''END''. Only the TITLE-block at the beginning of
   each file is mandatory, others blocks can be chosen depending on the task.

   The trajectory format as well as all permitted blocks and their data are
   documented in the GROMOS Manual Vol. 4, chapter 2 and 4.
   The manual can be downloaded here:
   https://gromos.net/gromos11_pdf_manuals/vol4.pdf

This reader is designed to read the blocks "TIMESTEP", "POSITIONRED" and
"GENBOX" from the trajectory which covers the standard trajectories of most
simulations . If a trajectory includes other blocks, a warning is served
and the blocks are ignored.

.. Note::
   MDAnalysis requires the blocks to be in the same order for each frame
   and ignores non-supported blocks.



Classes
-------

.. autoclass:: TRCReader
   :members:

"""

import pathlib
import errno
import warnings
import numpy as np

from . import base
from .timestep import Timestep
from ..lib import util
from ..lib.util import cached, store_init_arguments

import logging

logger = logging.getLogger("MDAnalysis.coordinates.GROMOS11")


class TRCReader(base.ReaderBase):
    """Coordinate reader for the GROMOS11 format"""

    format = "TRC"
    units = {"time": "ps", "length": "nm"}
    _Timestep = Timestep

    SUPPORTED_BLOCKS = ["TITLE", "TIMESTEP", "POSITIONRED", "GENBOX"]
    NOT_SUPPORTED_BLOCKS = [
        "POSITION",
        "REFPOSITION",
        "VELOCITY",
        "VELOCITYRED",
        "FREEFORCE",
        "FREEFORCERED",
        "CONSFORCE",
        "CONSFORCERED",
        "SHAKEFAILPOSITION",
        "SHAKEFAILPREVPOSITION",
        "LATTICESHIFTS",
        "COSDISPLACEMENTS",
        "STOCHINT",
        "NHCVARIABLES",
        "ROTTRANSREFPOS",
        "PERTDATA",
        "DISRESEXPAVE",
        "JVALUERESEXPAVE",
        "JVALUERESEPS",
        "JVALUEPERSCALE",
        "ORDERPARAMRESEXPAVE",
        "XRAYRESEXPAVE",
        "LEUSBIAS",
        "LEUSBIASBAS",
        "ENERGY03",
        "VOLUMEPRESSURE03",
        "FREEENERGYDERIVS03",
        "BFACTOR",
        "AEDSS",
    ]

    @store_init_arguments
    def __init__(self, filename, convert_units=True, **kwargs):
        super(TRCReader, self).__init__(filename, **kwargs)

        # GROMOS11 trajectories are usually either *.trc or *.trc.gz.
        # (trj suffix can come up when trajectory obtained from clustering)
        ext = pathlib.Path(self.filename).suffix
        if (ext[1:] == "trc") or (ext[1:] == "trj"):
            self.compressed = None
        else:
            self.compressed = ext[1:]

        self.trcfile = util.anyopen(self.filename)
        self.convert_units = convert_units

        # Read and calculate some information about the trajectory
        self.traj_properties = self._read_traj_properties()

        self._cache = {}
        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)

        self._reopen()
        self.ts.dt = self.traj_properties["dt"]

        self._read_frame(0)

    @property
    @cached("n_atoms")
    def n_atoms(self):
        """The number of atoms in one frame."""
        return self.traj_properties["n_atoms"]

    @property
    @cached("n_frames")
    def n_frames(self):
        """The number of frames in the trajectory."""
        return self.traj_properties["n_frames"]

    def _frame_to_ts(self, frameDat, ts):
        """Convert a frame to a :class:`TimeStep`"""
        ts.frame = self._frame
        ts.time = frameDat["time"]

        ts.data["time"] = frameDat["time"]
        ts.data["step"] = frameDat["step"]

        ts.dimensions = frameDat["dimensions"]
        ts.positions = frameDat["positions"]

        # Convert the units
        if self.convert_units:
            if ts.has_positions:
                self.convert_pos_from_native(ts.positions)
            if ts.dimensions is not None:
                self.convert_pos_from_native(ts.dimensions[:3])

        return ts

    def _read_traj_properties(self):
        """
        * Reads the number of atoms per frame (n_atoms)
        * Reads the number of frames (n_frames)
        * Reads the startposition of the timestep block
          for each frame (l_blockstart_offset)
        """

        traj_properties = {}

        #
        # Check which of the supported blocks comes first in the trajectory
        #
        first_block = None
        with util.anyopen(self.filename) as f:
            for line in iter(f.readline, ""):
                stripped_line = line.strip()
                if stripped_line in TRCReader.SUPPORTED_BLOCKS and (
                    stripped_line != "TITLE"
                ):
                    # Save the name of the first non-title block
                    # in the trajectory file
                    first_block = stripped_line
                    break

        #
        # Calculate meta-data of the trajectory
        #
        n_atoms = 0
        frame_counter = 0
        block_size = {}
        block_size["POSITIONRED"] = None
        inconsistent_size = False

        l_blockstart_offset = []
        l_timestep_timevalues = []

        with util.anyopen(self.filename) as f:
            for line in iter(f.readline, ""):
                stripped_line = line.strip()
                #
                # First block of frame
                #
                if first_block == stripped_line:
                    l_blockstart_offset.append(f.tell() - len(line))
                    frame_counter += 1

                #
                # Timestep-block
                #
                if "TIMESTEP" == stripped_line:
                    l_timestep_timevalues.append(
                        float(f.readline().split()[1])
                    )
                    while stripped_line != "END":
                        stripped_line = f.readline().strip()
                #
                # Coordinates-block
                #
                if "POSITIONRED" == stripped_line:
                    # Only count the number of atoms in the first block
                    # Also stores the size of the POSITIONRED block which
                    # should be consistent in most cases.
                    if n_atoms == 0:
                        start = f.tell() - len(line)
                        while stripped_line != "END":
                            line = f.readline()
                            stripped_line = line.strip()
                            if len(stripped_line.split()) == 3:
                                n_atoms += 1
                        end = f.tell() - len(line)
                        block_size["POSITIONRED"] = end - start

                    elif inconsistent_size:
                        while stripped_line != "END":
                            stripped_line = f.readline().strip()

                    else:
                        # While it is allowed for trailing whitespace to exist
                        # in the trajectories, none of the GROMOS trajectory
                        # writers do actually produce trailing whitespace.
                        # By default we try to assume a consistent size of the
                        # POSITIONRED block. We then try to skip reading subsequent
                        # occurences of this block.
                        # If we end up somewhere unexpected, we will throw
                        # a warning and fall back to iterating over the file.

                        # instead of looping over the file, looking for an END
                        # we can seek to where the end of the block should be.
                        current_pos = f.tell() - len(line)
                        f.seek(
                            f.tell() - len(line) + block_size["POSITIONRED"]
                        )

                        # Check if we are at the correct position
                        # If not, set inconsistent_size to true and seek back
                        # to where we were before
                        if f.readline().strip() != "END":
                            inconsistent_size = True
                            warnmsg = f"Inconsistent POSITIONRED block size in file {f.name}. Falling back to slow reader."
                            warnings.warn(warnmsg, UserWarning)
                            logger.warning(warnmsg)
                            f.seek(current_pos)

        if frame_counter == 0:
            errormsg = (
                "No supported blocks were found within the GROMOS trajectory!"
            )
            logger.error(errormsg)
            raise ValueError(errormsg)

        traj_properties["n_atoms"] = n_atoms
        traj_properties["n_frames"] = frame_counter
        traj_properties["l_blockstart_offset"] = l_blockstart_offset

        if len(l_timestep_timevalues) >= 2:
            traj_properties["dt"] = (
                l_timestep_timevalues[1] - l_timestep_timevalues[0]
            )
        else:
            traj_properties["dt"] = 0
            warnmsg = "The trajectory does not contain TIMESTEP blocks!"
            warnings.warn(warnmsg, UserWarning)
            logger.warning(warnmsg)

        return traj_properties

    def _read_GROMOS11_trajectory(self):
        frameDat = {}
        frameDat["step"] = int(self._frame)
        frameDat["time"] = float(0.0)
        frameDat["positions"] = np.zeros(
            (self.traj_properties["n_atoms"], 3), dtype=np.float64
        )
        frameDat["dimensions"] = None
        self.periodic = False

        # Read the trajectory
        f = self.trcfile
        for line in iter(f.readline, ""):
            stripped_line = line.strip()
            if "TIMESTEP" == stripped_line:
                tmp_step, tmp_time = f.readline().split()
                frameDat["step"] = int(tmp_step)
                frameDat["time"] = float(tmp_time)

            elif "POSITIONRED" == stripped_line:
                i = 0
                buffer = []
                while True:
                    coords_str = f.readline()
                    if "#" in coords_str:
                        continue
                    elif "END" in coords_str:
                        break
                    else:
                        buffer.append(coords_str)
                        i += 1

                data = np.fromstring(
                    "".join(buffer),
                    sep=" ",
                    dtype=np.float64,
                    count=self.n_atoms * 3,
                )

                frameDat["positions"] = data.reshape(-1, 3)

                if i != self.traj_properties["n_atoms"]:
                    errormsg = f"Found {i} atoms in step {self._frame}, but expected {self.traj_properties['n_atoms']}"
                    logger.error(errormsg)
                    raise ValueError(errormsg)

            elif "GENBOX" == stripped_line:
                ntb_setting = int(f.readline())
                if ntb_setting == 0:
                    frameDat["dimensions"] = None
                    self.periodic = False

                elif ntb_setting in [1, 2]:
                    tmp_a, tmp_b, tmp_c = map(float, f.readline().split())
                    tmp_alpha, tmp_beta, tmp_gamma = map(
                        float, f.readline().split()
                    )
                    frameDat["dimensions"] = [
                        tmp_a,
                        tmp_b,
                        tmp_c,
                        tmp_alpha,
                        tmp_beta,
                        tmp_gamma,
                    ]
                    self.periodic = True

                    # gb_line3
                    if (
                        sum(abs(float(v)) for v in f.readline().split())
                        > 1e-10
                    ):
                        errormsg = (
                            "This reader doesnt't support a shifted origin!"
                        )
                        logger.error(errormsg)
                        raise ValueError(errormsg)

                    # gb_line4
                    if (
                        sum(abs(float(v)) for v in f.readline().split())
                        > 1e-10
                    ):
                        errormsg = "This reader doesnt't support yawed, pitched or rolled boxes!"
                        logger.error(errormsg)
                        raise ValueError(errormsg)

                else:
                    errormsg = "This reader does't support truncated-octahedral periodic boundary conditions"
                    logger.error(errormsg)
                    raise ValueError(errormsg)
                break

            elif any(
                non_supp_bn in line
                for non_supp_bn in TRCReader.NOT_SUPPORTED_BLOCKS
            ):
                for non_supp_bn in TRCReader.NOT_SUPPORTED_BLOCKS:
                    if non_supp_bn == stripped_line:
                        warnmsg = non_supp_bn + " block is not supported!"
                        warnings.warn(warnmsg, UserWarning)
                        logger.warning(warnmsg)

        return frameDat

    def _read_frame(self, i):
        """read frame i"""
        self._frame = i - 1

        # Move position in file just (-2 byte) before the start of the block
        self.trcfile.seek(
            self.traj_properties["l_blockstart_offset"][i] - 2, 0
        )

        return self._read_next_timestep()

    def _read_next_timestep(self):
        self._frame += 1
        if self._frame >= self.n_frames:
            errormsg = "Trying to go over trajectory limit"
            logger.error(errormsg)
            raise IOError(errormsg)

        raw_framedata = self._read_GROMOS11_trajectory()
        self._frame_to_ts(raw_framedata, self.ts)

        return self.ts

    def _reopen(self):
        """Close and reopen the trajectory"""
        self.close()
        self.open_trajectory()

    def open_trajectory(self):
        if self.trcfile is not None:
            errormsg = errno.EALREADY, "TRC file already opened", self.filename
            logger.error(errormsg)
            raise IOError(errormsg)

        # Reload trajectory file
        self.trcfile = util.anyopen(self.filename)

        # Reset ts
        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)

        # Set _frame to -1, so next timestep is zero
        self._frame = -1

        return self.trcfile

    def close(self):
        """Close the trc trajectory file if it was open."""
        if self.trcfile is not None:
            self.trcfile.close()
            self.trcfile = None
