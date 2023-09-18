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


"""GROMOS11 trajectory reader --- :mod:`MDAnalysis.coordinates.TRC`
====================================================================

Reads coordinates, timesteps and box-sizes from GROMOS11 TRC trajectories.

To load the trajectory into :class:`~MDAnalysis.core.universe.Universe`,
you need to provide topology information using a pdb::

    import MDAnalysis as mda
    u = mda.Universe("topology.pdb", ["md_1.trc.gz","md_2.trc.gz"], 
    format="TRC", continuous=True)

.. Note::
   -----
   The reader is capable to read the blocks "TIMESTEP", "POSITIONRED" and 
   "GENBOX" from the trajectory.

Classes
-------

.. autoclass:: TRCReader
   :members:

"""

import os
import errno
import numpy as np

from . import base
from .timestep import Timestep
from ..lib import util
from ..lib.util import cached, store_init_arguments
from ..exceptions import NoDataError
from ..version import __version__

import logging
logger = logging.getLogger("MDAnalysis.coordinates.GROMOS11")


class TRCReader(base.ReaderBase):
    """Coordinate reader for the GROMOS11 format
    """

    format = 'TRC'
    units = {'time': 'ps', 'length': 'nm'}
    _Timestep = Timestep
    
    @store_init_arguments
    def __init__(self, filename, **kwargs):
        super(TRCReader, self).__init__(filename, **kwargs)

        # GROMOS11 trajectories can be either *.trc or *.trc.
        root, ext = os.path.splitext(self.filename)
        self.trcfile = util.anyopen(self.filename)
        self.compression = ext[1:] if ext[1:] != "trj" else None

        # Read and calculate some information about the trajectory
        self.traj_properties = self._read_traj_properties()
                
        self._cache = {}
        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
        
        self._reopen()
        self.ts.dt = self.traj_properties["dt"]
        
        self._read_frame(0)
        
    @property
    @cached('n_atoms')
    def n_atoms(self):
        """The number of atoms in one frame."""
        try:
            return self._read_atom_count()
        except IOError:
            return 0

    def _read_atom_count(self):
        n_atoms = self.traj_properties["n_atoms"] 
        return n_atoms
        
    @property
    @cached('n_frames')
    def n_frames(self):
        """The number of frames in the trajectory."""
        try:
            return self._read_frame_count()
        except IOError:
            return 0

    def _read_frame_count(self):
        n_frames = self.traj_properties["n_frames"] 
        return n_frames
                
    def _frame_to_ts(self, frameDat, ts):
        """Convert a frame to a :class: TimeStep"""
        ts.frame = self._frame
        ts.time = frameDat["time"]
        
        ts.data['time'] = frameDat["time"]
        ts.data['step'] = frameDat["step"]
        
        ts.dimensions = frameDat["dimensions"]
        ts.positions = frameDat["positions"]

        #
        # Convert the units
        #        
        if self.convert_units:
            if ts.has_positions:
                self.convert_pos_from_native(ts.positions)
            if ts.dimensions is not None:
                self.convert_pos_from_native(ts.dimensions[:3])
            if ts.has_velocities:
                self.convert_velocities_from_native(ts.velocities)
      
        return ts
    
    def _read_traj_properties(self):
        """
        * Reads the number of atoms per frame (n_atoms)
        * Reads the number of frames (n_frames) 
        * Reads the startposition of the timestep block 
          for each frame (l_timestep_offset) 
        """
        
        traj_properties = {}
        
        in_positionred_block = False
        in_timestep_block = False
        lastblock_was_timestep = False
        
        atom_counter = 0
        atom_len = 0
        frame_counter = 0  
        frame_len = 0   
 
        l_timestep_offset = []
        l_timestep_timevalues = []
        
        #
        # Loop through the file and save positions of datablocks
        #        
        with util.anyopen(self.filename) as f:
            while (True): 
                line = f.readline() 
                if (line == '') or (line == b''): break  # EOF

                # 
                # Timestep-Block
                # 
                if "TIMESTEP" in line:
                    in_timestep_block = True    
                    lastblock_was_timestep = True
                    l_timestep_offset.append(f.tell())     
                
                elif (lastblock_was_timestep == True):
                    l_timestep_timevalues.append(float(line.split()[1]))
                    lastblock_was_timestep = False       
                
                if ("END" in line) and (in_timestep_block == True):                
                    in_timestep_block = False 

                #
                # Coordinates-Block
                #
                if "POSITIONRED" in line:
                    in_positionred_block = True
                    # coordblock_start = abs_line_counter
                    frame_counter += 1
                    
                if (in_positionred_block == True) and (atom_len == 0):
                    if (len(line.split()) == 3):
                        atom_counter += 1  
                
                if ("END" in line) and (in_positionred_block == True):
                    atom_len = atom_counter
                    in_positionred_block = False
                    
                frame_len = frame_counter

        traj_properties["n_atoms"] = atom_len
        traj_properties["n_frames"] = frame_len
        traj_properties["l_timestep_offset"] = l_timestep_offset
        
        traj_properties["dt"] = (l_timestep_timevalues[1] 
                                 - l_timestep_timevalues[0])

        return traj_properties

    def _read_GROMOS11_trajectory(self, _frame):
    
        frameDat = {}
        f = self.trcfile
        
        l_timestep_offset = self.traj_properties["l_timestep_offset"]

        try:
            
            #
            # Read time
            #           
            while (True):
                line = f.readline()
                if (line == '') or (line == b''): break  # EOF
                if (f.tell() in l_timestep_offset):
                    tmp_step, tmp_time = f.readline().split()
                    frameDat["step"] = int(round(float(tmp_step)))
                    frameDat["time"] = float(tmp_time)
                    
                elif ("POSITIONRED" in line): 
                    tmp_buf = []
                    while (True): 
                        coords_str = f.readline()
                        if '#' in coords_str:
                            continue
                        elif "END" in coords_str:
                            break
                        else:
                            tmp_buf.append(coords_str.split())
                        
                    if (np.array(tmp_buf).shape[0] == self.n_atoms):
                        frameDat["positions"] = np.asarray(tmp_buf, 
                                                           dtype=np.float64)
                    else:
                        raise ValueError("The trajectory contains \
                                         the wrong number of atoms!")                 

                elif ("GENBOX" in line): 
                    ntb_setting = int(f.readline())
                    if (ntb_setting == 0):
                        frameDat["dimensions"] = None
                        self.periodic = False
                                        
                    elif (ntb_setting in [-1, 1]):  
                        tmp_a, tmp_b, tmp_c = f.readline().split()
                        tmp_alpha, tmp_beta, tmp_gamma = f.readline().split()
                        frameDat["dimensions"] = [float(tmp_a), 
                                                  float(tmp_b), 
                                                  float(tmp_c), 
                                                  float(tmp_alpha),
                                                  float(tmp_beta), 
                                                  float(tmp_gamma)]
                        self.periodic = True
                                    
                        line3 = f.readline().split()
                        line4 = f.readline().split()          
                        for v in (line3 + line4):
                            if (float(v) != 0.0):
                                raise NotImplementedError("This reader \
                                                          supports neither \
                                                          triclinic and/or \
                                                          (yawed,pitched, \
                                                          rolled) boxes!")
                    else:
                        raise NotImplementedError("This reader does not \
                                                  support this box type!")                             
                    break
                            
            return frameDat
            
        except (ValueError, IndexError) as err:
            raise EOFError(err) from None
    
    def _read_frame(self, i):
        """read frame i"""
        self._frame = i - 1

        # Move position in file just (-2 byte) before the start of the block 
        self.trcfile.seek(self.traj_properties["l_timestep_offset"][i]-2, 0)

        return self._read_next_timestep()

    def _read_next_timestep(self, ts=None):

        if ts is None:
            ts = self.ts
        
        self._frame += 1
        if (self._frame >= self.n_frames):
            raise EOFError('Trying to go over trajectory limit')

        raw_framedata = self._read_GROMOS11_trajectory(self._frame)        
        self._frame_to_ts(raw_framedata, ts)
        self.ts = ts
        
        return ts
        
    def _reopen(self):
        """Close and reopen the trajectory"""
        self.close()
        self.open_trajectory()

    def open_trajectory(self):
        if self.trcfile is not None:
            raise IOError(
                errno.EALREADY, 'TRC file already opened', self.filename)

        # Reload trajectory file
        self.trcfile = util.anyopen(self.filename)

        # Reset ts
        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
        
        # Set frame to -1, so next timestep is zero
        self._frame = -1 

        return self.trcfile

    def close(self):
        """Close the trc trajectory file if it was open."""
        if self.trcfile is None:
            return
        self.trcfile.close()
        self.trcfile = None
        
