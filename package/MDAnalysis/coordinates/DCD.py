
#
# based on MDAnalysis.coordinates.DCD from
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2011 Naveen Michaud-Agrawal,
#               Elizabeth J. Denning, Oliver Beckstein,
#               and contributors (see website for details)
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
#     N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
#     O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
#     Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
#     doi:10.1002/jcc.21787
#

from __future__ import absolute_import, division, print_function, unicode_literals

import os
import errno
import types
import numpy as np

from ..core import flags
from .. import units as mdaunits  # use mdaunits instead of units to avoid a clash
from ..exceptions import NoDataError
from . import base
from . import core
# dcdtimeseries is implemented with Pyrex - hopefully all dcd reading functionality can move to pyrex
from . import _dcdmodule
from . import dcdtimeseries

class Timestep(base.Timestep):
    #: Indices into :attr:`Timestep._unitcell` (``[A, gamma, B, beta, alpha,
    #: C]``, provided by the :class:`DCDReader` C code) to pull out
    #: ``[A, B, C, alpha, beta, gamma]``.
    _ts_order = [0, 2, 5, 4, 3, 1]

    @property
    def dimensions(self):
        """unitcell dimensions (*A*, *B*, *C*, *alpha*, *beta*, *gamma*)
        lengths *A*, *B*, *C* are in the MDAnalysis length unit (Å), and
        angles are in degrees.
        :attr:`dimensions` is read-only because it transforms the actual format
        of the unitcell (which differs between different trajectory formats) to
        the representation described here, which is used everywhere in
        MDAnalysis.
        The ordering of the angles in the unitcell is the same as in recent
        versions of VMD's DCDplugin_ (2013), namely the `X-PLOR DCD format`_:
        The original unitcell is read as ``[A, gamma, B, beta, alpha, C]`` from
        the DCD file (actually, the direction cosines are stored instead of the
        angles but the underlying C code already does this conversion); if any
        of these values are < 0 or if any of the angles are > 180 degrees then
        it is assumed it is a new-style CHARMM unitcell (at least since c36b2)
        in which box vectors were recorded.
        .. warning:: The DCD format is not well defined. Check your unit cell
           dimensions carefully, especially when using triclinic
           boxes. Different software packages implement different conventions
           and MDAnalysis is currently implementing the newer NAMD/VMD convention
           and tries to guess the new CHARMM one. Old CHARMM trajectories might
           give wrong unitcell values. For more details see `Issue 187`_.
        .. versionchanged:: 0.9.0
           Unitcell is now interpreted in the newer NAMD DCD format as ``[A,
           gamma, B, beta, alpha, C]`` instead of the old MDAnalysis/CHARMM
           ordering ``[A, alpha, B, beta, gamma, C]``. We attempt to detect the
           new CHARMM DCD unitcell format (see `Issue 187`_ for a discussion).
        .. _`X-PLOR DCD format`: http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html
        .. _Issue 187: https://github.com/MDAnalysis/mdanalysis/issues/187
        .. _DCDplugin: http://www.ks.uiuc.edu/Research/vmd/plugins/doxygen/dcdplugin_8c-source.html#l00947
        """

        # Layout of unitcell is [A, alpha, B, beta, gamma, C] --- (originally CHARMM DCD)
        # override for other formats; this strange ordering is kept for historical reasons
        # (the user should not need concern themselves with this)
        ## orig MDAnalysis 0.8.1/dcd.c (~2004)
        ##return np.take(self._unitcell, [0,2,5,1,3,4])

        # MDAnalysis 0.9.0 with recent dcd.c (based on 2013 molfile
        # DCD plugin, which implements the ordering of recent NAMD
        # (>2.5?)). See Issue 187.
        uc = np.take(self._unitcell, self._ts_order)
        # heuristic sanity check: uc = A,B,C,alpha,beta,gamma
        # XXX: should we worry about these comparisons with floats?
        if np.any(uc < 0.) or np.any(uc[3:] > 180.):
            # might be new CHARMM: box matrix vectors
            H = self._unitcell
            e1, e2, e3 = H[[0,1,3]],  H[[1,2,4]], H[[3,4,5]]
            uc = core.triclinic_box(e1, e2, e3)
        return uc

    @dimensions.setter
    def dimensions(self, box):
        """Set unitcell with (*A*, *B*, *C*, *alpha*, *beta*, *gamma*)
        .. versionadded:: 0.9.0
        """
        # note that we can re-use self._ts_order with put!
        np.put(self._unitcell, self._ts_order, box)

class DCDWriter(base.Writer):
    """Writes to a DCD file
    Typical usage::
       with DCDWriter("new.dcd", u.atoms.n_atoms) as w:
           for ts in u.trajectory
               w.write_next_timestep(ts)
    Keywords are available to set some of the low-level attributes of the DCD.
    :Methods:
       ``d = DCDWriter(dcdfilename, n_atoms, start, step, delta, remarks)``
    .. Note::
       The Writer will write the **unit cell information** to the DCD in a
       format compatible with NAMD and older CHARMM versions, namely the unit
       cell lengths in Angstrom and the angle cosines (see
       :class:`Timestep`). Newer versions of CHARMM (at least c36b2) store the
       matrix of the box vectors. Writing this matrix to a DCD is currently not
       supported (although reading is supported with the
       :class:`DCDReader`); instead the angle cosines are written,
       which *might make the DCD file unusable in CHARMM itself*. See
       `Issue 187`_ for further information.
       The writing behavior of the :class:`DCDWriter` is identical to
       that of the DCD molfile plugin of VMD with the exception that
       by default it will use AKMA time units.
    .. _Issue 187: https://github.com/MDAnalysis/mdanalysis/issues/187
    """
    format = 'DCD'
    multiframe = True
    flavor = 'CHARMM'
    units = {'time': 'AKMA', 'length': 'Angstrom'}

    def __init__(self, filename, n_atoms, start=0, step=1,
                 delta=mdaunits.convert(1., 'ps', 'AKMA'), dt=None,
                 remarks="Created by DCDWriter", convert_units=None):
        """Create a new DCDWriter
        :Arguments:
         *filename*
           name of output file
         *n_atoms*
           number of atoms in dcd file
         *start*
           starting timestep
         *step*
           skip between subsequent timesteps (indicate that *step* MD
           integrator steps (!) make up one trajectory frame); default is 1.
         *delta*
           timestep (MD integrator time step (!), in AKMA units); default is
           20.45482949774598 (corresponding to 1 ps).
         *remarks*
           comments to annotate dcd file
         *dt*
           **Override** *step* and *delta* so that the DCD records that *dt* ps
           lie between two frames. (It sets *step* = 1 and *delta* = ``AKMA(dt)``.)
           The default is ``None``, in which case *step* and *delta* are used.
         *convert_units*
           units are converted to the MDAnalysis base format; ``None`` selects
           the value of :data:`MDAnalysis.core.flags` ['convert_lengths'].
           (see :ref:`flags-label`)
       .. Note::
          The keyword arguments set the low-level attributes of the DCD
          according to the CHARMM format. The time between two frames would be
          *delta* * *step* ! For convenience, one can alternatively supply the
          *dt* keyword (see above) to just tell the writer that it should
          record "There are dt ps between each frame".
        """
        if n_atoms == 0:
            raise ValueError("DCDWriter: no atoms in output trajectory")
        elif n_atoms is None:
            # probably called from MDAnalysis.Writer() so need to give user a gentle heads up...
            raise ValueError("DCDWriter: REQUIRES the number of atoms in the 'n_atoms' argument\n" +
                             " " * len("ValueError: ") +
                             "For example: n_atoms=universe.atoms.n_atoms")
        self.filename = filename
        # convert length and time to base units on the fly?
        self.convert_units = flags['convert_lengths'] if convert_units is None \
            else convert_units
        self.n_atoms = n_atoms

        self.frames_written = 0
        self.start = start
        if dt is not None:
            if dt > 0:
                # ignore step and delta
                self.step = 1
                self.delta = mdaunits.convert(dt, 'ps', 'AKMA')
            else:
                raise ValueError("DCDWriter: dt must be > 0, not {0}".format(dt))
        else:
            self.step = step
            self.delta = delta
        self.dcdfile = open(self.filename, 'wb')
        self.remarks = remarks
        self._write_dcd_header(self.n_atoms, self.start, self.step, self.delta, self.remarks)

    def _dcd_header(self):
        """Returns contents of the DCD header C structure::
             typedef struct {
               fio_fd fd;                 // FILE *
               fio_size_t header_size;    // size_t == sizeof(int)
               int natoms;
               int nsets;
               int setsread;
               int istart;
               int nsavc;
               double delta;
               int nfixed;
               int *freeind;
               float *fixedcoords;
               int reverse;
               int charmm;
               int first;
               int with_unitcell;
             } dcdhandle;
        .. deprecated:: 0.7.5
           This function only exists for debugging purposes and might
           be removed without notice. Do not rely on it.
        """
        # was broken (no idea why [orbeckst]), see Issue 27
        # 'PiiiiiidiPPiiii' should be the unpack string according to the struct.
        #    struct.unpack("LLiiiiidiPPiiii",self._dcd_C_str)
        # seems to do the job on Mac OS X 10.6.4 ... but I have no idea why,
        # given that the C code seems to define them as normal integers
        desc = [
            'file_desc', 'header_size', 'natoms', 'nsets', 'setsread', 'istart',
            'nsavc', 'delta', 'nfixed', 'freeind_ptr', 'fixedcoords_ptr',
            'reverse', 'charmm', 'first', 'with_unitcell']
        return dict(zip(desc, struct.unpack("LLiiiiidiPPiiii", self._dcd_C_str)))

    def write_next_timestep(self, ts=None):
        ''' write a new timestep to the dcd file
        *ts* - timestep object containing coordinates to be written to dcd file
        .. versionchanged:: 0.7.5
           Raises :exc:`ValueError` instead of generic :exc:`Exception`
           if wrong number of atoms supplied and :exc:`~MDAnalysis.NoDataError`
           if no coordinates to be written.
        '''
        if ts is None:
            try:
                ts = self.ts
            except AttributeError:
                raise NoDataError("DCDWriter: no coordinate data to write to trajectory file")
        if not ts.n_atoms == self.n_atoms:
            raise ValueError("DCDWriter: Timestep does not have the correct number of atoms")
        unitcell = self.convert_dimensions_to_unitcell(ts).astype(np.float32)  # must be float32 (!)
        if not ts._pos.flags.f_contiguous:  # Not in fortran format
            ts = Timestep.from_timestep(ts)  # wrap in a new fortran formatted Timestep
        if self.convert_units:
            pos = self.convert_pos_to_native(ts._pos,
                                             inplace=False)  # possibly make a copy to avoid changing the trajectory
        self._write_next_frame(pos[:, 0], pos[:, 1], pos[:, 2], unitcell)
        self.frames_written += 1

    def convert_dimensions_to_unitcell(self, ts, _ts_order=Timestep._ts_order):
        """Read dimensions from timestep *ts* and return appropriate native unitcell.
        .. SeeAlso:: :attr:`Timestep.dimensions`
        """
        # unitcell is A,B,C,alpha,beta,gamma - convert to order expected by low level DCD routines
        lengths, angles = ts.dimensions[:3], ts.dimensions[3:]
        self.convert_pos_to_native(lengths)
        unitcell = np.zeros_like(ts.dimensions)
        # __write_DCD_frame() wants uc_array = [A, gamma, B, beta, alpha, C] and
        # will write the lengths and the angle cosines to the DCD. NOTE: It
        # will NOT write CHARMM36+ box matrix vectors. When round-tripping a
        # C36+ DCD, we loose the box vector information. However, MDAnalysis
        # itself will not detect this because the DCDReader contains the
        # heuristic to deal with either lengths/angles or box matrix on the fly.
        #
        # use np.put so that we can re-use ts_order (otherwise we
        # would need a ts_reverse_order such as _ts_reverse_order = [0, 5, 1, 4, 3, 2])
        np.put(unitcell, _ts_order, np.concatenate([lengths, angles]))
        return unitcell

    def close(self):
        """Close trajectory and flush buffers."""
        if hasattr(self, 'dcdfile') and self.dcdfile is not None:
            self._finish_dcd_write()
            self.dcdfile.close()
            self.dcdfile = None

class DCDReader(base.Reader):
    format = 'DCD'
    flavor = 'CHARMM'
    units = {'time': 'AKMA', 'length': 'Angstrom'}
    _Timestep = Timestep
    def __init__( self, dcdfilename, **kwargs ):
        super(DCDReader, self).__init__(dcdfilename, **kwargs)

        self.dcdfilename = self.filename # dcdfilename is legacy
        self.dcdfile = None  # set right away because __del__ checks

        # Issue #32 (MDanalysis): segfault if dcd is 0-size
        # Hack : test here... (but should be fixed in dcd.c)
        stats = os.stat( dcdfilename )
        if stats.st_size == 0:
            raise IOError( errno.ENODATA, "DCD file is zero size", dcdfilename )

        self.dcdfile = dcdfilename
        self.n_atoms = 0
        self.n_frames = 0
        self.fixed = 0
        self.skip = 1
        self.periodic = False

        self._read_dcd_header()

        # Convert delta to ps
        delta = mdaunits.convert(self.delta, self.units['time'], 'ps')

        self._ts_kwargs.setdefault('dt', self.skip_timestep * delta)
        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
        # Read in the first timestep
        self._read_next_timestep()


    def _dcd_header(self):  # pragma: no cover
        """Returns contents of the DCD header C structure::
             typedef struct {
               fio_fd fd;                 // FILE *
               fio_size_t header_size;    // size_t == sizeof(int)
               int natoms;
               int nsets;
               int setsread;
               int istart;
               int nsavc;
               double delta;
               int nfixed;
               int *freeind;
               float *fixedcoords;
               int reverse;
               int charmm;
               int first;
               int with_unitcell;
             } dcdhandle;
        .. deprecated:: 0.7.5
           This function only exists for debugging purposes and might
           be removed without notice. Do not rely on it.
        """
        # was broken (no idea why [orbeckst]), see Issue 27
        # 'PiiiiiidiPPiiii' should be the unpack string according to the struct.
        #    struct.unpack("LLiiiiidiPPiiii",self._dcd_C_str)
        # seems to do the job on Mac OS X 10.6.4 ... but I have no idea why,
        # given that the C code seems to define them as normal integers
        desc = [
            'file_desc', 'header_size', 'natoms', 'nsets', 'setsread', 'istart',
            'nsavc', 'delta', 'nfixed', 'freeind_ptr', 'fixedcoords_ptr', 'reverse',
            'charmm', 'first', 'with_unitcell']
        return dict(zip(desc, struct.unpack("LLiiiiidiPPiiii", self._dcd_C_str)))

    def _reopen(self):
        self.ts.frame = -1
        self._reset_dcd_read()

    def _read_next_timestep(self, ts=None):
        """Read the next frame
        .. versionchanged 0.11.0::
           Native frame read into ts._frame, ts.frame naively iterated
        """
        if ts is None:
            ts = self.ts
        ts._frame = self._read_next_frame(ts._x, ts._y, ts._z, ts._unitcell, 1)
        ts.frame += 1
        return ts

    def _read_frame(self, frame):
        """Skip to frame and read
        .. versionchanged:: 0.11.0
           Native frame read into ts._frame, ts.frame naively set to frame
        """
        self._jump_to_frame(frame)
        ts = self.ts
        ts._frame = self._read_next_frame(ts._x, ts._y, ts._z, ts._unitcell, 1)
        ts.frame = frame
        return ts

    def timeseries(self, asel, start=0, stop=-1, skip=1, format='afc'):
        """Return a subset of coordinate data for an AtomGroup
        :Arguments:
            *asel*
               :class:`~MDAnalysis.core.AtomGroup.AtomGroup` object
            *start, stop, skip*
               range of trajectory to access, start and stop are inclusive
            *format*
               the order/shape of the return data array, corresponding
               to (a)tom, (f)rame, (c)oordinates all six combinations
               of 'a', 'f', 'c' are allowed ie "fac" - return array
               where the shape is (frame, number of atoms,
               coordinates)
        """
        start, stop, skip = self.check_slice_indices(start, stop, skip)
        if len(asel) == 0:
            raise NoDataError("Timeseries requires at least one atom to analyze")
        if len(format) != 3 and format not in ['afc', 'acf', 'caf', 'cfa', 'fac', 'fca']:
            raise ValueError("Invalid timeseries format")
        atom_numbers = list(asel.indices)
        # Check if the atom numbers can be grouped for efficiency, then we can read partial buffers
        # from trajectory file instead of an entire timestep
        # XXX needs to be implemented
        return self._read_timeseries(atom_numbers, start, stop, skip, format)

    def correl(self, timeseries, start=0, stop=-1, skip=1):
        """Populate a TimeseriesCollection object with timeseries computed from the trajectory
        :Arguments:
            *timeseries*
               :class:`MDAnalysis.core.Timeseries.TimeseriesCollection`
            *start, stop, skip*
               subset of trajectory to use, with start and stop being inclusive
        """
        start, stop, skip = self.check_slice_indices(start, stop, skip)
        atomlist = timeseries._getAtomList()
        format = timeseries._getFormat()
        lowerb, upperb = timeseries._getBounds()
        sizedata = timeseries._getDataSize()
        atomcounts = timeseries._getAtomCounts()
        auxdata = timeseries._getAuxData()
        return self._read_timecorrel(atomlist, atomcounts, format, auxdata,
                                     sizedata, lowerb, upperb, start, stop, skip)

    def close(self):
        if self.dcdfile is not None:
            self._finish_dcd_read()
            #self.dcdfile.close()
            self.dcdfile = None

    def Writer(self, filename, **kwargs):
            """Returns a DCDWriter for *filename* with the same parameters as this DCD.
            All values can be changed through keyword arguments.
            :Arguments:
              *filename*
                  filename of the output DCD trajectory
            :Keywords:
              *n_atoms*
                  number of atoms
              *start*
                  number of the first recorded MD step
              *step*
                  indicate that *step* MD steps (!) make up one trajectory frame
              *delta*
                  MD integrator time step (!), in AKMA units
              *dt*
                 **Override** *step* and *delta* so that the DCD records that *dt* ps
                 lie between two frames. (It sets *step* = 1 and *delta* = ``AKMA(dt)``.)
                 The default is ``None``, in which case *step* and *delta* are used.
              *remarks*
                  string that is stored in the DCD header [XXX -- max length?]
            :Returns: :class:`DCDWriter`
            .. Note::
               The keyword arguments set the low-level attributes of the DCD
               according to the CHARMM format. The time between two frames would be
               *delta* * *step* !
            .. SeeAlso:: :class:`DCDWriter` has detailed argument description
            """
            n_atoms = kwargs.pop('n_atoms', self.n_atoms)
            kwargs.setdefault('start', self.start_timestep)
            kwargs.setdefault('step', self.skip_timestep)
            kwargs.setdefault('delta', self.delta)
            kwargs.setdefault('remarks', self.remarks)
            # dt keyword is simply passed through if provided
            return DCDWriter(filename, n_atoms, **kwargs)

    @property
    def dt(self):
        """Time between two trajectory frames in picoseconds."""
        return self.ts.dt



# Add the c functions to their respective classes so they act as class methods

try:
    import new
    DCDReader._read_dcd_header = new.instancemethod( _dcdmodule._read_dcd_header, None, DCDReader )
    DCDReader._read_next_frame = new.instancemethod( _dcdmodule._read_next_frame, None, DCDReader )
    DCDReader._jump_to_frame = new.instancemethod( _dcdmodule._jump_to_frame, None, DCDReader )
    DCDReader._reset_dcd_read = new.instancemethod( _dcdmodule._reset_dcd_read, None, DCDReader )
    DCDReader._finish_dcd_read = new.instancemethod( _dcdmodule._finish_dcd_read, None, DCDReader )
    DCDReader._read_timeseries = new.instancemethod( _dcdmodule._read_timeseries, None, DCDReader )

    DCDWriter._write_dcd_header = new.instancemethod(_dcdmodule._write_dcd_header, None, DCDWriter)
    DCDWriter._write_next_frame = new.instancemethod(_dcdmodule._write_next_frame, None, DCDWriter)
    DCDWriter._finish_dcd_write = new.instancemethod(_dcdmodule._finish_dcd_write, None, DCDWriter)

    DCDReader._read_timecorrel = new.instancemethod(dcdtimeseries.__read_timecorrel, None, DCDReader)

    del( _dcdmodule )
except ImportError:
    DCDReader._read_dcd_header = lambda self: _dcdmodule._read_dcd_header( self )
    DCDReader._read_next_frame = lambda self, x, y, z, unitcell, skip: _dcdmodule._read_next_frame( self, x, y, z, unitcell, skip )
    DCDReader._jump_to_frame = lambda self, frame: _dcdmodule._jump_to_frame( self, frame )
    DCDReader._reset_dcd_read = lambda self: _dcdmodule._reset_dcd_read( self )
    DCDReader._finish_dcd_read = lambda self: _dcdmodule._finish_dcd_read( self )
    DCDReader._read_timeseries = lambda self: _dcdmodule._read_timeseries( self )

    #DCDWriter._write_dcd_header = new.instancemethod(_dcdmodule._write_dcd_header, None, DCDWriter)
    DCDWriter._write_dcd_header = lambda self, n_atoms, start, step, delta, remarks: _dcdmodule._write_dcd_header( self, n_atoms, start, step, delta, remarks)
    DCDWriter._write_next_frame = lambda self, x, y , z, unitcell: _dcdmodule._write_next_frame( self , x, y, z, unitcell)
    DCDWriter._finish_dcd_write = lambda self: _dcdmodule._finish_dcd_write( self )

    #DCDReader._read_timecorrel = lambda self, atomlist, atomcounts, format, auxdata, sizedata, lowerb, upperb, start, stop, skip : dcdtimeseries.__read_timecorrel( self, atomlist, atomcounts, format, auxdata,sizedata, lowerb, upperb, start, stop, skip)
