import errno

from . import base
from ..lib.formats.trr import TRRFile
from ..lib.mdamath import triclinic_box, triclinic_vectors


class TRRReader(base.Reader):
    format = 'TRR'
    units = {'time': 'ps', 'length': 'nm', 'velocity': 'nm/ps',
             'force': 'kJ/(mol*nm)'}

    def __init__(self, filename, convert_units=True, sub=None, **kwargs):
        super(TRRReader, self).__init__(filename, convert_units=convert_units,
                                        **kwargs)
        self._trr = TRRFile(filename)

        self.n_frames = len(self._trr)
        self._sub = sub
        trr_frame = self._trr.read()
        try:
            trr_frame_2 = self._trr.read()
            time_2 = trr_frame_2.time
        except:
            time_2 = trr_frame.time

        if sub is not None:
            self.n_atoms = len(sub)
        else:
            self.n_atoms = self._trr.n_atoms
        self._trr.seek(0)

        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
        self._frame = 0
        self.ts = self.frame_to_ts(trr_frame, self.ts)

        self.ts.dt = time_2 - trr_frame.time
        self.ts.dimensions = triclinic_box(*trr_frame.box)

        if self.convert_units:
            self.convert_pos_from_native(self.ts.dimensions[:3])

        self.rewind()

    def frame_to_ts(self, frame, ts):
        ts.has_positions = frame.hasx
        ts.has_velocities = frame.hasv
        ts.has_forces = frame.hasf

        if ts.has_positions:
            if self._sub is not None:
                ts.positions = frame.x[self._sub]
            else:
                ts.positios = frame.x
            if self.convert_units:
                self.convert_pos_from_native(ts.positions)

        if ts.has_velocities:
            if self._sub is not None:
                ts.velocities = frame.v[self._sub]
            else:
                ts.velocities = frame.v
            if self.convert_units:
                self.convert_velocities_from_native(ts.velocities)

        if ts.has_forces:
            if self._sub is not None:
                ts.forces = frame.f[self._sub]
            else:
                ts.forces = frame.f
            if self.convert_units:
                self.convert_forces_from_native(ts.forces)

        ts.time = frame.time
        ts.frame = self._frame
        ts.data['step'] = frame.step

        return ts

    def close(self):
        self._trr.close()

    def Writer(self, filename, n_atoms=None, **kwargs):
        if n_atoms is None:
            n_atoms = self.n_atoms
        return TRRWriter(filename, n_atoms=n_atoms, **kwargs)

    def rewind(self):
        self._read_frame(0)

    def _reopen(self):
        self.ts.frame = 0
        self._frame = -1
        self._trr.close()
        self._trr.open(self._trr.fname, 'r')

    def _read_frame(self, i):
        self._trr.seek(i)
        self._frame = i - 1
        return self._read_next_timestep()

    def _read_next_timestep(self, ts=None):
        if self._frame == self.n_frames - 1:
            raise IOError(errno.EIO, 'trying to go over trajectory limit')
        if ts is None:
            ts = self.ts
        trr_frame = self._trr.read()
        self._frame += 1
        ts = self.frame_to_ts(trr_frame, ts)

        return ts


class TRRWriter(base.Writer):
    format = 'TRR'
    units = {'time': 'ps', 'length': 'nm', 'velocity': 'nm/ps',
             'force': 'kJ/(mol*nm)'}

    def __init__(self, filename, n_atoms, convert_units=True, **kwargs):
        self._convert_units = convert_units
        self._trr = TRRFile(filename, 'w')
        self.n_atoms = n_atoms

    def close(self):
        self._trr.close()

    def __del__(self):
        self.close()

    def write_next_timestep(self, ts):
        xyz = None
        if ts.has_positions:
            xyz = ts.positions.copy()
            if self._convert_units:
                self.convert_pos_to_native(xyz)

        velo = None
        if ts.has_velocities:
            velo = ts.velocities.copy()
            if self._convert_units:
                self.convert_velocities_to_native(velo)

        forces = None
        if ts.has_forces:
            forces = ts.forces.copy()
            if self._convert_units:
                self.convert_forces_to_native(forces)

        time = ts.time
        step = ts.frame

        if self._convert_units:
            dimensions = self.convert_dimensions_to_unitcell(ts, inplace=False)

        box = triclinic_vectors(dimensions)

        self._trr.write(xyz, velo, forces, box, step, time, 1, self.n_atoms)
