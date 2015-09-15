import errno

from . import base
from ..lib.formats.xtc import XTCFile
from ..lib.mdamath import triclinic_box, triclinic_vectors


class XTCReader(base.Reader):
    format = 'XTC'
    units = {'time': 'ps', 'length': 'nm'}

    def __init__(self, filename, convert_units=True, sub=None, **kwargs):
        super(XTCReader, self).__init__(filename, convert_units=convert_units,
                                        **kwargs)
        self._xtc = XTCFile(filename)

        self.n_frames = len(self._xtc)
        xyz, box, step, time, prec = self._xtc.read()
        try:
            xtc_frame = self._xtc.read()
            time_2 = xtc_frame[3]
        except:
            time_2 = 2 * time
        self._sub = sub
        if sub is not None:
            self.n_atoms = len(sub)
        else:
            self.n_atoms = self._xtc.n_atoms

        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
        self.ts.frame = 0  # 0-based frame number as starting frame
        self._frame = -1
        if self._sub is not None:
            xyz = xyz[self._sub]
        self.ts._pos = xyz
        self.ts.time = time
        self.ts.dt = time_2 - time
        self.ts.dimensions = triclinic_box(*box)

        if self.convert_units:
            self.convert_pos_from_native(self.ts._pos)
            self.convert_pos_from_native(self.ts._unitcell[:3])

        self.rewind()

    def close(self):
        self._xtc.close()

    def Writer(self, filename, n_atoms=None, **kwargs):
        if n_atoms is None:
            n_atoms = self.n_atoms
        return XTCWriter(filename, n_atoms=n_atoms, **kwargs)

    def rewind(self):
        self._read_frame(0)

    def _reopen(self):
        self.ts.frame = 0
        self._frame = -1
        self._xtc.close()
        self._xtc.open(self._xtc.fname, 'r')

    def _read_frame(self, i):
        self._xtc.seek(i)
        self.ts.frame = i - 1
        self._frame = i - 1
        return self._read_next_timestep()

    def _read_next_timestep(self, ts=None):
        if self._frame == self.n_frames - 1:
            raise IOError(errno.EIO, 'trying to go over trajectory limit')
        xyz, box, step, time, prec = self._xtc.read()
        self._frame += 1

        if ts is None:
            ts = self.ts

        if self._sub is not None:
            xyz = xyz[self._sub]

        ts._pos = xyz
        ts.frame = self._frame
        ts.time = time
        if self.convert_units:
            self.convert_pos_from_native(ts._pos)

        return ts


class XTCWriter(base.Writer):
    format = 'XTC'
    units = {'time': None, 'length': 'nm'}

    def __init__(self, filename, n_atoms, convert_units=True, **kwargs):
        self._convert_units = convert_units
        self._xtc = XTCFile(filename, 'w')
        self.n_atoms = n_atoms

    def close(self):
        self._xtc.close()

    def __del__(self):
        self.close()

    def write_next_timestep(self, ts):
        xyz = ts._pos.copy()
        time = ts.time
        step = ts.frame
        dimensions = ts.dimensions

        if self._convert_units:
            xyz = self.convert_pos_to_native(xyz, inplace=False)
            dimensions = self.convert_dimensions_to_unitcell(ts, inplace=False)

        box = triclinic_vectors(dimensions)

        self._xtc.write(xyz, box, step, time)
