import errno
import numpy as np
from os.path import getctime, getsize, isfile, split, join
import warnings

from . import base
from ..lib.mdamath import triclinic_box


def offsets_filename(filename, ending='npz'):
    head, tail = split(filename)
    return join(head, '.{}_offsets.{}'.format(tail, ending))


class XDRBaseReader(base.Reader):
    def __init__(self, filename, convert_units=True, sub=None,
                 refresh_offsets=False, **kwargs):
        super(XDRBaseReader, self).__init__(filename,
                                            convert_units=convert_units,
                                            **kwargs)
        self._xdr = self._file(self.filename)

        frame = self._xdr.read()
        try:
            xdr_frame = self._xdr.read()
            dt = xdr_frame.time - frame.time
            # avoid triggering a offset calculation!
            self._xdr.seek(0)
            self._xdr.read()
        except:
            dt = 0

        self._sub = sub
        if self._sub is not None:
            self.n_atoms = len(self._sub)
        else:
            self.n_atoms = self._xdr.n_atoms

        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
        self._frame = 0
        self._frame_to_ts(frame, self.ts)
        # these should only be initialized once
        self.ts.dt = dt
        self.ts.dimensions = triclinic_box(*frame.box)
        if self.convert_units:
            self.convert_pos_from_native(self.ts._unitcell[:3])

        if not refresh_offsets:
            self._load_offsets()
        else:
            self._read_offsets(store=True)

    def close(self):
        self._xdr.close()

    def _load_offsets(self):
        fname = offsets_filename(self.filename)

        if not isfile(fname):
            self._read_offsets(store=True)
            return

        with open(fname) as f:
            data = np.load(f)
            offsets = data['offsets']
            ctime = data['ctime']
            size = data['size']

        ctime_ok = getctime(self.filename) == ctime
        size_ok = getsize(self.filename) == size

        if not (ctime_ok and size_ok):
            warnings.warn("Aborted loading offsets from file\n "
                          "ctime or size did not match")
            self._read_offsets(store=True)
        else:
            self._xdr.set_offsets(offsets)

    def _read_offsets(self, store=False):
        offsets = self._xdr.offsets
        if store:
            ctime = getctime(self.filename)
            size = getsize(self.filename)
            np.savez(offsets_filename(self.filename),
                     offsets=offsets, size=size, ctime=ctime)

    def rewind(self):
        """Read the first frame again"""
        self._read_frame(0)

    @property
    def n_frames(self):
        return len(self._xdr)

    def _reopen(self):
        self.ts.frame = 0
        self._frame = -1
        self._xdr.close()
        self._xdr.open(self.filename, 'r')

    def _read_frame(self, i):
        self._xdr.seek(i)
        self.ts.frame = i - 1
        self._frame = i - 1
        return self._read_next_timestep()

    def _read_next_timestep(self, ts=None):
        if self._frame == self.n_frames - 1:
            raise IOError(errno.EIO, 'trying to go over trajectory limit')
        if ts is None:
            ts = self.ts
        frame = self._xdr.read()
        self._frame += 1
        self._frame_to_ts(frame, ts)
        return ts

    def Writer(self, filename, n_atoms=None, **kwargs):
        """Return writer for trajectory format"""
        if n_atoms is None:
            n_atoms = self.n_atoms
        return self._writer(filename, n_atoms=n_atoms, **kwargs)


class XDRBaseWriter(base.Writer):

    def __init__(self, filename, n_atoms, convert_units=True, **kwargs):
        self.filename = filename
        self._convert_units = convert_units
        self.n_atoms = n_atoms
        self._xdr = self._file(self.filename, 'w')

    def close(self):
        self._xdr.close()

    def __del__(self):
        self.close()
