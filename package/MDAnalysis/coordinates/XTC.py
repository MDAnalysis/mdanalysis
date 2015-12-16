from .XDR import XDRBaseReader, XDRBaseWriter
from ..lib.formats.xdrlib import XTCFile
from ..lib.mdamath import triclinic_vectors


class XTCWriter(XDRBaseWriter):
    format = 'XTC'
    units = {'time': 'ps', 'length': 'nm'}
    _file = XTCFile

    def write_next_timestep(self, ts):
        """Write timestep object into trajectory.

        Parameters
        ----------
        ts: TimeStep

        See Also
        --------
        <FormatWriter>.write(AtomGroup/Universe/TimeStep)
        The normal write() method takes a more general input
        """
        xyz = ts.positions.copy()
        time = ts.time
        step = ts.frame
        dimensions = ts.dimensions

        if self._convert_units:
            xyz = self.convert_pos_to_native(xyz, inplace=False)
            dimensions = self.convert_dimensions_to_unitcell(ts, inplace=False)

        box = triclinic_vectors(dimensions)

        self._xdr.write(xyz, box, step, time)


class XTCReader(XDRBaseReader):
    format = 'XTC'
    units = {'time': 'ps', 'length': 'nm'}
    _writer = XTCWriter
    _file = XTCFile

    def _frame_to_ts(self, frame, ts):
        """convert a xtc-frame to a mda TimeStep"""
        ts.frame = self._frame
        ts.time = frame.time
        ts.data['step'] = frame.step

        if self._sub is not None:
            ts.positions = frame.x[self._sub]
        else:
            ts.positions = frame.x
        if self.convert_units:
            self.convert_pos_from_native(ts.positions)

        return ts
