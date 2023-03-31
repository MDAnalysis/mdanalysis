from ctypes import c_char_p, c_uint64

from .frame import Frame, Topology
from .misc import ChemfilesError
from .utils import CxxPointer, _call_with_growing_buffer


class BaseTrajectory(CxxPointer):
    def __init__(self, ptr):
        self.__closed = False
        super(BaseTrajectory, self).__init__(ptr, is_const=False)

    def __check_opened(self):
        if self.__closed:
            raise ChemfilesError("Can not use a closed Trajectory")

    def __del__(self):
        if not self.__closed:
            self.close()

    def __enter__(self):
        self.__check_opened()
        return self

    def __exit__(self, *args):
        self.close()

    def __iter__(self):
        self.__check_opened()
        for step in range(self.nsteps):
            yield self.read_step(step)

    def read(self):
        """
        Read the next step of this :py:class:`Trajectory` and return the
        corresponding :py:class:`Frame`.
        """
        self.__check_opened()
        frame = Frame()
        self.ffi.chfl_trajectory_read(self.mut_ptr, frame.mut_ptr)
        return frame

    def read_step(self, step):
        """
        Read a specific ``step`` in this :py:class:`Trajectory` and return the
        corresponding :py:class:`Frame`.
        """
        self.__check_opened()
        frame = Frame()
        self.ffi.chfl_trajectory_read_step(self.mut_ptr, c_uint64(step), frame.mut_ptr)
        return frame

    def write(self, frame):
        """Write a :py:class:`Frame` to this :py:class:`Trajectory`."""
        self.__check_opened()
        self.ffi.chfl_trajectory_write(self.mut_ptr, frame.ptr)

    def set_topology(self, topology, format=""):
        """
        Set the :py:class:`Topology` associated with this :py:class:`Trajectory`.

        The new topology will be used when reading and writing the files,
        replacing any topology in the frames or files.

        If the ``topology`` parameter is a :py:class:`Topology` instance, it is
        used directly. If the ``topology`` parameter is a string, the first
        :py:class:`Frame` of the corresponding file is read, and the topology of
        this frame is used.

        When reading from a file, if ``format`` is not the empty string, it is
        used as the file format instead of guessing it from the file extension.
        """
        self.__check_opened()
        if isinstance(topology, Topology):
            self.ffi.chfl_trajectory_set_topology(self.mut_ptr, topology.ptr)
        else:
            self.ffi.chfl_trajectory_topology_file(
                self.mut_ptr, topology.encode("utf8"), format.encode("utf8")
            )

    def set_cell(self, cell):
        """
        Set the :py:class:`UnitCell` associated with this :py:class:`Trajectory`
        to a copy of ``cell``.

        This :py:class:`UnitCell` will be used when reading and writing the
        files, replacing any unit cell in the frames or files.
        """
        self.__check_opened()
        self.ffi.chfl_trajectory_set_cell(self.mut_ptr, cell.ptr)

    @property
    def nsteps(self):
        """Get the current number of steps in this :py:class:`Trajectory`."""
        self.__check_opened()
        nsteps = c_uint64()
        self.ffi.chfl_trajectory_nsteps(self.mut_ptr, nsteps)
        return nsteps.value

    @property
    def path(self):
        """Get the path used to open this :py:class:`Trajectory`."""
        self.__check_opened()
        return _call_with_growing_buffer(
            lambda buffer, size: self.ffi.chfl_trajectory_path(self.ptr, buffer, size),
            initial=256,
        )

    def close(self):
        """
        Close this :py:class:`Trajectory` and write any buffered content to the
        file.
        """
        self.__check_opened()
        self.__closed = True
        self.ffi.chfl_trajectory_close(self.ptr)


class Trajectory(BaseTrajectory):
    """
    A :py:class:`Trajectory` represent a physical file from which we can read
    :py:class:`Frame`.
    """

    def __init__(self, path, mode="r", format=""):
        """
        Open the file at the given ``path`` using the given ``mode`` and
        optional file ``format``.

        Valid modes are ``'r'`` for read, ``'w'`` for write and ``'a'`` for
        append.

        The ``format`` parameter is needed when the file format does not match
        the extension, or when there is not standard extension for this format.
        If `format` is an empty string, the format will be guessed from the
        file extension.
        """
        ptr = self.ffi.chfl_trajectory_with_format(
            path.encode("utf8"), mode.encode("utf8"), format.encode("utf8")
        )
        # Store mode and format for __repr__
        self.__mode = mode
        self.__format = format
        super(Trajectory, self).__init__(ptr)

    def __repr__(self):
        return f"Trajectory('{self.path}', '{self.__mode}', '{self.__format}')"


class MemoryTrajectory(BaseTrajectory):
    """
    A :py:class:`MemoryTrajectory` allow to read/write in-memory data as though
    it was a formatted file.
    """

    def __init__(self, data="", mode="r", format=""):
        """
        The ``format`` parameter is always required.

        When reading (``mode`` is ``'r'``), the ``data`` parameter will be used
        as the formatted file.

        When writing (``mode`` is ``'w'``), the ``data`` parameter is ignored.
        To get the memory buffer containing everything already written, use the
        :py:func:`buffer` function.
        """

        if not format:
            raise ChemfilesError(
                "'format' is required when creating a MemoryTrajectory"
            )

        if mode == "r":
            if isinstance(data, str):
                data = data.encode("utf8")
            elif not isinstance(data, bytes):
                raise ChemfilesError("the 'data' parameter must be a string")

            ptr = self.ffi.chfl_trajectory_memory_reader(
                data, len(data), format.encode("utf8")
            )
        elif mode == "w":
            ptr = self.ffi.chfl_trajectory_memory_writer(format.encode("utf8"))
        else:
            raise ChemfilesError(f"invalid mode '{mode}' passed to MemoryTrajectory")

        super(MemoryTrajectory, self).__init__(ptr)

    def __repr__(self):
        return f"MemoryTrajectory({self.__mode}', '{self.__format}')"

    def buffer(self):
        """
        Get the data written to this in-memory trajectory. This is not valid to
        call when reading in-memory data.
        """
        buffer = c_char_p()
        size = c_uint64()
        self.ffi.chfl_trajectory_memory_buffer(self.ptr, buffer, size)

        return buffer.value
