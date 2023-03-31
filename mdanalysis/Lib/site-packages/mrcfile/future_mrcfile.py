# Copyright (c) 2018, Science and Technology Facilities Council
# This software is distributed under a BSD licence. See LICENSE.txt.

"""
future_mrcfile
--------------

Module which exports the :class:`FutureMrcFile` class.

Classes:
    :class:`FutureMrcFile`: An object which represents an MRC file being
        opened asynchronously.

"""

# Import Python 3 features for future-proofing
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import threading


class FutureMrcFile(object):
    """Object representing an MRC file being opened asynchronously.

    This API deliberately mimics a :class:`~concurrent.futures.Future` object
    from the :mod:`concurrent.futures` module in Python 3.2+ (which we do
    not use directly because this code still needs to run in Python 2.7).

    """

    def __init__(self, open_function, args=(), kwargs={}):
        """Initialise a new :class:`FutureMrcFile` object.

        This constructor starts a new thread which will invoke the callable
        given in ``open_function`` with the given arguments.

        Args:
            open_function: The callable to use to open the MRC file. (This will
                normally be :func:`mrcfile.open`, but could also be
                :class:`~mrcfile.mrcfile.MrcFile` or any of its subclasses.)
            args: A tuple of positional arguments to use when ``open_function``
                is called. (Normally a 1-tuple containing the name of the file
                to open.)
            kwargs: A dictionary of keyword arguments to use when
                ``open_function`` is called.
        """
        self._result_holder = [None]
        self._open_function = open_function
        self._thread = threading.Thread(target=self._run,
                                        args=args,
                                        kwargs=kwargs)
        self._thread.start()

    def _run(self, *args, **kwargs):
        """Call the open function and store the result in the holder list.

        (For internal use only.)
        """
        try:
            mrc = self._open_function(*args, **kwargs)
            self._result_holder[0] = mrc
        except Exception as ex:
            self._result_holder[0] = ex

    def cancel(self):
        """Return :data:`False`.

        (See :meth:`concurrent.futures.Future.cancel` for more details. This
        implementation does not allow jobs to be cancelled.)
        """
        return False

    def cancelled(self):
        """Return :data:`False`.

        (See :meth:`concurrent.futures.Future.cancelled` for more details.
        This implementation does not allow jobs to be cancelled.)
        """
        return False

    def running(self):
        """Return :data:`True` if the :class:`~mrcfile.mrcfile.MrcFile` is
        currently being opened.

        (See :meth:`concurrent.futures.Future.running` for more details.)
        """
        return self._thread.is_alive()

    def done(self):
        """Return :data:`True` if the file opening has finished.

        (See :meth:`concurrent.futures.Future.done` for more details.)
        """
        return not self.running()

    def result(self, timeout=None):
        """Return the :class:`~mrcfile.mrcfile.MrcFile` that has been opened.

        (See :meth:`concurrent.futures.Future.result` for more details.)

        Args:
            timeout: Time to wait (in seconds) for the file opening to finish.
                If ``timeout`` is not specified or is :data:`None`, there is no
                limit to the wait time.

        Returns:
            An :class:`~mrcfile.mrcfile.MrcFile` object (or one of its
            subclasses).

        Raises:
            :exc:`RuntimeError`: If the operation has not finished within the
                 time limit set by ``timeout``. (Note that the type of this
                 exception will change in future if this class is replaced by
                 :class:`concurrent.futures.Future`.)
            :exc:`Exception`: Any exception raised by the
                :class:`~mrcfile.mrcfile.MrcFile` opening operation will be
                re-raised here.
        """
        result = self._get_result(timeout)
        if isinstance(result, Exception):
            raise result
        else:
            return result

    def exception(self, timeout=None):
        """Return the exception raised by the file opening operation.

        (See :meth:`concurrent.futures.Future.exception` for more details.)

        Args:
            timeout: Time to wait (in seconds) for the operation to finish. If
                ``timeout`` is not specified or is :data:`None`, there is no
                limit to the wait time.

        Returns:
            An :exc:`Exception`, if one was raised by the file opening
            operation, or :data:`None` if no exception was raised.

        Raises:
            :exc:`RuntimeError`: If the operation has not finished within the
                time limit set by ``timeout``. (Note that the type of this
                exception will change in future if this class is replaced by
                :class:`concurrent.futures.Future`.)
        """
        result = self._get_result(timeout)
        if isinstance(result, Exception):
            return result
        else:
            return None

    def _get_result(self, timeout):
        """Return the result or exception from the file opening operation.

        (For internal use only.)
        """
        self._thread.join(timeout=timeout)
        if self._thread.is_alive():
            raise RuntimeError('Timed out waiting for result')
        return self._result_holder[0]

    def add_done_callback(self, fn):
        """Not implemented.

        (See :meth:`concurrent.futures.Future.add_done_callback` for more details.)
        """
        raise NotImplementedError
