# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
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
#     Molecular Dynamics Simulations. J. Comput. Chem. (2011),
#     in press.
#

"""
Setting up logging --- :mod:`MDAnalysis.core.log`
====================================================

Configure logging for MDAnalysis. Import this module if logging is
desired in application code.

Logging to a file and the console is set up by default as described
under `logging to multiple destinations`_.

The top level logger of the library is named *MDAnalysis* by
convention; a simple logger that writes to the console and logfile can
be created with the :func:`create` function. This only has to be done
*once*. For convenience, the default MDAnalysis logger can be created
with :func:`MDAnalysis.start_logger`::

 import MDAnalysis
 MDAnalysis.start_logger()

Once this has been done, MDAnalysis will write messages to the logfile
(named `MDAnalysis.log` by default but this can be changed with the
optional argument to :func:`~MDAnalysis.start_logger`).

Any code can log to the MDAnalysis logger by using ::

 import logging
 logger = logging.getLogger('MDAnalysis.MODULENAME')

 # use the logger, for example at info level:
 logger.info("Starting task ...")

The important point is that the name of the logger begins with
"MDAnalysis.".

.. _logging to multiple destinations:
   http://docs.python.org/library/logging.html?#logging-to-multiple-destinations

.. SeeAlso:: The :mod:`logging` module in the standard library contains
             in depth documentation about using logging.


Convenience functions
---------------------

Two convenience functions at the top level make it easy to start and
stop the default *MDAnalysis* logger.

.. autofunction:: MDAnalysis.start_logging
.. autofunction:: MDAnalysis.stop_logging


Functions and classes
---------------------

"""

import logging

def create(logger_name="MDAnalysis", logfile="MDAnalysis.log"):
    """Create a top level logger.

    - The file logger logs everything (including DEBUG).
    - The console logger only logs INFO and above.

    Logging to a file and the console as described under `logging to
    multiple destinations`_.

    The top level logger of MDAnalysis is named *MDAnalysis*.  Note
    that we are configuring this logger with console output. If a root
    logger also does this then we will get two output lines to the
    console.

    .. _logging to multiple destinations:
       http://docs.python.org/library/logging.html?#logging-to-multiple-destinations
    """

    logger = logging.getLogger(logger_name)

    logger.setLevel(logging.DEBUG)

    # handler that writes to logfile
    logfile_handler = logging.FileHandler(logfile)
    logfile_formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    logfile_handler.setFormatter(logfile_formatter)
    logger.addHandler(logfile_handler)

    # define a Handler which writes INFO messages or higher to the sys.stderr
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    return logger

def clear_handlers(logger):
    """clean out handlers in the library top level logger

    (only important for reload/debug cycles...)
    """
    for h in logger.handlers:
        logger.removeHandler(h)

class NullHandler(logging.Handler):
    """Silent Handler.

    Useful as a default::

      h = NullHandler()
      logging.getLogger("MDAnalysis").addHandler(h)
      del h

    """
    def emit(self, record):
        pass


