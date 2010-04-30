"""
Status codes and symbols
========================

The module makes all the status symbols available that are used in
:mod:`libxdrfile`. The value of each symbol is an integer (as defined
in ``xdrfile.h``).

The dictionary :data:`xdrfile.errno.errorcode` maps numeric codes to
symbol names.

.. data:: errorcode
.. data:: errorsymbols
"""

import libxdrfile

#: List of all error symbols ``exdr*`` extracted from :mod:`libxdrfile`.
errorsymbols = [k for k in libxdrfile.__dict__.keys() if k[:4] == 'exdr']

#: Dictionary that maps error codes to symbol names.
errorcode = dict(((libxdrfile.__dict__[k], k) for k in errorsymbols))

globals().update(dict((errorcode[n], n) for n in errorcode))

