# emacs: at the end of the file
# ex: set sts=4 ts=4 sw=4 et:
# ## ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### #
"""Stub file for a guaranteed safe import of duecredit constructs:  if duecredit
is not available.

To use it, place it into your project codebase to be imported, e.g. copy as

    cp stub.py /path/tomodule/module/due.py

Note that it might be better to avoid naming it duecredit.py to avoid shadowing
installed duecredit.

Then use in your code as

    from .due import due, Doi, BibTeX

See  https://github.com/duecredit/duecredit/blob/master/README.md for examples.

Origin:     Originally a part of the duecredit
Copyright:  2015-2016  DueCredit developers
License:    BSD-2

Modified for MDAnalysis to avoid calls to fork which raises cryptic
warning under MPI(see PR #1794 for rationale)

"""

from __future__ import absolute_import
__version__ = '0.0.5'


class InactiveDueCreditCollector(object):
    """Just a stub at the Collector which would not do anything"""
    def _donothing(self, *args, **kwargs):
        """Perform no good and no bad"""
        pass # pylint: disable=unnecessary-pass

    def dcite(self, *args, **kwargs):
        """If I could cite I would"""
        def nondecorating_decorator(func):
            return func
        return nondecorating_decorator

    cite = load = add = _donothing

    def __repr__(self):
        return self.__class__.__name__ + '()'


def _donothing_func(*args, **kwargs):
    """Perform no good and no bad"""
    pass # pylint: disable=unnecessary-pass


try:
    # Avoid call of fork inside duecredit; see
    # https://github.com/MDAnalysis/mdanalysis/pull/1822#issuecomment-373009050
    import sys
    import os
    if sys.version_info >= (3, 7):
        import duecredit
    else:
        from mock import patch
        if sys.version_info <= (2, ):
            from contextlib import nested
            with nested(patch('os.fork'), patch('os.popen')) \
                 as (os_dot_fork, os_dot_popen):
                import duecredit
        else:
            if not os.name == 'nt':
                with patch('os.fork') as os_dot_fork, patch('os.popen') as os_dot_popen:
                    import duecredit
            else:
                # Windows doesn't have os.fork
                import duecredit

    from duecredit import due, BibTeX, Doi, Url
    if 'due' in locals() and not hasattr(due, 'cite'):
        raise RuntimeError(
            "Imported due lacks .cite. DueCredit is now disabled")
except Exception as err:
    if not isinstance(err, ImportError):
        import logging
        import warnings
        errmsg = "Failed to import duecredit due to {}".format(str(err))
        warnings.warn(errmsg)
        logging.getLogger("duecredit").error(
            "Failed to import duecredit due to {}".format(str(err)))
    # else:
    #   Do not issue any warnings if duecredit is not installed;
    #   this is the user's choice (Issue #1872)

    # Initiate due stub
    due = InactiveDueCreditCollector()
    BibTeX = Doi = Url = _donothing_func

# Emacs mode definitions
# Local Variables:
# mode: python
# py-indent-offset: 4
# tab-width: 4
# indent-tabs-mode: nil
# End:
