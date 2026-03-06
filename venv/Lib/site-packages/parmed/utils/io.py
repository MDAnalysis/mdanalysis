"""
Tools to aid in input/output within the parmed package
"""

__all__ = ['genopen']

import os
from io import TextIOWrapper, BytesIO
from pathlib import Path
from typing import Union
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

from ..constants import DEFAULT_ENCODING

def genopen(name: Union[Path, str], mode: str = 'r'):
    """
    Opens a file, automatically detecting compression schemes by filename
    extension. Note, these files are opened in a way that *always* returns a
    string. This is an important distinction in Python 3 where many file-like
    objects return bytes instead of strings. This is detected and handled
    properly so that the object returned from this function is always
    string-based (so you cannot write or read bytes directly from a file opened
    via ``genopen``).

    This routine also recognizes URLs and will read remote files when given a
    URL starting with either http:// or https://. Like with standard local file
    names, compression is automatically detected by filename extension, and both
    gzip and bzip2 files are supported.

    Parameters
    ----------
    name : str or Path
        Name of the file to open or URL to a remote file to access
    mode : str, optional
        Whether to open the file to 'r'ead, 'w'rite, or 'a'ppend. Default is 'r'

    Returns
    -------
    file : file-like
        A file-like object in the requested mode

    Notes
    -----
    Python's BZ2File does not support writing to ``append`` mode (mode='a'), so
    it is faked here. The entire file contents are read into memory and then
    written into a 'new' file with the same name as the original. As such, it is
    noticeably slower and more resource-intensive (particularly for large files)
    than using gzipped files.
    """
    if mode not in ['w', 'r', 'a']:
        raise ValueError('open mode must be "w", "r", or "a"')

    # Handle arbitrary online files. file:// is just an alias for a local file
    is_url = False
    if name.startswith('file:///'):
        name = name[7:]
    elif name.startswith('http://') or name.startswith('https://') or name.startswith('ftp://'):
        is_url = True
        if mode in ['w', 'a']:
            raise ValueError('Cannot write or append a webpage')
        try:
            open_url = urlopen(name)
        except (HTTPError, URLError) as err:
            raise IOError(f"Could not open {name}") from err

    if name.endswith('.bz2'):
        import bz2
        # BZ2File cannot open in append mode, so we have to fake it. Read the
        # entire existing contents into memory, open a new file, write the
        # contents back, and return the file that is now open for writing
        if mode == 'a':
            tmp = BytesIO()
            if os.path.exists(name):
                with bz2.BZ2File(name, 'rb') as f:
                    tmp.write(f.read())
                tmp.seek(0)
            f = bz2.BZ2File(name, 'wb')
            f.write(tmp.read())
            return TextIOWrapper(f)
        # If it is a URL, just pass in the urlopen object as a filename
        if is_url:
            name = open_url
        return TextIOWrapper(bz2.BZ2File(name, mode+'b'))
    elif name.endswith('.gz'):
        import gzip
        if is_url:
            return TextIOWrapper(gzip.GzipFile(fileobj=open_url, mode='r'))
        else:
            return TextIOWrapper(gzip.open(name, mode+'b'))

    if is_url:
        return TextIOWrapper(open_url)
    else:
        return open(name, mode, encoding=DEFAULT_ENCODING)
