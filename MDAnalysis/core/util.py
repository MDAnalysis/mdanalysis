# $Id$
"""
Helper functions
================

Small helper functions that don't fit anywhere else.

.. autofunction:: filename
.. autofunction:: iterable
.. autofunction:: asiterable
"""


import os.path
def filename(name,ext=None,keep=False):
    """Return a new name that has suffix attached; replaces other extensions.

    name        filename; extension is replaced unless keep=True
    ext         extension
    keep        False: replace existing extension; True: keep if exists
    """ 
    name = str(name)
    if ext is None:
        return name
    if not ext.startswith(os.path.extsep):
        ext = os.path.extsep + ext
    #if name.find(ext) > 0:    # normally >= 0 but if we start with '.' then we keep it
    #    return name
    root, origext = os.path.splitext(name)
    if len(origext) == 0 or not keep:
        return root + ext
    return name


def iterable(obj):
    """Returns ``True`` if *obj* can be iterated over and is *not* a  string."""
    if type(obj) is str:
        return False    # avoid iterating over characters of a string

    if hasattr(obj, 'next'):
        return True    # any iterator will do 
    try: 
        len(obj)       # anything else that might work
    except TypeError: 
        return False
    return True

def asiterable(obj):
    """Returns obj so that it can be iterated over; a string is *not* treated as iterable"""
    if not iterable(obj):
        obj = [obj]
    return obj
