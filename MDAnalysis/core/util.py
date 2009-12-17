# $Id$
"""Small helper functions that don't fit anywhere else."""


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
