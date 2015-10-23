from pkg_resources import resource_filename


def get_file(fname):
    """return full path to trajectory named 'fname'"""
    return resource_filename(__name__, fname)
