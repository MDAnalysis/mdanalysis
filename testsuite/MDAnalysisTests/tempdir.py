# License
# -------
# MIT License
#
# Copyright (c) 2010-2016 Thomas Fenzl
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
import os
import tempfile
import shutil
from functools import wraps


class TempDir(object):
    """ class for temporary directories
creates a (named) directory which is deleted after use.
All files created within the directory are destroyed
Might not work on windows when the files are still opened
"""
    def __init__(self, suffix="", prefix="tmp", basedir=None):
        self.name = tempfile.mkdtemp(suffix=suffix, prefix=prefix, dir=basedir)

    def __del__(self):
        try:
            if self.name:
                self.dissolve()
        except AttributeError:
            pass

    def __enter__(self):
        return self.name

    def __exit__(self, *errstuff):
        self.dissolve()

    def dissolve(self):
        """remove all files and directories created within the tempdir"""
        if self.name:
            shutil.rmtree(self.name)
        self.name = ""

    def __str__(self):
        if self.name:
            return "temporary directory at: {}".format(self.name,)
        else:
            return "dissolved temporary directory"


class in_tempdir(object):
    """Create a temporary directory and change to it.  """

    def __init__(self, *args, **kwargs):
        self.tmpdir = TempDir(*args, **kwargs)

    def __enter__(self):
        self.old_path = os.getcwd()
        os.chdir(self.tmpdir.name)
        return self.tmpdir.name

    def __exit__(self, *errstuff):
        os.chdir(self.old_path)
        self.tmpdir.dissolve()


def run_in_tempdir(*args, **kwargs):
    """Make a function execute in a new tempdir.
    Any time the function is called, a new tempdir is created and destroyed.
    """
    def change_dird(fnc):
        @wraps(fnc)
        def wrapper(*funcargs, **funckwargs):
            with in_tempdir(*args, **kwargs):
                return fnc(*funcargs, **funckwargs)
        return wrapper
    return change_dird
