"""
This package contains code for reading TINKER-style parameter and system-setup
files (like output from "analyze" and the xyz-format file).
"""

__all__ = ['system', 'parameterfile', 'tinkerfiles', 'XyzFile']

from .tinkerfiles import XyzFile
