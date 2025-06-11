from transformations import *

from .mdamath import rotaxis

from .util import deprecate

rotaxis = deprecate(
    rotaxis,
    old_name="MDAnalysis.lib.transformations.rotaxis",
    new_name="MDAnalysis.lib.mdamath.rotaxis",
    release="2.10.0",
    remove="3.0.0",
)
