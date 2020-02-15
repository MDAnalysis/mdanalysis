
from __future__ import absolute_import
from ...due import due, Doi
from .hole import hole, HoleAnalysis

due.cite(Doi("10.1016/S0006-3495(93)81293-1"),
         description="HOLE program",
         path="MDAnalysis.analysis.hole",
         cite_module=True)
due.cite(Doi("10.1016/S0263-7855(97)00009-X"),
         description="HOLE program",
         path="MDAnalysis.analysis.hole",
         cite_module=True)
due.cite(Doi("10.1016/j.jmb.2013.10.024"),
         description="HOLE trajectory analysis with orderparameters",
         path="MDAnalysis.analysis.hole",
         cite_module=True)
del Doi
