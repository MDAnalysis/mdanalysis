from __future__ import absolute_import

from pkg_resources import resource_filename

__all__ = ["DX", "CCP4", "gOpenMol"]

DX = resource_filename(__name__, 'test.dx')
DXGZ = resource_filename(__name__, 'test.dx.gz')
CCP4 = resource_filename(__name__, 'test.ccp4')
# from http://www.ebi.ac.uk/pdbe/coordinates/files/1jzv.ccp4
# (see issue #57)
CCP4_1JZV = resource_filename(__name__, '1jzv.ccp4')
# from https://github.com/ccpem/mrcfile/blob/master/tests/test_data/EMD-3001.map.bz2
MRC_EMD3001 = resource_filename(__name__, 'EMD-3001.map.bz2')
# water density around M2 TM helices of nAChR from MD simulations
# [O. Beckstein and M. S. P. Sansom. Physical Biology 3(2):147-159, 2006]
gOpenMol = resource_filename(__name__, 'nAChR_M2_water.plt')
