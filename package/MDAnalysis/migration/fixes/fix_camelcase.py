'''
run with: python ten2eleven.py -f camelcase test_dummy_old_MDA_code.py 
Author: Tyler Reddy
'''

from lib2to3.fixer_base import BaseFix
from lib2to3.fixer_util import Name, Call, LParen, RParen, ArgList, Dot
from lib2to3 import pytree
import re


class FixCamelcase(BaseFix):

    PATTERN = """
                trailer< dot = '.' method=('totalMass'|'totalCharge'|'centerOfGeometry'|
                                         'radiusOfGyration'|'shapeParameter'|'momentOfInertia'|
                                         'principalAxes'|'packIntoBox'|'asUniverse'|
                                         'applyPBC' | 'align_principalAxis' | 'centerOfMass')>
    """

    def transform(self, node, results):
        method = results['method'][0]
        method_name = method.value
        #with the exception of applyPBC, all camelcase changes here involve replacing a capital letter with an underscore and the lower case version of that letter
        if not method_name == 'applyPBC':
            method_name = re.sub(r'([A-Z]{1})', r'_\1', method_name).lower()
        else:
            method_name = 'apply_PBC'

        syms = self.syms
        args = [pytree.Node(syms.trailer, [Dot(), Name(method_name)])]
        new = pytree.Node(syms.power, args)
        return new

