'''
run with: python ten2eleven.py -f torsions test_dummy_old_MDA_code.py 
Author: Tyler Reddy
'''
from __future__ import absolute_import

from lib2to3.fixer_base import BaseFix
from lib2to3.fixer_util import Name, Call, LParen, RParen, ArgList, Dot
from lib2to3 import pytree


class FixTorsions(BaseFix):

    PATTERN = """
        power< head =any+
                trailer< dot = '.' method='torsions'>
                                tail=any*>
    """

    def transform(self, node, results):
        head = results['head']
        method = results['method']
        tail = results['tail']
        syms = self.syms
        method_name = 'dihedrals'
        head = [n.clone() for n in head]
        tail = [n.clone() for n in tail]
        args = head + [pytree.Node(syms.trailer, [Dot(), Name(method_name, prefix = method.prefix)])]
        new = pytree.Node(syms.power, args)
        return new
