'''
run with: python ten2eleven.py -f agcountmethods test_dummy_old_MDA_code.py 
Author: Tyler Reddy
'''
from __future__ import absolute_import

from lib2to3.fixer_base import BaseFix
from lib2to3.fixer_util import Name, Call, LParen, RParen, ArgList, Dot
from lib2to3 import pytree


class FixAgcountmethods(BaseFix):

    PATTERN = """
        power< head =any+
                trailer< dot = '.' method=('numberOfAtoms'|'numberOfResidues'|'numberOfSegments')>
                                parens=trailer< '(' ')' >
                                tail=any*>
    """

    def transform(self, node, results):
        head = results['head']
        method = results['method'][0]
        tail = results['tail']
        syms = self.syms
        method_name = method.value
        head = [n.clone() for n in head]
        tail = [n.clone() for n in tail]
        replacement_dict = {'numberOfAtoms': 'n_atoms', 'numberOfResidues': 'n_residues', 'numberOfSegments':'n_segments'}
        method_name = replacement_dict[method_name]
        args = head + [pytree.Node(syms.trailer, [Dot(), Name(method_name, prefix = method.prefix)])]
        new = pytree.Node(syms.power, args)
        return new
