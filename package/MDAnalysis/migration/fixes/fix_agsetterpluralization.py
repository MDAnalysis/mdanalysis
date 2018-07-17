'''
run with: python ten2eleven.py -f agsetterpluralization test_dummy_old_MDA_code.py 
Author: Tyler Reddy
'''
from __future__ import absolute_import

from lib2to3.fixer_base import BaseFix
from lib2to3.fixer_util import Name, Call, LParen, RParen, ArgList, Dot
from lib2to3 import pytree


class FixAgsetterpluralization(BaseFix):

    PATTERN = """
        power< head =any+
                trailer< dot = '.' method=('set_mass'|'set_charge'|'set_name'|
                                            'set_type'|'set_radius'|'set_bfactor'|
                                            'set_altloc'|'set_serial'|'set_resid'|
                                            'set_resname'|'set_resnum'|'set_segid')>
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
        if method_name == 'set_radius': #different plural conversion 
            method_name = 'set_radii'
        elif method_name == 'set_mass': #another different plural form
            method_name = 'set_masses'
        else:
            method_name += 's' #standard plural conversion for all others
        args = head + [pytree.Node(syms.trailer, [Dot(), Name(method_name, prefix = method.prefix)])] + tail
        new = pytree.Node(syms.power, args)
        return new
