'''
run with: python ten2eleven.py -f numframes test_dummy_old_MDA_code.py 
Author: Tyler Reddy
'''

from lib2to3.fixer_base import BaseFix
from lib2to3.fixer_util import Name, Dot
from lib2to3 import pytree

class FixNumframes(BaseFix):

    PATTERN = """
                trailer< dot = '.' method='numframes'>
    """

    def transform(self, node, results):
        method = results['method']
        method_name = method.value
        if method_name == 'numframes':
            method_name = 'n_frames'
        syms = self.syms
        args = [pytree.Node(syms.trailer, [Dot(), Name(method_name)])]
        new = pytree.Node(syms.power, args)
        return new

