'''
run with: python ten2eleven.py -f selectatoms test_dummy_old_MDA_code.py 
Author: Tyler Reddy
'''

from lib2to3.fixer_base import BaseFix
from lib2to3.pgen2 import token


class FixSelectatoms(BaseFix):

    _accept_type = token.NAME

    def match(self, node):
        if node.value == 'selectAtoms':
            return True
        return False

    def transform(self, node, results):
        node.value = 'select_atoms'
        node.changed()

