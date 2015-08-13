'''
Copy to /Users/treddy/anaconda/lib/python2.7/lib2to3/fixes/
Then run with: 2to3 -f selectatoms test_dummy_old_MDA_code.py 
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

