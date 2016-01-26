'''
run with: python ten2eleven.py -f writer test_dummy_old_MDA_code.py 
Author: Tyler Reddy
'''

from lib2to3.fixer_base import BaseFix
from lib2to3.fixer_util import Name

class FixWriter(BaseFix):

    PATTERN = """
                power< any+ trailer< '.' method = ('Writer' | 'writer') > trailer< '(' arglist< any+ argument< argname = 'numatoms' '=' any*> any* > ')' > >
    """

    def transform(self, node, results):
        argname = results['argname']
        argname.replace(Name('n_atoms', prefix=argname.prefix))
