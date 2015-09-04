'''
run with: python ten2eleven.py -f framenumberingzerobased test_dummy_old_MDA_code.py 
Author: Tyler Reddy
'''

from lib2to3.fixer_base import BaseFix
from lib2to3.fixer_util import Name, Call, LParen, RParen, ArgList, Dot
from lib2to3 import pytree


class FixFramenumberingzerobased(BaseFix):

    PATTERN = """
    power< head = any+
    trailer< dot = '.' method = 'frame' >
    tail = any*>
    """

    def transform(self, node, results):
        #add a warning comment for 0-based frame numbering rather than attempting to fix; this is a rather complicated issue for which user intervention is preferable
        head = results['head'][0]
        comment_string = '\n#ten2eleven.py detected a possible incompatibility between this code and MDAnalysis >= 0.11.0\n#Frame numbering is now 0-based\n#Please manually review the following lines (and remove these comments afterwards):\n'
        head.prefix = comment_string + head.prefix 
