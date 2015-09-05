'''
run with: python ten2eleven.py -f framenumberingzerobased test_dummy_old_MDA_code.py 
Author: Tyler Reddy
'''

from lib2to3.fixer_base import BaseFix
from lib2to3.fixer_util import Name, Call, LParen, RParen, ArgList, Dot
from lib2to3 import pytree


class FixFramenumberingzerobased(BaseFix):

    PATTERN = """
    power<head=any+
    trailer< dot = '.' method = 'frame' >
    tail=any*>
    """

    def transform(self, node, results):
        #add a warning comment for 0-based frame numbering rather than attempting to fix; this is a rather complicated issue for which user intervention is preferable
        syms = self.syms
        method = results['method']
        method_name = method.value
        head = results['head']
        head = [n.clone() for n in head]
        comment_string = '\n#ten2eleven.py detected a possible incompatibility between this code and MDAnalysis >= 0.11.0\n#Frame numbering is now 0-based\n#Please manually review the following lines (and remove these comments afterwards):\n'
        #args = [Name('#dummy')] + head + [pytree.Node(syms.trailer, [Dot(), Name(method_name, prefix = method.prefix), Dot()])]
        #print 'node:', node
        if node.prefix:
            node.set_prefix(comment_string + node.prefix)
            node.changed()
        else: #more problematic cases treated in another fixer
            pass
