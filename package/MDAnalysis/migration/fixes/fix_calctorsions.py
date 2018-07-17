'''
run with: python ten2eleven.py -f calctorsions test_dummy_old_MDA_code.py 
Author: Tyler Reddy
'''
from __future__ import absolute_import

from lib2to3.fixer_base import BaseFix
from lib2to3.fixer_util import Name, Call, LParen, RParen, ArgList, Dot
from lib2to3 import pytree


class FixCalctorsions(BaseFix):

    PATTERN = """
                trailer< '.' method='calc_torsions' >
                |
                import_from< 'from' dotted_name< 'MDAnalysis' '.' 'lib' '.' 'distances' > 'import' import_name='calc_torsions' >
    """

    def transform(self, node, results):
        if 'method' in results.keys():
            method = results['method']
            syms = self.syms
            method_name = method.value
            replacement_dict = {'calc_torsions': 'calc_dihedrals'}
            method_name = replacement_dict[method_name]
            args = [pytree.Node(syms.trailer, [Dot(), Name(method_name, prefix = method.prefix)])]
            new = pytree.Node(syms.power, args)
            return new
        elif 'import_name' in results.keys():
            import_name = results['import_name']
            import_name.replace(Name(' calc_dihedrals'))
