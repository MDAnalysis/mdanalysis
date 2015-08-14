'''
run with: python ten2eleven.py -f mdaimports test_dummy_old_MDA_code.py 
Author: Tyler Reddy
'''


from lib2to3.fixer_base import BaseFix
from lib2to3.fixer_util import Name

class FixMdaimports(BaseFix):

    PATTERN = """
import_name< start_statement='import' dotted_name< import_1 = 'MDAnalysis' dot = '.' second_statement='KDTree' > > |
import_name< start_statement='import' dotted_name< import_1 = 'MDAnalysis' dot = '.' second_statement='core' '.' core_module={core_modules} > > |
import_from< start_statement='from' import_1 = 'MDAnalysis' 'import' second_statement='KDTree' > |
import_from< start_statement='from' dotted_name< import_1 = 'MDAnalysis' dot = '.' second_statement='core' > 'import' core_module={core_modules}>
""".format(core_modules = "('transformations'|'util'|'log'|'units'|'distances'|'parallel')")

    def transform(self, node, results):
        start_statement = results['start_statement']
        second_statement = results['second_statement']
        import_1 = results['import_1']

        if second_statement.value == 'core':
            core_module = results['core_module']
            dot = results['dot']
            if core_module[0].value == 'units': #not placed in lib
                second_statement.replace(Name('', prefix = second_statement.prefix))
                dot.replace(Name('', prefix = second_statement.prefix))
            else:
                second_statement.replace(Name('lib', prefix = second_statement.prefix))
        elif second_statement.value == 'KDTree' and start_statement.value == 'import':
            second_statement.replace(Name('lib.KDTree', prefix = second_statement.prefix))
        elif second_statement.value == 'KDTree' and start_statement.value == 'from':
            import_1.replace(Name(' MDAnalysis.lib'))
