'''
run with: python ten2eleven.py -f torsionclasses test_dummy_old_MDA_code.py 
Author: Tyler Reddy (but effectively a hack of lib2to3/fixes/fix_metaclass.py)
'''
from __future__ import absolute_import

from lib2to3.fixer_base import BaseFix
from lib2to3.fixer_util import Name, Dot, syms, Node, Leaf
from lib2to3 import pytree
from lib2to3.fixes import fix_metaclass
from lib2to3.fixes.fix_metaclass import has_metaclass, fixup_parse_tree, find_metas, fixup_indent
from lib2to3.pygram import token


class FixTorsionclasses(fix_metaclass.FixMetaclass):
    def transform(self, node, results):
        fixup_parse_tree(node)

        text_type = node.children[0].type # always Leaf(nnn, 'class')

        # figure out what kind of classdef we have
        if len(node.children) == 7:
            # Node(classdef, ['class', 'name', '(', arglist, ')', ':', suite])
            #                 0        1       2    3        4    5    6
            if node.children[3].type == syms.arglist:
                arglist = node.children[3]
            # Node(classdef, ['class', 'name', '(', 'Parent', ')', ':', suite])
            else:
                parent = node.children[3].clone()
                arglist = Node(syms.arglist, [parent])
                node.set_child(3, arglist)
        elif len(node.children) == 6:
            # Node(classdef, ['class', 'name', '(',  ')', ':', suite])
            #                 0        1       2     3    4    5
            arglist = Node(syms.arglist, [])
            node.insert_child(3, arglist)
        elif len(node.children) == 4:
            # Node(classdef, ['class', 'name', ':', suite])
            #                 0        1       2    3
            arglist = Node(syms.arglist, [])
            node.insert_child(2, Leaf(token.RPAR, u')'))
            node.insert_child(2, arglist)
            node.insert_child(2, Leaf(token.LPAR, u'('))
        else:
            raise ValueError("Unexpected class definition")

        for child in arglist.children:
            if child.value == 'Torsion':
                child.value = 'Dihedral'
            if child.value == 'Improper_Torsion':
                child.value = 'ImproperDihedral'
