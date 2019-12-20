#!/usr/bin/env python

from __future__ import print_function
import sys
from pprint import pprint
import collections
import os
import textwrap
import re
import inspect
import tabulate

from sphinx.ext.napoleon import NumpyDocstring

# Make sure we use the same version of MDAnalysis as sphinx
sys.path.insert(0, os.path.abspath('../../..'))
import MDAnalysis as mda

INDENT = ' ' * 8

def clear_citations(doc):
    citation_re = re.compile(r'^ *\.\. \[[^]]+\]')

    result = []
    in_citation = False
    for line in doc.splitlines():
        match = citation_re.match(line)
        if match is not None:
            in_citation = True
        elif in_citation and not line.strip():
            in_citation = False
        elif not in_citation:
            result.append(line)

    return '\n'.join(result)


def clean_doc(doc):
    """Remove whitespace and citations"""
    if not doc:
        return ''
    dedent_doc = textwrap.dedent(INDENT + doc)
    numpy_doc = NumpyDocstring(dedent_doc)
    # doc_clear = clear_citations(str(numpy_doc))
    doc_clear = str(numpy_doc)
    return doc_clear


class TransplantedMethod:

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        return self.name == other.name
    
    def __lt__(self, other):
        return self.name < other.name

    def __init__(self, name, method):
        self.method = method
        self.name = name

        try:
            # We may be dealing with a property; then we need to get the
            # actual method out of it.
            self.method = method.fget
            self.is_property = True
        except AttributeError:
            # Well, it was not a property
            self.is_property = False

        self.is_private = self.name.startswith('_') or self.name.endswith('_')

        self.doc = clean_doc(self.method.__doc__)
        indented = textwrap.indent(self.doc, prefix=INDENT)

        raw_sig = str(inspect.signature(self.method))
        no_group = re.sub(r'\(group *,? *', '(', raw_sig)
        signature = re.sub(r'\(self *,? *', '(', no_group)

        try:
            self.short_desc = self.doc.split('.\n')[0].strip() + '.'
        except IndexError:
            self.short_desc = ''

        if self.is_property:
            self.short_fmt = ':attr:`{}`'.format(self.name)
            self.formatted = '.. method:: {}{}\n    :property:\n\n{}\n\n'.format(self.name,
                                                                signature,
                                                                indented)
        else:
            self.short_fmt = ':meth:`{}`'.format(self.name)
            self.formatted = '.. method:: {}{}\n\n{}\n\n'.format(self.name,
                                                                signature,
                                                                indented)

TRANSPLANT_EXPLANATION = """
Some methods and properties in this class only exist when certain topology attributes are available. They are first listed with the required topology attribute below, and then documented under that.
"""

class GroupTable:

    REQUIRES = '**Requires {}**'

    def __init__(self, group_target, attrs):
        self.attrs = {k:sorted(list(set(v))) for k, v in attrs.items()}

        try:
            name = group_target.__name__
        except AttributeError:
            # For some reason, some target are not classes but str
            name = group_target
        
        file = '{}_methods_table.txt'.format(name)
        file2 = '{}_methods_docs.txt'.format(name)
        self.filename = os.path.join('documentation_pages', 'core', file)
        self.filename2 = os.path.join('documentation_pages', 'core', file2)
        self.table = []
        self.all_methods = []

    def make_table(self):
        for attrname, methods in self.attrs.items():
            self.table.append([self.REQUIRES.format(attrname), ''])
            
            for method in methods:
                self.table.append([method.short_fmt, method.short_desc])
        
            self.all_methods.extend(methods)
                
    
    def write(self):
        with open(self.filename, 'w') as f:

            print(TRANSPLANT_EXPLANATION, file=f)

            print(tabulate.tabulate(self.table, tablefmt='rst'), file=f)
        
            print(file=f)

        with open(self.filename2, 'w') as f:
            for method in sorted(self.all_methods):
                print(method.formatted, file=f)
        
        print('Wrote ', self.filename)



# {Group name: {attribute name: [method]}}
TRANSPLANTS = collections.defaultdict(lambda: collections.defaultdict(list))

# collect transplanted methods
for attr in mda._TOPOLOGY_ATTRS.values():
    for group_target, methods in attr.transplants.items():
        for name, method in methods:
            fn = TransplantedMethod(name, method)
            if not fn.is_private:
                TRANSPLANTS[group_target][attr.attrname].append(fn)

# copy parent methods into subclasses
for group_target, attrmethods in TRANSPLANTS.items():
    try:
        for parent in group_target.__mro__:
            for attrname, parent_methods in TRANSPLANTS.get(parent, {}).items():
                for method in parent_methods:
                    if method not in attrmethods[attrname]:
                        attrmethods[attrname].append(method)
    except AttributeError:
        pass

# write docs
for group_target, attrs in TRANSPLANTS.items():
    tablewriter = GroupTable(group_target, attrs)
    tablewriter.make_table()
    tablewriter.write()
