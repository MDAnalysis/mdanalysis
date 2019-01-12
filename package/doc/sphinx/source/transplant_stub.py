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

class TransplantedMethod:
    def __init__(self, method):
        self.method = method
        try:
            # We may be dealing with a property; then we need to get the
            # actual method out of it.
            self.method = method.fget
        except AttributeError:
            # Well, it was not a property
            pass

    @property
    def name(self):
        return self.method.__name__

    @property
    def doc(self):
        dedent_doc = textwrap.dedent('        ' + self.method.__doc__)
        numpy_doc = NumpyDocstring(dedent_doc)
        doc_clear = clear_citations(str(numpy_doc))
        return doc_clear

    @property
    def signature(self):
        return get_signature(self.method)

    @property
    def short_desc(self):
        return self.doc.splitlines()[0].strip()

    @property
    def is_private(self):
        return self.name.startswith('_') or self.name.endswith('_')

    @property
    def formatted(self):
        text = '.. method:: {}{}\n\n{}\n\n'.format(
            self.name,
            self.signature,
            textwrap.indent(self.doc, prefix=' ' * 8)
        )
        return text


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


def get_signature(method):
    signature = str(inspect.signature(method))
    return re.sub(r'\(self,? *', '(', signature)


# Collect the transplanted functions from the topopoly attributes
targets = collections.defaultdict(lambda : collections.defaultdict(list))
for attribute_key, attribute in mda.core.topologyattrs._TOPOLOGY_ATTRS.items():
    for target, methods in attribute.transplants.items():
        all_methods = []
        for method in methods:
            function = TransplantedMethod(method[1])
            if not function.is_private:
                all_methods.append(function)
        if all_methods:
            targets[target][attribute.attrname] = all_methods


for target_key, target_dict in targets.items():
    try:
        target_name = target_key.__name__
    except AttributeError:
        # For some reason, some target are not classes but str
        target_name = target_key
    if hasattr(target_key, '__mro__'):
        for parent in target_key.__mro__:
            for attribute_key, method_list in targets.get(parent, {}).items():
                if attribute_key not in target_dict:
                    target_dict[attribute_key] = []
                for method in method_list:
                    if method not in target_dict[attribute_key]:
                        target_dict[attribute_key].append(method)


for target_key, target_dict in targets.items():
    try:
        target_name = target_key.__name__
    except AttributeError:
        # For some reason, some target are not classes but str
        target_name = target_key
    table = []
    for attribute_key, method_list in target_dict.items():
        table.append(['**Requires {}**'.format(attribute_key), ''])
        for method in method_list:
            table.append([method.name, method.short_desc])
    print(tabulate.tabulate(table, tablefmt='grid'))
        
    
for target_key, target_dict in targets.items():
    try:
        target_name = target_key.__name__
    except AttributeError:
        # For some reason, some target are not classes but str
        target_name = target_key
    file_name = os.path.join(
        'documentation_pages',
        'core',
        '{}.txt'.format(target_name)
    )
    with open(file_name, 'w') as outfile:
        table = []
        for attribute_key, method_list in target_dict.items():
            table.append(['**Requires {}**'.format(attribute_key), ''])
            for method in method_list:
                table.append([':meth:`{}`'.format(method.name), method.short_desc])
        print(tabulate.tabulate(table, tablefmt='grid'), file=outfile)

        for attribute_key, method_list in target_dict.items():
            print(file=outfile)

            for method in method_list:
                print(method.formatted, file=outfile)


