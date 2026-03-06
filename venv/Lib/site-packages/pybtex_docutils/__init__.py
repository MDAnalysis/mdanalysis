"""
API
~~~

The backend renders :class:`pybtex.richtext.Text` instances
into a list of :class:`docutils.nodes.Node` instances.
For typical use cases, all you need to care about are the methods
:meth:`Backend.paragraph`,
:meth:`Backend.citation`, and
:meth:`Backend.citation_reference`
which are to be called on *formatted* entries,
as in the :ref:`minimal example <minimal-example>`.

Unless you are subclassing :class:`Backend` to create a new backend,
you should normally not import the :mod:`pybtex_docutils` module directly.
Instead, use pybtex's plugin system to get the :class:`Backend` class,
as in the :ref:`minimal example <minimal-example>`.

.. autoclass:: Backend
   :show-inheritance:
   :members: RenderType, paragraph, citation, citation_reference,
             footnote, footnote_reference

.. autoclass:: SimpleBibliography
   :show-inheritance:
   :members: run
"""

import docutils.nodes
import itertools

import docutils.parsers.rst.directives as directives
from docutils.parsers.rst import Directive
from pybtex.backends import BaseBackend
import pybtex.database
import pybtex.plugin
import os.path
from typing import TYPE_CHECKING, List, Type

from pybtex.database.input.bibtex import Parser

if TYPE_CHECKING:
    from pybtex.style import FormattedEntry


class Backend(BaseBackend):
    name = 'docutils'

    symbols = {
        'ndash': [docutils.nodes.Text('\u2013')],
        'newblock': [docutils.nodes.Text(' ')],
        'nbsp': [docutils.nodes.Text('\u00a0')],
    }
    tags = {
        'emph': docutils.nodes.emphasis,  # note: deprecated
        'em': docutils.nodes.emphasis,
        'strong': docutils.nodes.strong,
        'i': docutils.nodes.emphasis,
        'b': docutils.nodes.strong,
        'tt': docutils.nodes.literal,
        'sup': docutils.nodes.superscript,
        'sub': docutils.nodes.subscript,
    }

    RenderType: Type[List[docutils.nodes.Node]] = list

    # for compatibility only
    def format_text(self, text: str) -> List[docutils.nodes.Node]:
        return self.format_str(text)

    def format_str(self, str_: str) -> List[docutils.nodes.Node]:
        return [docutils.nodes.Text(str_)]

    def format_tag(self, tag_name: str, text: List[docutils.nodes.Node]
                   ) -> List[docutils.nodes.Node]:
        if tag_name in self.tags:
            tag = self.tags[tag_name]
            return [tag('', '', *text)]
        else:
            return text

    def format_href(self, url: str, text: List[docutils.nodes.Node],
                    external: bool = False) -> List[docutils.nodes.Node]:
        node = docutils.nodes.reference('', '', *text, refuri=url)
        return [node]

    def write_entry(self, key: str, label: str, text: str) -> None:
        raise NotImplementedError("use Backend.citation() instead")

    def render_sequence(self, rendered_list: List[List[docutils.nodes.Node]]
                        ) -> List[docutils.nodes.Node]:
        return list(itertools.chain(*rendered_list))

    def paragraph(self, entry: "FormattedEntry") -> docutils.nodes.paragraph:
        """Return a docutils.nodes.paragraph
        containing the rendered text for *entry* (without label).

        .. versionadded:: 0.2.0
        """
        return docutils.nodes.paragraph('', '', *entry.text.render(self))

    def citation(self, entry: "FormattedEntry",
                 document: docutils.nodes.document, use_key_as_label=True
                 ) -> docutils.nodes.citation:
        """Return citation node, with key as name, label as first
        child, and paragraph with entry text as second child. The citation is
        expected to be inserted into *document* prior to any docutils
        transforms.
        """
        # see docutils.parsers.rst.states.Body.citation()
        if use_key_as_label:
            label = entry.key
        else:
            label = entry.label
        name = docutils.nodes.fully_normalize_name(entry.key)
        citation = docutils.nodes.citation()
        citation['names'].append(name)
        citation += docutils.nodes.label('', label)
        citation += self.paragraph(entry)
        document.note_citation(citation)
        document.note_explicit_target(citation, citation)
        return citation

    def citation_reference(
            self, entry: "FormattedEntry",
            document: docutils.nodes.document, use_key_as_label=True
            ) -> docutils.nodes.citation_reference:
        """Return citation_reference node to the given citation. The
        citation_reference is expected to be inserted into *document*
        prior to any docutils transforms.
        """
        # see docutils.parsers.rst.states.Body.footnote_reference()
        if use_key_as_label:
            label = entry.key
        else:
            label = entry.label
        refname = docutils.nodes.fully_normalize_name(entry.key)
        refnode = docutils.nodes.citation_reference(
            '[%s]_' % label, refname=refname)
        refnode += docutils.nodes.Text(label)
        document.note_citation_ref(refnode)
        return refnode

    def footnote(self, entry: "FormattedEntry",
                 document: docutils.nodes.document) -> docutils.nodes.footnote:
        """Return footnote node, with key as name, and paragraph with
        entry text as child. The footnote is expected to be
        inserted into *document* prior to any docutils transforms.

        .. versionadded:: 0.2.2
        """
        # see docutils.parsers.rst.states.Body.footnote()
        name = docutils.nodes.fully_normalize_name(entry.key)
        footnote = docutils.nodes.footnote(auto=1)
        footnote['names'].append(name)
        footnote += self.paragraph(entry)
        document.note_autofootnote(footnote)
        document.note_explicit_target(footnote, footnote)
        return footnote

    def footnote_reference(self, entry: "FormattedEntry",
                           document: docutils.nodes.document
                           ) -> docutils.nodes.footnote_reference:
        """Return footnote_reference node to the given citation. The
        footnote_reference is expected to be inserted into *document*
        prior to any docutils transforms.

        .. versionadded:: 0.2.2
        """
        # see docutils.parsers.rst.states.Body.footnote_reference()
        refname = docutils.nodes.fully_normalize_name(entry.key)
        refnode = docutils.nodes.footnote_reference(
            '[#%s]_' % entry.key, refname=refname, auto=1)
        document.note_autofootnote_ref(refnode)
        document.note_footnote_ref(refnode)
        return refnode


class SimpleBibliography(Directive):
    name = "simplebibliography"
    required_arguments = 1
    optional_arguments = 0
    final_argument_whitespace = True
    has_content = False
    option_spec = {
        'encoding': directives.unchanged,
        'style': directives.unchanged,
    }

    def run(self):
        parser = Parser(self.options.get("encoding", "utf-8-sig"))
        for filename_raw in self.arguments[0].split():
            filename = os.path.join(
                os.path.dirname(self.state_machine.document['source']),
                filename_raw)
            if not os.path.isfile(filename):
                raise self.error(f"could not open bibtex file {filename}")
            else:
                try:
                    parser.parse_file(filename)
                except pybtex.database.BibliographyDataError as exc:
                    raise self.error(
                        f"bibliography data error in {filename}: {exc}")
        style = pybtex.plugin.find_plugin(
            "pybtex.style.formatting", self.options.get("style", "unsrt"))()
        backend = Backend()
        document = self.state_machine.document
        return [
            backend.citation(style.format_entry(entry.key, entry), document)
            for entry in style.sort(parser.data.entries.values())]
