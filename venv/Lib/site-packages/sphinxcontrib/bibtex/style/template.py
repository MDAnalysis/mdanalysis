"""
.. autofunction:: sphinxcontrib.bibtex.style.template.join(\
        sep='', sep2=None, last_sep=None, other=None)

.. autofunction:: sphinxcontrib.bibtex.style.template.sentence(\
        capfirst=False, capitalize=False, add_period=True, \
        sep=', ', sep2=None, last_sep=None, other=None)

.. autofunction:: sphinxcontrib.bibtex.style.template.names(\
        role, sep='', sep2=None, last_sep=None, other=None)

.. autofunction:: sphinxcontrib.bibtex.style.template.entry_label()

.. autofunction:: sphinxcontrib.bibtex.style.template.reference()

.. autofunction:: sphinxcontrib.bibtex.style.template.footnote_reference()
"""

import docutils.nodes
import pybtex_docutils
from pybtex.richtext import Text
from pybtex.style.template import (
    Node, _format_list, FieldIsMissing, field, first_of, optional, tag
)
from sphinx.util.nodes import make_refnode

from sphinxcontrib.bibtex.nodes import raw_latex
from sphinxcontrib.bibtex.richtext import BaseReferenceText

from typing import TYPE_CHECKING, Dict, Any, cast, NamedTuple, List

if TYPE_CHECKING:
    from pybtex.backends import BaseBackend
    from pybtex.richtext import BaseText
    from pybtex.style import FormattedEntry
    from sphinx.builders import Builder


# extended from pybtex: also copies the docstring into the wrapped object
def node(f):
    n = Node(f.__name__, f)
    n.__doc__ = f.__doc__
    return n


# copied from pybtex join but extended to allow "et al" formatting
@node
def join(children, data, sep='', sep2=None, last_sep=None, other=None):
    """Join text fragments together."""

    if sep2 is None:
        sep2 = sep
    if last_sep is None:
        last_sep = sep
    parts = [part for part in _format_list(children, data) if part]
    if len(parts) <= 1:
        return Text(*parts)
    elif len(parts) == 2:
        return Text(sep2).join(parts)
    elif other is None:
        return Text(last_sep).join([Text(sep).join(parts[:-1]), parts[-1]])
    else:
        return Text(parts[0], other)


# copied from pybtex names but using the new join
@node
def sentence(children, data, capfirst=False, capitalize=False, add_period=True,
             sep=', ', sep2=None, last_sep=None, other=None):
    """Join text fragments, capitalize the first letter,
    and add a period to the end.
    """
    text = join(sep=sep, sep2=sep2, last_sep=last_sep, other=other)[
        children
    ].format_data(data)
    if capfirst:
        text = text.capfirst()
    if capitalize:
        text = text.capitalize()
    if add_period:
        text = text.add_period()
    return text


# copied from pybtex names but using the new join allowing "et al" formatting
@node
def names(children, data, role, **kwargs):
    """Return formatted names."""
    assert not children
    try:
        persons = data['entry'].persons[role]
    except KeyError:
        raise FieldIsMissing(role, data['entry'])
    style = data['style']
    formatted_names = [
        style.person.style_plugin.format(person, style.person.abbreviate)
        for person in persons]
    return join(**kwargs)[formatted_names].format_data(data)


@node
def entry_label(children, data) -> "BaseText":
    """Node for inserting the label of a formatted entry."""
    assert not children
    entry = cast("FormattedEntry", data['formatted_entry'])
    return Text(entry.label)


class SphinxReferenceInfo(NamedTuple):
    """Tuple containing reference info to enable sphinx to resolve a reference
    to a citation.
    """
    builder: "Builder"  #: The Sphinx builder.
    fromdocname: str    #: Document name of the citation reference.
    todocname: str      #: Document name of the bibliography.
    citation_id: str    #: Unique id of the citation within the bibliography.
    title: str          #: Title attribute for reference node.


class SphinxReferenceText(BaseReferenceText[SphinxReferenceInfo]):
    """Pybtex rich text class generating
    a docutils reference node to a citation
    for use with :class:`SphinxReferenceInfo`.
    """

    def render(self, backend: "BaseBackend") -> List[docutils.nodes.Element]:
        assert isinstance(backend, pybtex_docutils.Backend), \
               "SphinxReferenceText only supports the docutils backend"
        info = self.info[0]
        if info.builder.name == 'latex':
            key = f'cite.{info.todocname}:{info.citation_id}'
            return (
                [raw_latex(f'\\hyperlink{{{key}}}{{')]
                + super().render(backend)
                + [raw_latex('}')]
            )
        elif info.builder.name == 'rinoh':
            children = super().render(backend)
            refid = f"%{info.todocname}#{info.citation_id}"
            refnode = docutils.nodes.citation_reference(
                text=children[0], refid=refid, reftitle=info.title)
            refnode.extend(children[1:])
            return [refnode]
        else:
            children = super().render(backend)
            # make_refnode only takes a single child
            refnode2 = make_refnode(
                builder=info.builder,
                fromdocname=info.fromdocname,
                todocname=info.todocname,
                targetid=info.citation_id,
                child=children[0],
                title=info.title,
            )
            refnode2.extend(children[1:])  # type: ignore
            return [refnode2]


@node
def reference(children, data: Dict[str, Any]):
    """Pybtex node for inserting a docutils reference node to a citation.
    The children of the node
    comprise the content of the reference, and any referencing information
    is stored in the *reference_info* key of the *data*.
    The data must also contain a *style* key pointing to the corresponding
    :class:`~sphinxcontrib.bibtex.style.referencing.BaseReferenceStyle`.
    """
    parts = _format_list(children, data)
    info = data['reference_info']
    assert isinstance(info, SphinxReferenceInfo)
    return SphinxReferenceText(info, *parts)


class FootReferenceInfo(NamedTuple):
    """Tuple containing reference info to enable sphinx to resolve a footnote
    reference.
    """
    key: str                             #: Citation key.
    document: "docutils.nodes.document"  #: Current docutils document.
    refname: str                         #: Citation reference name.


class FootReferenceText(BaseReferenceText[FootReferenceInfo]):
    """Pybtex rich text class generating
    a docutils footnote_reference node to a citation
    for use with :class:`FootReferenceInfo`.
    """

    def render(self, backend: "BaseBackend"):
        assert isinstance(backend, pybtex_docutils.Backend), \
               "FootReferenceText only supports the docutils backend"
        info = self.info[0]
        # see docutils.parsers.rst.states.Body.footnote_reference()
        refnode = docutils.nodes.footnote_reference(
            '[#%s]_' % info.key, refname=info.refname, auto=1)
        info.document.note_autofootnote_ref(refnode)
        info.document.note_footnote_ref(refnode)
        return [refnode]


@node
def footnote_reference(children, data: Dict[str, Any]):
    """Pybtex node for inserting a footnote_reference docutils node.
    Any referencing information
    is stored in the *reference_info* key of the *data*.
    The data must also contain a *style* key pointing to the corresponding
    :class:`~sphinxcontrib.bibtex.style.referencing.BaseReferenceStyle`.
    """
    assert not children
    info = data['reference_info']
    assert isinstance(info, FootReferenceInfo)
    # we need to give the footnote text some fake content
    # otherwise pybtex richtext engine will mess things up
    return FootReferenceText(info, '#')


@node
def year(children, data: Dict[str, Any]) -> "BaseText":
    assert not children
    return first_of[optional[field('year')], 'n.d.'].format_data(data)


@node
def author_or_editor_or_title(children, data, **kwargs):
    assert not children
    return first_of[
        optional[names('author', **kwargs)],
        optional[names('editor', **kwargs)],
        tag('em')[field('title')]].format_data(data)
