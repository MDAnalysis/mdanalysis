"""
    Classes and methods to maintain any bibtex information that is stored
    outside the doctree.

    .. autoclass:: Citation
        :members:

    .. autoclass:: BibtexDomain
        :members:
"""

import ast
from typing import TYPE_CHECKING
from typing import List, Dict, NamedTuple, cast, Iterable, Tuple, Set, Optional

import docutils.frontend
import docutils.nodes
import docutils.parsers.rst
import docutils.utils
import pybtex_docutils
import pybtex.plugin
from pybtex.richtext import Tag
from pybtex.style.template import FieldIsMissing
from pybtex.style import FormattedEntry

import sphinxcontrib.bibtex.plugin
import sphinx.util
import re

from sphinx.domains import Domain, ObjType
from sphinx.errors import ExtensionError
from sphinx.locale import _

from .roles import CiteRole
from .bibfile import normpath_filename, process_bibdata, BibData
from .style.referencing import BaseReferenceStyle, format_references
from .style.template import SphinxReferenceInfo

if TYPE_CHECKING:
    from pybtex.database import Entry
    from pybtex.style.formatting import BaseStyle
    from sphinx.addnodes import pending_xref
    from sphinx.application import Sphinx
    from sphinx.builders import Builder
    from sphinx.environment import BuildEnvironment
    from .directives import BibliographyKey, BibliographyValue
    from .roles import CitationRef

logger = sphinx.util.logging.getLogger(__name__)


def _raise_invalid_node(node):
    """Helper method to raise an exception when an invalid node is
    visited.
    """
    raise ValueError("invalid node %s in filter expression" % node)


class _FilterVisitor(ast.NodeVisitor):

    """Visit the abstract syntax tree of a parsed filter expression."""

    entry = None
    """The bibliographic entry to which the filter must be applied."""

    cited_docnames = False
    """The documents where the entry is cited (empty if not cited)."""

    def __init__(self, entry, docname, cited_docnames):
        self.entry = entry
        self.docname = docname
        self.cited_docnames = cited_docnames

    def visit_Module(self, node):
        if len(node.body) != 1:
            raise ValueError(
                "filter expression cannot contain multiple expressions")
        return self.visit(node.body[0])

    def visit_Expr(self, node):
        return self.visit(node.value)

    def visit_BoolOp(self, node):
        outcomes = (self.visit(value) for value in node.values)
        if isinstance(node.op, ast.And):
            return all(outcomes)
        elif isinstance(node.op, ast.Or):
            return any(outcomes)
        else:  # pragma: no cover
            # there are no other boolean operators
            # so this code should never execute
            assert False, "unexpected boolean operator %s" % node.op

    def visit_UnaryOp(self, node):
        if isinstance(node.op, ast.Not):
            return not self.visit(node.operand)
        else:
            _raise_invalid_node(node)

    def visit_BinOp(self, node):
        left = self.visit(node.left)
        op = node.op
        right = self.visit(node.right)
        if isinstance(op, ast.Mod):
            # modulo operator is used for regular expression matching
            if not isinstance(left, str):
                raise ValueError(
                    "expected a string on left side of %s" % node.op)
            if not isinstance(right, str):
                raise ValueError(
                    "expected a string on right side of %s" % node.op)
            return re.search(right, left, re.IGNORECASE)
        elif isinstance(op, ast.BitOr):
            return left | right
        elif isinstance(op, ast.BitAnd):
            return left & right
        else:
            _raise_invalid_node(node)

    def visit_Compare(self, node):
        # keep it simple: binary comparators only
        if len(node.ops) != 1:
            raise ValueError("syntax for multiple comparators not supported")
        left = self.visit(node.left)
        op = node.ops[0]
        right = self.visit(node.comparators[0])
        if isinstance(op, ast.Eq):
            return left == right
        elif isinstance(op, ast.NotEq):
            return left != right
        elif isinstance(op, ast.Lt):
            return left < right
        elif isinstance(op, ast.LtE):
            return left <= right
        elif isinstance(op, ast.Gt):
            return left > right
        elif isinstance(op, ast.GtE):
            return left >= right
        elif isinstance(op, ast.In):
            return left in right
        elif isinstance(op, ast.NotIn):
            return left not in right
        else:
            # not used currently: ast.Is | ast.IsNot
            _raise_invalid_node(op)

    def visit_Name(self, node):
        """Calculate the value of the given identifier."""
        id_ = node.id
        if id_ == 'type':
            return self.entry.type.lower()
        elif id_ == 'key':
            return self.entry.key.lower()
        elif id_ == 'cited':
            return bool(self.cited_docnames)
        elif id_ == 'docname':
            return self.docname
        elif id_ == 'docnames':
            return self.cited_docnames
        elif id_ == 'author' or id_ == 'editor':
            if id_ in self.entry.persons:
                return u' and '.join(
                    str(person)  # XXX needs fix in pybtex?
                    for person in self.entry.persons[id_])
            else:
                return u''
        else:
            return self.entry.fields.get(id_, "")

    def visit_Set(self, node):
        return frozenset(self.visit(elt) for elt in node.elts)

    # NameConstant is Python 3.4 only
    def visit_NameConstant(self, node):
        return node.value  # pragma: no cover

    # Constant is Python 3.6+ only
    # Since 3.8 Num, Str, Bytes, NameConstant and Ellipsis are just Constant
    def visit_Constant(self, node):
        return node.value

    # Not used on 3.8+
    def visit_Str(self, node):
        return node.s  # pragma: no cover

    def generic_visit(self, node):
        _raise_invalid_node(node)


def get_docnames(env):
    """Get document names in order."""
    rel = env.collect_relations()
    docname = (
        env.config.master_doc if sphinx.version_info < (4, 0)
        else env.config.root_doc)
    docnames = set()
    while docname is not None:
        docnames.add(docname)
        yield docname
        parent, prevdoc, nextdoc = rel[docname]
        docname = nextdoc
    for docname in sorted(env.found_docs - docnames):
        yield docname


class Citation(NamedTuple):
    """Information about a citation."""
    citation_id: str                     #: Unique id of this citation.
    bibliography_key: "BibliographyKey"  #: Key of its bibliography directive.
    key: str                             #: Key (with prefix).
    entry: "Entry"                       #: Entry from pybtex.
    formatted_entry: "FormattedEntry"    #: Formatted entry for bibliography.
    tooltip_entry: Optional["FormattedEntry"]  #: Formatted entry for tooltip.


def env_updated(app: "Sphinx", env: "BuildEnvironment") -> Iterable[str]:
    domain = cast(BibtexDomain, env.get_domain('cite'))
    return domain.env_updated()


def parse_header(header: str, source_path: str):
    parser = docutils.parsers.rst.Parser()
    # note: types stub for docutils doesn't know about components argument
    settings = docutils.frontend.OptionParser(
        components=(docutils.parsers.rst.Parser,)  # type: ignore
    ).get_default_values()
    document = docutils.utils.new_document(source_path, settings)
    parser.parse(header, document)
    return document[0]


class BibtexDomain(Domain):
    """Sphinx domain for the bibtex extension."""

    name = 'cite'
    label = 'BibTeX Citations'
    data_version = 4
    initial_data = dict(
        bibdata=BibData(
            encoding='',
            bibfiles={},
            data=pybtex.database.BibliographyData()),
        bibliography_header=docutils.nodes.container(),
        bibliographies={},
        citations=[],
        citation_refs=[],
    )
    backend = pybtex_docutils.Backend()
    reference_style: BaseReferenceStyle

    @property
    def bibdata(self) -> BibData:
        """Information about the bibliography files."""
        return self.data['bibdata']

    @property
    def bibliography_header(self) -> docutils.nodes.Element:
        return self.data['bibliography_header']

    @property
    def bibliographies(self) -> Dict["BibliographyKey", "BibliographyValue"]:
        """Map storing information about each bibliography directive."""
        return self.data['bibliographies']

    @property
    def citations(self) -> List[Citation]:
        """Citation data."""
        return self.data['citations']

    @property
    def citation_refs(self) -> List["CitationRef"]:
        """Citation reference data."""
        return self.data['citation_refs']

    def __init__(self, env: "BuildEnvironment"):
        # set up referencing style
        style = sphinxcontrib.bibtex.plugin.find_plugin(
                'sphinxcontrib.bibtex.style.referencing',
                env.app.config.bibtex_reference_style)
        self.reference_style = style()
        # set up object types and roles for referencing style
        role_names = self.reference_style.role_names()
        self.object_types = dict(
            citation=ObjType(_('citation'), *role_names, searchprio=-1),
        )
        self.roles = dict((name, CiteRole()) for name in role_names)
        # initialize the domain
        super().__init__(env)
        # connect env-updated
        env.app.connect('env-updated', env_updated)
        # check config
        if env.app.config.bibtex_bibfiles is None:
            raise ExtensionError(
                "You must configure the bibtex_bibfiles setting")
        # update bib file information in the cache
        bibfiles = [
            normpath_filename(env, "/" + bibfile)
            for bibfile in env.app.config.bibtex_bibfiles]
        self.data['bibdata'] = process_bibdata(
            self.bibdata, bibfiles, env.app.config.bibtex_encoding)
        # parse bibliography header
        header = getattr(env.app.config, "bibtex_bibliography_header")
        if header:
            self.data["bibliography_header"] += \
                parse_header(header, "bibliography_header")

    def clear_doc(self, docname: str) -> None:
        self.data['citations'] = [
            citation for citation in self.citations
            if citation.bibliography_key.docname != docname]
        self.data['citation_refs'] = [
            ref for ref in self.citation_refs if ref.docname != docname]
        for bib_key in list(self.bibliographies.keys()):
            if bib_key.docname == docname:
                del self.bibliographies[bib_key]

    def merge_domaindata(self, docnames: List[str], otherdata: Dict) -> None:
        for bib_key, bib_value in otherdata['bibliographies'].items():
            if bib_key.docname in docnames:
                self.bibliographies[bib_key] = bib_value
        for citation_ref in otherdata['citation_refs']:
            if citation_ref.docname in docnames:
                self.citation_refs.append(citation_ref)
        # 'citations' domain data calculated in env_updated

    def env_updated(self) -> Iterable[str]:
        # This function is called when all doctrees are parsed,
        # but before any post transforms are applied. We use it to
        # determine which citations will be added to which bibliography
        # directive, and also to format the labels. We need to format
        # the labels here because they must be known when resolve_xref is
        # called.
        self.citations.clear()  # might have been restored from pickle
        docnames = list(get_docnames(self.env))
        # we keep track of this to quickly check for duplicates
        used_keys: Set[str] = set()
        used_labels: Dict[str, str] = {}
        for bibliography_key, bibliography in self.bibliographies.items():
            for entry, formatted_entry, tooltip_entry in \
                    self.get_formatted_entries(
                        bibliography_key, docnames,
                        self.env.app.config.bibtex_tooltips,
                        self.env.app.config.bibtex_tooltips_style):
                key = bibliography.keyprefix + formatted_entry.key
                if bibliography.list_ == 'citation' and key in used_keys:
                    logger.warning(
                        'duplicate citation for key "%s"' % key,
                        location=(bibliography_key.docname, bibliography.line),
                        type="bibtex", subtype="duplicate_citation")
                self.citations.append(Citation(
                    citation_id=bibliography.citation_nodes[key]['ids'][0],
                    bibliography_key=bibliography_key,
                    key=key,
                    entry=entry,
                    formatted_entry=formatted_entry,
                    tooltip_entry=tooltip_entry,
                ))
                if bibliography.list_ == 'citation':
                    used_keys.add(key)
                    if formatted_entry.label not in used_labels:
                        used_labels[formatted_entry.label] = key
                    elif used_labels[formatted_entry.label] != key:
                        # if used_label[label] == key then already
                        # duplicate key warning
                        logger.warning(
                            'duplicate label "%s" for keys "%s" and "%s"' % (
                                formatted_entry.label,
                                used_labels[formatted_entry.label], key),
                            location=(bibliography_key.docname,
                                      bibliography.line),
                            type="bibtex", subtype="duplicate_label")
        return []  # expects list of updated docnames

    def resolve_xref(self, env: "BuildEnvironment", fromdocname: str,
                     builder: "Builder", typ: str, target: str,
                     node: "pending_xref", contnode: docutils.nodes.Element
                     ) -> docutils.nodes.Element:
        """Replace node by list of citation references (one for each key)."""
        keys = [key.strip() for key in target.split(',')]
        citations: Dict[str, Citation] = {
            cit.key: cit for cit in self.citations
            if cit.key in keys
            and self.bibliographies[cit.bibliography_key].list_ == 'citation'}
        for key in keys:
            if key not in citations:
                logger.warning('could not find bibtex key "%s"' % key,
                               location=node, type="bibtex",
                               subtype="key_not_found")
        plaintext = pybtex.plugin.find_plugin('pybtex.backends', 'plaintext')()
        references = [
            (citation.entry, citation.formatted_entry, SphinxReferenceInfo(
                builder=builder,
                fromdocname=fromdocname,
                todocname=citation.bibliography_key.docname,
                citation_id=citation.citation_id,
                title=(
                    citation.tooltip_entry.text.render(plaintext).replace(
                        "\\url ", "")
                    if citation.tooltip_entry else None
                )
            ))
            for citation in citations.values()]
        formatted_references = format_references(
            self.reference_style, typ, references)
        result_node = docutils.nodes.inline(rawsource=target)
        result_node += formatted_references.render(self.backend)
        return result_node

    def resolve_any_xref(self, env: "BuildEnvironment", fromdocname: str,
                         builder: "Builder", target: str,
                         node: "pending_xref", contnode: docutils.nodes.Element
                         ) -> List[Tuple[str, docutils.nodes.Element]]:
        """Replace node by list of citation references (one for each key),
        provided that the target has citation keys.
        """
        keys = [key.strip() for key in target.split(',')]
        citations: Set[str] = {
            cit.key for cit in self.citations
            if cit.key in keys
            and self.bibliographies[cit.bibliography_key].list_ == 'citation'}
        if any(key in citations for key in keys):
            result_node = self.resolve_xref(
                env, fromdocname, builder, 'p', target, node, contnode)
            return [('p', result_node)]
        else:
            return []

    def get_all_cited_keys(self, docnames):
        """Yield all citation keys for given *docnames* in order, then
        ordered by citation order.
        """
        for citation_ref in sorted(
                self.citation_refs, key=lambda c: docnames.index(c.docname)):
            for key in citation_ref.keys:
                yield key

    def get_entries(
            self, bibfiles: List[str]) -> Iterable["Entry"]:
        """Return all bibliography entries from the bib files, unsorted (i.e.
        in order of appearance in the bib files.
        """
        for bibfile in bibfiles:
            for key in self.bibdata.bibfiles[bibfile].keys:
                yield self.bibdata.data.entries[key]

    def get_filtered_entries(
            self, bibliography_key: "BibliographyKey"
            ) -> Iterable[Tuple[str, "Entry"]]:
        """Return unsorted bibliography entries filtered by the filter
        expression.
        """
        bibliography = self.bibliographies[bibliography_key]
        for entry in self.get_entries(bibliography.bibfiles):
            key = bibliography.keyprefix + entry.key
            cited_docnames = {
                citation_ref.docname
                for citation_ref in self.citation_refs
                if key in citation_ref.keys
            }
            visitor = _FilterVisitor(
                entry=entry,
                docname=bibliography_key.docname,
                cited_docnames=cited_docnames)
            try:
                success = visitor.visit(bibliography.filter_)
            except ValueError as err:
                logger.warning(
                    "syntax error in :filter: expression; %s" % err,
                    location=(bibliography_key.docname, bibliography.line),
                    type="bibtex", subtype="filter_syntax_error")
                # recover by falling back to the default
                success = bool(cited_docnames)
            if success or entry.key in bibliography.keys:
                yield key, entry

    def get_sorted_entries(
            self, bibliography_key: "BibliographyKey", docnames: List[str]
            ) -> Iterable[Tuple[str, "Entry"]]:
        """Return filtered bibliography entries sorted by citation order."""
        entries = dict(
            self.get_filtered_entries(bibliography_key))
        for key in self.get_all_cited_keys(docnames):
            try:
                entry = entries.pop(key)
            except KeyError:
                pass
            else:
                yield key, entry
        # then all remaining keys, in order of bibliography file
        for key, entry in entries.items():
            yield key, entry

    def get_formatted_entries(
            self, bibliography_key: "BibliographyKey", docnames: List[str],
            tooltips: bool, tooltips_style: str
            ) -> Iterable[Tuple["Entry", "FormattedEntry",
                                Optional["FormattedEntry"]]]:
        """Get sorted bibliography entries along with their pybtex labels,
        with additional sorting and formatting applied from the pybtex style.
        """
        bibliography = self.bibliographies[bibliography_key]
        entries = dict(
            self.get_sorted_entries(bibliography_key, docnames))
        style: BaseStyle = cast("BaseStyle", pybtex.plugin.find_plugin(
            'pybtex.style.formatting', bibliography.style)())
        style2: Optional[BaseStyle] = (
            cast("BaseStyle", pybtex.plugin.find_plugin(
                'pybtex.style.formatting', tooltips_style)())
            if tooltips_style else style) if tooltips else None
        sorted_entries: Iterable[Entry] = style.sort(entries.values())
        labels = style.format_labels(sorted_entries)
        for label, entry in zip(labels, sorted_entries):
            try:
                yield (
                    entry,
                    style.format_entry(
                        bibliography.labelprefix + label, entry),
                    style2.format_entry(
                        bibliography.labelprefix + label, entry)
                    if style2 else None,
                )
            except FieldIsMissing as exc:
                logger.warning(
                    str(exc),
                    location=(bibliography_key.docname, bibliography.line),
                    type="bibtex", subtype="missing_field")
                formatted_error_entry = FormattedEntry(
                    entry.key, Tag('b', str(exc)),
                    bibliography.labelprefix + label)
                yield entry, formatted_error_entry, None
