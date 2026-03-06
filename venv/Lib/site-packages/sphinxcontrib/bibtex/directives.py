"""
    .. autoclass:: BibliographyKey
        :members:

    .. autoclass:: BibliographyValue
        :members:

    .. autoclass:: BibliographyDirective

        .. automethod:: run
"""

from typing import TYPE_CHECKING, cast, NamedTuple, List, Dict
from docutils.parsers.rst import Directive
import docutils.parsers.rst.directives as directives

import ast  # parse(), used for filter
import docutils.nodes
import sphinx.util

from .bibfile import normpath_filename, _make_ids
from .nodes import bibliography as bibliography_node

if TYPE_CHECKING:
    from sphinx.environment import BuildEnvironment
    from .domain import BibtexDomain


logger = sphinx.util.logging.getLogger(__name__)


class BibliographyKey(NamedTuple):
    """Unique key for each bibliography directive."""
    docname: str  #: Name of the document where the bibliography resides.
    id_: str      #: The id of the bibliography node in the document.


class BibliographyValue(NamedTuple):
    """Contains information about a bibliography directive."""
    line: int            #: Line number of the directive in the document.
    bibfiles: List[str]  #: List of bib files for this directive.
    style: str           #: The pybtex style.
    list_: str           #: The list type.
    enumtype: str        #: The sequence type (for enumerated lists).
    start: int           #: The start of the sequence (for enumerated lists).
    labelprefix: str     #: String prefix for pybtex generated labels.
    keyprefix: str       #: String prefix for citation keys.
    filter_: ast.AST     #: Parsed filter expression.
    citation_nodes: Dict[str, docutils.nodes.Element]  #: key -> citation node
    keys: List[str]      #: Keys listed as content of the directive.


class BibliographyDirective(Directive):

    """Class for processing the :rst:dir:`bibliography` directive.

    Produces a
    :class:`~sphinxcontrib.bibtex.nodes.bibliography` node,
    along with (empty) citation nodes that will be formatted later in the
    *env-updated* stage, and inserted into the document in a post-transform.
    We cannot insert the citation nodes here because we do not yet know
    which keys have been cited.

    .. seealso::

       Further processing of the resulting
       :class:`~sphinxcontrib.bibtex.nodes.bibliography` node is done
       by
       :class:`~sphinxcontrib.bibtex.transforms.BibliographyTransform`.
    """

    required_arguments = 0
    optional_arguments = 1
    final_argument_whitespace = True
    has_content = True
    option_spec = {
        'cited': directives.flag,
        'notcited': directives.flag,
        'all': directives.flag,
        'filter': directives.unchanged,
        'style': directives.unchanged,
        'list': directives.unchanged,
        'enumtype': directives.unchanged,
        'start': (
            lambda value:
            directives.positive_int(value) if value != 'continue' else -1),
        'labelprefix': directives.unchanged,
        'keyprefix': directives.unchanged,
    }

    def _get_filter(self):
        """Get parsed filter from options."""
        env = cast("BuildEnvironment", self.state.document.settings.env)
        if "filter" in self.options:
            if "all" in self.options:
                logger.warning(":filter: overrides :all:",
                               location=(env.docname, self.lineno),
                               type="bibtex", subtype="filter_overrides")
            if "notcited" in self.options:
                logger.warning(":filter: overrides :notcited:",
                               location=(env.docname, self.lineno),
                               type="bibtex",
                               subtype="filter_overrides")
            if "cited" in self.options:
                logger.warning(":filter: overrides :cited:",
                               location=(env.docname, self.lineno),
                               type="bibtex", subtype="filter_overrides")
            try:
                return ast.parse(self.options["filter"])
            except SyntaxError:
                logger.warning(
                    "syntax error in :filter: expression" +
                    " (" + self.options["filter"] + "); "
                    "the option will be ignored",
                    location=(env.docname, self.lineno),
                    type="bibtex", subtype="filter_syntax_error")
                return ast.parse("cited")
        elif "all" in self.options:
            return ast.parse("True")
        elif "notcited" in self.options:
            return ast.parse("not cited")
        else:
            # the default filter: include only cited entries
            return ast.parse("cited")

    def run(self):
        """Process .bib files, set file dependencies, and create a
        node that is to be transformed to the entries of the
        bibliography.
        """
        env = cast("BuildEnvironment", self.state.document.settings.env)
        domain = cast("BibtexDomain", env.get_domain('cite'))
        filter_ = self._get_filter()
        if self.arguments:
            bibfiles = []
            for bibfile in self.arguments[0].split():
                normbibfile = normpath_filename(env, bibfile)
                if normbibfile not in domain.bibdata.bibfiles:
                    logger.warning(
                        "{0} not found or not configured"
                        " in bibtex_bibfiles".format(bibfile),
                        location=(env.docname, self.lineno),
                        type="bibtex", subtype="bibfile_error")
                else:
                    bibfiles.append(normbibfile)
        else:
            bibfiles = list(domain.bibdata.bibfiles.keys())
        for bibfile in bibfiles:
            env.note_dependency(bibfile)
        # generate nodes and ids
        keyprefix: str = self.options.get("keyprefix", "")
        list_: str = self.options.get("list", "citation")
        if list_ not in {"bullet", "enumerated", "citation"}:
            logger.warning(
                "unknown bibliography list type '{0}'.".format(list_),
                location=(env.docname, self.lineno),
                type="bibtex", subtype="list_type_error")
            list_ = "citation"
        if list_ in {"bullet", "enumerated"}:
            citation_node_class = docutils.nodes.list_item
        else:
            citation_node_class = docutils.nodes.citation
        bibliography_count = env.temp_data["bibtex_bibliography_count"] = \
            env.temp_data.get("bibtex_bibliography_count", 0) + 1
        ids = set(self.state.document.ids.keys())
        node = bibliography_node(
            '', docname=env.docname, ids=_make_ids(
                docname=env.docname, lineno=self.lineno,
                ids=ids,
                raw_id=env.app.config.bibtex_bibliography_id.format(
                    bibliography_count=bibliography_count)))
        self.state.document.note_explicit_target(node, node)
        # we only know which citations to included at resolve stage
        # but we need to know their ids before resolve stage
        # so for now we generate a node, and thus, an id, for every entry
        citation_nodes: Dict[str, docutils.nodes.Element] = {
            keyprefix + entry.key:
                citation_node_class(ids=_make_ids(
                    docname=env.docname,
                    lineno=self.lineno,
                    ids=ids,
                    raw_id=env.app.config.bibtex_cite_id.format(
                        bibliography_count=bibliography_count,
                        key=keyprefix + entry.key)))
            for entry in domain.get_entries(bibfiles)}
        for citation_node in citation_nodes.values():
            self.state.document.note_explicit_target(
                citation_node, citation_node)
        # check and get keys
        keys = []
        for key in self.content:
            if keyprefix + key not in citation_nodes:
                logger.warning('could not find bibtex key "%s"' % key,
                               location=(env.docname, self.lineno),
                               type="bibtex", subtype="key_not_found")
            else:
                keys.append(key)
        # create bibliography object
        bibliography = BibliographyValue(
            line=self.lineno,
            list_=list_,
            enumtype=self.options.get("enumtype", "arabic"),
            start=self.options.get("start", 1),
            style=self.options.get(
                "style", env.app.config.bibtex_default_style),
            filter_=filter_,
            labelprefix=self.options.get("labelprefix", ""),
            keyprefix=keyprefix,
            bibfiles=bibfiles,
            citation_nodes=citation_nodes,
            keys=keys,
        )
        bib_key = BibliographyKey(docname=env.docname, id_=node['ids'][0])
        assert bib_key not in domain.bibliographies
        domain.bibliographies[bib_key] = bibliography
        return [node]
