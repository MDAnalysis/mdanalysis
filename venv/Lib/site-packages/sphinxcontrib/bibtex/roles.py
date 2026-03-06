"""
    .. autoclass:: CitationRef
        :members:

    .. autoclass:: CiteRole
        :show-inheritance:

        .. automethod:: result_nodes
"""

import docutils.nodes

from typing import TYPE_CHECKING, cast, NamedTuple, List
from pybtex.plugin import find_plugin
from sphinx.roles import XRefRole

if TYPE_CHECKING:
    from .domain import BibtexDomain


class CitationRef(NamedTuple):
    """Information about a citation reference."""
    citation_ref_id: str  #: Unique id of this citation reference.
    docname: str          #: Document name.
    line: int             #: Line number.
    keys: List[str]       #: Citation keys (including key prefix).


class CiteRole(XRefRole):

    """Class for processing the :rst:role:`cite` role."""
    backend = find_plugin('pybtex.backends', 'docutils')()
    innernodeclass = docutils.nodes.inline

    def result_nodes(self, document, env, node, is_ref):
        """Associate the pending_xref with the cite domain,
        and note the cited citation keys.
        """
        if not node.get('refdomain'):
            assert node['reftype'] == 'cite'
            node['refdomain'] = 'cite'
            node['reftype'] = 'p'
        document.note_explicit_target(node, node)  # for backrefs
        domain = cast("BibtexDomain", env.get_domain('cite'))
        domain.citation_refs.append(CitationRef(
            citation_ref_id=node['ids'][0],
            docname=env.docname,
            line=document.line,
            keys=[key.strip() for key in self.target.split(',')],
        ))
        return [node], []
