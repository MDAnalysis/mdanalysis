"""
    Domain for footnote citations.

    .. autoclass:: BibtexFootDomain
        :members:
"""

from typing import TYPE_CHECKING, List, Dict, Tuple

import docutils.nodes
import docutils.utils

import sphinxcontrib.bibtex.plugin
import sphinx.util

from sphinx.domains import Domain, ObjType
from sphinx.locale import _

from .foot_roles import FootCiteRole
from .domain import parse_header
from .style.referencing import BaseReferenceStyle

if TYPE_CHECKING:
    from sphinx.environment import BuildEnvironment
    from sphinx.addnodes import pending_xref
    from sphinx.builders import Builder

logger = sphinx.util.logging.getLogger(__name__)


class BibtexFootDomain(Domain):
    """Sphinx domain for footnote citations."""

    name = 'footcite'
    label = 'BibTeX Footnote Citations'
    data_version = 0
    initial_data = dict(
        bibliography_header=docutils.nodes.container(),
    )
    reference_style: BaseReferenceStyle

    @property
    def bibliography_header(self) -> docutils.nodes.Element:
        return self.data['bibliography_header']

    def __init__(self, env: "BuildEnvironment"):
        # set up referencing style
        style = sphinxcontrib.bibtex.plugin.find_plugin(
                'sphinxcontrib.bibtex.style.referencing',
                env.app.config.bibtex_foot_reference_style)
        self.reference_style = style()
        # set up object types and roles for referencing style
        role_names = self.reference_style.role_names()
        self.object_types = dict(
            citation=ObjType(_('citation'), *role_names, searchprio=-1),
        )
        self.roles = dict((name, FootCiteRole()) for name in role_names)
        # initialize the domain
        super().__init__(env)
        # parse bibliography header
        header = getattr(env.app.config, "bibtex_footbibliography_header")
        if header:
            self.data["bibliography_header"] += \
                parse_header(header, "foot_bibliography_header")

    def merge_domaindata(self, docnames: List[str], otherdata: Dict) -> None:
        """Merge in data regarding *docnames* from domain data
        inventory *otherdata*.

        As there is no document specific data for this domain, this function
        does nothing.
        """
        pass

    def resolve_any_xref(self, env: "BuildEnvironment", fromdocname: str,
                         builder: "Builder", target: str,
                         node: "pending_xref", contnode: docutils.nodes.Element
                         ) -> List[Tuple[str, docutils.nodes.Element]]:
        """Resolve the pending reference *node* with the given *target*,
        where the reference comes from an "any" role.

        Since citation references are resolved to regular citations,
        and not to footnote citations,
        this implementation simply returns an empty list.
        """
        return []
