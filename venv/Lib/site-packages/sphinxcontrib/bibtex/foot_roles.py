"""
    .. autoclass:: FootCiteRole
        :show-inheritance:

        .. automethod:: result_nodes
"""

from typing import TYPE_CHECKING, cast, Tuple, List

import docutils.nodes
from docutils.nodes import make_id
from pybtex.plugin import find_plugin
from sphinx.roles import XRefRole
from sphinx.util.logging import getLogger

from .style.referencing import format_references
from .style.template import FootReferenceInfo
from .transforms import node_text_transform

if TYPE_CHECKING:
    from sphinx.environment import BuildEnvironment
    from .domain import BibtexDomain
    from .foot_domain import BibtexFootDomain

logger = getLogger(__name__)


class FootCiteRole(XRefRole):
    """Class for processing the :rst:role:`footcite` role."""

    def result_nodes(self, document: "docutils.nodes.document",
                     env: "BuildEnvironment", node: "docutils.nodes.Element",
                     is_ref: bool
                     ) -> Tuple[List["docutils.nodes.Node"],
                                List["docutils.nodes.system_message"]]:
        """Transform node into footnote references, and
        add footnotes to a node stored in the environment's temporary data
        if they are not yet present.

        .. seealso::

           The node containing all footnotes is inserted into the document by
           :meth:`.foot_directives.FootBibliographyDirective.run`.
        """
        if not node.get('refdomain'):
            assert node['reftype'] == 'footcite'
            node['refdomain'] = 'footcite'
            node['reftype'] = 'p'
        foot_domain = cast("BibtexFootDomain", self.env.get_domain('footcite'))
        keys = [key.strip() for key in self.target.split(',')]  # type: ignore
        try:
            foot_bibliography = env.temp_data["bibtex_foot_bibliography"]
        except KeyError:
            env.temp_data["bibtex_foot_bibliography"] = foot_bibliography = \
                foot_domain.bibliography_header.deepcopy()
        foot_old_refs = env.temp_data.setdefault("bibtex_foot_old_refs", set())
        foot_new_refs = env.temp_data.setdefault("bibtex_foot_new_refs", set())
        style = find_plugin(
            'pybtex.style.formatting',
            self.config.bibtex_default_style)()
        references = []
        domain = cast("BibtexDomain", self.env.get_domain('cite'))
        # count only incremented at directive, see foot_directives run method
        footbibliography_count = env.temp_data.setdefault(
            "bibtex_footbibliography_count", 0)
        footcite_names = env.temp_data.setdefault(
            "bibtex_footcite_names", {})
        for key in keys:
            entry = domain.bibdata.data.entries.get(key)
            if entry is not None:
                formatted_entry = style.format_entry(label='', entry=entry)
                if key not in (foot_old_refs | foot_new_refs):
                    footnote = docutils.nodes.footnote(auto=1)
                    # no automatic ids for footnotes: force non-empty template
                    template: str = \
                        env.app.config.bibtex_footcite_id \
                        if env.app.config.bibtex_footcite_id \
                        else "footcite-{key}"
                    raw_id = template.format(
                            footbibliography_count=footbibliography_count + 1,
                            key=entry.key)
                    # format name with make_id for consistency with cite role
                    name = make_id(raw_id)
                    footnote['names'] += [name]
                    footcite_names[entry.key] = name
                    footnote += domain.backend.paragraph(formatted_entry)
                    document.note_autofootnote(footnote)
                    document.note_explicit_target(footnote, footnote)
                    node_text_transform(footnote)
                    foot_bibliography += footnote
                    foot_new_refs.add(key)
                references.append(
                    (entry, formatted_entry,
                     FootReferenceInfo(
                         key=entry.key, refname=footcite_names[entry.key],
                         document=document)))
            else:
                logger.warning('could not find bibtex key "%s"' % key,
                               location=(env.docname, self.lineno),
                               type="bibtex", subtype="key_not_found")
        ref_nodes = format_references(
            foot_domain.reference_style, node['reftype'], references
        ).render(domain.backend)
        return ref_nodes, []
