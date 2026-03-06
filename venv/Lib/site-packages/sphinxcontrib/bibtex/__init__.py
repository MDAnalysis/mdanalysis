"""
    .. autofunction:: setup
"""

from typing import Any, Dict

from sphinx.application import Sphinx

from .domain import BibtexDomain
from .foot_domain import BibtexFootDomain
from .nodes import bibliography, raw_latex, visit_raw_latex, depart_raw_latex
from .roles import CiteRole
from .directives import BibliographyDirective
from .transforms import BibliographyTransform
from .foot_roles import FootCiteRole
from .foot_directives import FootBibliographyDirective


def setup(app: Sphinx) -> Dict[str, Any]:
    """Set up the bibtex extension:

    * register config values
    * register directives
    * register nodes
    * register roles
    * register transforms
    * connect events to functions
    """
    app.add_config_value("bibtex_default_style", "alpha", "html")
    app.add_config_value("bibtex_tooltips", True, "html")
    app.add_config_value("bibtex_tooltips_style", "", "html")
    app.add_config_value("bibtex_bibfiles", None, "html")
    app.add_config_value("bibtex_encoding", "utf-8-sig", "html")
    app.add_config_value("bibtex_bibliography_header", "", "html")
    app.add_config_value("bibtex_footbibliography_header", "", "html")
    app.add_config_value("bibtex_reference_style", "label", "env")
    app.add_config_value("bibtex_foot_reference_style", "foot", "env")
    app.add_config_value("bibtex_cite_id", "", "html")
    app.add_config_value("bibtex_footcite_id", "", "html")
    app.add_config_value("bibtex_bibliography_id", "", "html")
    app.add_config_value("bibtex_footbibliography_id", "", "html")
    app.add_domain(BibtexDomain)
    app.add_directive("bibliography", BibliographyDirective)
    app.add_role("cite", CiteRole())
    app.add_node(bibliography, override=True)
    app.add_node(raw_latex, latex=(visit_raw_latex, depart_raw_latex),
                 override=True)
    app.add_post_transform(BibliographyTransform)
    app.add_domain(BibtexFootDomain)
    app.add_directive("footbibliography", FootBibliographyDirective)
    app.add_role("footcite", FootCiteRole())

    return {
        'version': '2.5.0',
        'env_version': 9,
        'parallel_read_safe': True,
        'parallel_write_safe': True,
        }
