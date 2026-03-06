"""
    .. autoclass:: bibliography
    .. autoclass:: raw_latex
    .. autofunction:: visit_raw_latex
    .. autofunction:: depart_raw_latex
"""

from docutils import nodes
from sphinx.writers.latex import LaTeXTranslator


class bibliography(nodes.General, nodes.Element):
    """Node for representing a bibliography. Replaced by a list of
    citations by
    :class:`~sphinxcontrib.bibtex.transforms.BibliographyTransform`.
    """
    pass


class raw_latex(
        nodes.Special, nodes.Inline, nodes.PreBibliographic,
        nodes.FixedTextElement):
    """Node for representing raw latex data."""
    pass


def visit_raw_latex(self: LaTeXTranslator, node: raw_latex):
    """Called when entering a raw_latex node. Appends the node's raw source
    to the latex body.
    """
    self.body.append(node.rawsource)


def depart_raw_latex(self: LaTeXTranslator, node: raw_latex):
    """Called when leaving a raw_latex node."""
    pass
