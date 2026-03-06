from dataclasses import dataclass

from typing import TYPE_CHECKING, List, Iterable
from sphinxcontrib.bibtex.style.template import join
from . import BaseReferenceStyle

if TYPE_CHECKING:
    from pybtex.richtext import BaseText
    from pybtex.style.template import Node


@dataclass
class ExtraEmptyReferenceStyle(BaseReferenceStyle):
    """A style which generates nothing, similar to LaTeX's nocite."""

    def role_names(self) -> Iterable[str]:
        return ['empty']

    def outer(self, role_name: str, children: List["BaseText"]) -> "Node":
        return join

    def inner(self, role_name: str) -> "Node":
        return join
