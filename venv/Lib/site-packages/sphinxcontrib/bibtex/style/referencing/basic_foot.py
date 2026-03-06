from dataclasses import dataclass, field
from typing import TYPE_CHECKING, List, Iterable, Union
from sphinxcontrib.bibtex.style.template import footnote_reference, join
from . import BaseReferenceStyle, PersonStyle, BracketStyle

if TYPE_CHECKING:
    from pybtex.richtext import BaseText
    from pybtex.style.template import Node


@dataclass
class BasicFootParentheticalReferenceStyle(BaseReferenceStyle):
    """Parenthetical footnote reference."""

    def role_names(self) -> Iterable[str]:
        return ['p', 'ps']

    def outer(self, role_name: str, children: List["BaseText"]) -> "Node":
        return join[children]

    def inner(self, role_name: str) -> "Node":
        return footnote_reference


@dataclass
class BasicFootTextualReferenceStyle(BaseReferenceStyle):
    """Textual footnote reference."""

    #: Bracket style.
    bracket: BracketStyle = field(default_factory=BracketStyle)

    #: Person style.
    person: PersonStyle = field(default_factory=PersonStyle)

    #: Separator between text and reference.
    text_reference_sep: Union["BaseText", str] = ' '

    def role_names(self) -> Iterable[str]:
        return [f'{capfirst}t{full_author}'
                for capfirst in ['', 'c'] for full_author in ['', 's']]

    def outer(self, role_name: str, children: List["BaseText"]) -> "Node":
        return self.bracket.outer(
            children,
            brackets=False,
            capfirst='c' in role_name)

    def inner(self, role_name: str) -> "Node":
        return join(sep=self.text_reference_sep)[
            self.person.author_or_editor_or_title(full='s' in role_name),
            footnote_reference,
        ]
