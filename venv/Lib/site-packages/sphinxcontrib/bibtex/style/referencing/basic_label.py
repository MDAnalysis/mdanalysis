from dataclasses import dataclass, field

from typing import TYPE_CHECKING, List, Iterable, Union
from sphinxcontrib.bibtex.style.template import reference, entry_label, join
from . import BracketStyle, PersonStyle, BaseReferenceStyle

if TYPE_CHECKING:
    from pybtex.richtext import BaseText
    from pybtex.style.template import Node


@dataclass
class BasicLabelParentheticalReferenceStyle(BaseReferenceStyle):
    """Reference by label if parenthetical,
    and by author and label if textual.
    """

    #: Bracket style.
    bracket: BracketStyle = field(default_factory=BracketStyle)

    def role_names(self) -> Iterable[str]:
        return [f'p{full_author}' for full_author in ['', 's']]

    def outer(self, role_name: str, children: List["BaseText"]) -> "Node":
        return self.bracket.outer(
            children,
            brackets=True,
            capfirst=False)

    def inner(self, role_name: str) -> "Node":
        return reference[entry_label]


@dataclass
class BasicLabelTextualReferenceStyle(BaseReferenceStyle):
    """Reference by label if parenthetical,
    and by author and label if textual.
    """

    #: Bracket style.
    bracket: BracketStyle = field(default_factory=BracketStyle)

    #: Person style.
    person: PersonStyle = field(default_factory=PersonStyle)

    #: Separator between text and reference.
    text_reference_sep: Union["BaseText", str] = ' '

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
            join[
                self.bracket.left,
                reference[entry_label],
                self.bracket.right
            ]
        ]
