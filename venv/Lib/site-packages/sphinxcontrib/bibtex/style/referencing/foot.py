from dataclasses import dataclass, field
from typing import Union, TYPE_CHECKING

from . import PersonStyle, GroupReferenceStyle, BracketStyle
from .basic_foot import (
    BasicFootParentheticalReferenceStyle,
    BasicFootTextualReferenceStyle,
)

if TYPE_CHECKING:
    from pybtex.richtext import BaseText


@dataclass
class FootReferenceStyle(GroupReferenceStyle):
    """Textual or parenthetical reference using footnotes."""

    #: Bracket style for textual citations (:cite:t: and variations).
    bracket_textual: BracketStyle = field(default_factory=BracketStyle)

    #: Person style (applies to all relevant citation commands).
    person: PersonStyle = field(default_factory=PersonStyle)

    #: Separator between text and reference for textual citations.
    text_reference_sep: Union["BaseText", str] = ''

    def __post_init__(self):
        self.styles.extend([
            BasicFootParentheticalReferenceStyle(),
            BasicFootTextualReferenceStyle(
                bracket=self.bracket_textual,
                person=self.person,
                text_reference_sep=self.text_reference_sep,
            ),
        ])
        super().__post_init__()
