from dataclasses import dataclass, field
from typing import Union, TYPE_CHECKING

from sphinxcontrib.bibtex.style.referencing import (
    BracketStyle, PersonStyle, GroupReferenceStyle
)
from .basic_super import (
    BasicSuperParentheticalReferenceStyle,
    BasicSuperTextualReferenceStyle,
)
from .extra_author import ExtraAuthorReferenceStyle
from .extra_label import ExtraLabelReferenceStyle
from .extra_year import ExtraYearReferenceStyle
from .extra_empty import ExtraEmptyReferenceStyle

if TYPE_CHECKING:
    from pybtex.richtext import BaseText


@dataclass
class SuperReferenceStyle(GroupReferenceStyle):
    """Textual or parenthetical reference by superscripted label,
    or just by author, label, or year.
    """

    #: Bracket style for textual citations (:cite:t: and variations).
    bracket_textual: BracketStyle = field(
        default_factory=lambda: BracketStyle(left='', right='', sep=', '))

    #: Bracket style for parenthetical citations
    #: (:cite:p: and variations).
    bracket_parenthetical: BracketStyle = field(
        default_factory=lambda: BracketStyle(left='', right='', sep=','))

    #: Bracket style for author citations
    #: (:cite:author: and variations).
    bracket_author: BracketStyle = field(default_factory=BracketStyle)

    #: Bracket style for label citations
    #: (:cite:label: and variations).
    bracket_label: BracketStyle = field(default_factory=BracketStyle)

    #: Bracket style for year citations
    #: (:cite:year: and variations).
    bracket_year: BracketStyle = field(default_factory=BracketStyle)

    #: Person style.
    person: PersonStyle = field(default_factory=PersonStyle)

    #: Separator between text and reference for textual citations.
    text_reference_sep: Union["BaseText", str] = ''

    def __post_init__(self):
        self.styles.extend([
            BasicSuperParentheticalReferenceStyle(
                bracket=self.bracket_parenthetical, person=self.person),
            BasicSuperTextualReferenceStyle(
                bracket=self.bracket_textual,
                person=self.person,
                text_reference_sep=self.text_reference_sep,
            ),
            ExtraAuthorReferenceStyle(
                bracket=self.bracket_author, person=self.person),
            ExtraLabelReferenceStyle(bracket=self.bracket_label),
            ExtraYearReferenceStyle(bracket=self.bracket_year),
            ExtraEmptyReferenceStyle(),
        ])
        super().__post_init__()
