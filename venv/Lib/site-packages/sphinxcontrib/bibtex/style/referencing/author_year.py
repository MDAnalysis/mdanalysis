from dataclasses import dataclass, field
from typing import Union, TYPE_CHECKING

from sphinxcontrib.bibtex.style.referencing import (
    BracketStyle, PersonStyle, GroupReferenceStyle
)
from .basic_author_year import (
    BasicAuthorYearParentheticalReferenceStyle,
    BasicAuthorYearTextualReferenceStyle,
)
from .extra_author import ExtraAuthorReferenceStyle
from .extra_label import ExtraLabelReferenceStyle
from .extra_year import ExtraYearReferenceStyle
from .extra_empty import ExtraEmptyReferenceStyle

if TYPE_CHECKING:
    from pybtex.richtext import BaseText


@dataclass
class AuthorYearReferenceStyle(GroupReferenceStyle):
    """Textual or parenthetical reference by author-year,
    or just by author, label, or year.
    """

    #: Bracket style for textual citations (:cite:t: and variations).
    bracket_textual: BracketStyle = field(default_factory=BracketStyle)

    #: Bracket style for parenthetical citations
    #: (:cite:p: and variations).
    bracket_parenthetical: BracketStyle = field(default_factory=BracketStyle)

    #: Bracket style for author citations
    #: (:cite:author: and variations).
    bracket_author: BracketStyle = field(default_factory=BracketStyle)

    #: Bracket style for label citations
    #: (:cite:label: and variations).
    bracket_label: BracketStyle = field(default_factory=BracketStyle)

    #: Bracket style for year citations
    #: (:cite:year: and variations).
    bracket_year: BracketStyle = field(default_factory=BracketStyle)

    #: Person style (applies to all relevant citation commands).
    person: PersonStyle = field(default_factory=PersonStyle)

    #: Separator between author and year for parenthetical citations.
    author_year_sep: Union["BaseText", str] = ', '

    #: Separator between text and reference for textual citations.
    text_reference_sep: Union["BaseText", str] = ' '

    def __post_init__(self):
        self.styles.extend([
            BasicAuthorYearParentheticalReferenceStyle(
                bracket=self.bracket_parenthetical,
                person=self.person,
                author_year_sep=self.author_year_sep,
            ),
            BasicAuthorYearTextualReferenceStyle(
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
