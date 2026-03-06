from dataclasses import dataclass, field
from abc import ABC

import pybtex.plugin
from pybtex.richtext import Text, Tag
from sphinxcontrib.bibtex.style.template import (
    names, sentence, join, author_or_editor_or_title
)
from typing import (
    TYPE_CHECKING, Tuple, List, Union, Iterable, Optional, Dict
)

if TYPE_CHECKING:
    from pybtex.database import Entry
    from pybtex.richtext import BaseText
    from pybtex.style import FormattedEntry
    from pybtex.style.names import BaseNameStyle
    from pybtex.style.template import Node
    from sphinxcontrib.bibtex.richtext import ReferenceInfo


@dataclass
class BaseReferenceStyle(ABC):
    """Base class for citation reference styles.

    For consistency, all subclasses of this class must be decorated
    as a :class:`dataclass`,
    and must provide a type annotation and default value for all attributes
    (unless ``init=False`` is used, in which case they can be
    initialized in :meth:`~dataclass.__post_init__`).
    This allows client code to instantiate any reference style
    without needing to specify any arguments through the constructor.
    """

    # see https://stackoverflow.com/a/59987363 as to why this is here
    def __post_init__(self):
        pass

    def role_names(self) -> Iterable[str]:
        """Get list of role names supported by this style."""
        raise NotImplementedError

    def outer(
            self, role_name: str, children: List["BaseText"]) -> "Node":
        """Returns outer template for formatting the references."""
        raise NotImplementedError

    def inner(self, role_name: str) -> "Node":
        """Returns inner template for formatting the references."""
        raise NotImplementedError


def format_references(
        style: BaseReferenceStyle,
        role_name: str,
        references: Iterable[Tuple[
            "Entry", "FormattedEntry", "ReferenceInfo"]],
        ) -> "BaseText":
    """Format the list of references according to the given role.

    First formats each reference using the style's
    :meth:`~BaseReferenceStyle.get_inner` method,
    then joins all these formatted references together using
    the style's :meth:`~BaseReferenceStyle.get_outer` method.
    """
    children = [
        style.inner(role_name).format_data(
            data=dict(
                entry=entry,
                formatted_entry=formatted_entry,
                reference_info=info,
                style=style))
        for entry, formatted_entry, info in references]
    return style.outer(role_name, children).format()


@dataclass
class BracketStyle:
    """A class which provides brackets, as well as separators
    and a function to facilitate formatting of the outer template.
    """

    #: Left bracket.
    left: Union["BaseText", str] = '['

    #: Right bracket.
    right: Union["BaseText", str] = ']'

    #: Separators used for outer template (i.e. in between references
    #: if multiple keys are referenced in a single citation).
    sep: Union["BaseText", str] = ', '

    #: Separator for outer template, if only two items.
    sep2: Optional[Union["BaseText", str]] = None

    #: Separator for outer template, for last item if three or more items.
    last_sep: Optional[Union["BaseText", str]] = None

    def outer(
            self, children: List["BaseText"],
            brackets=False, capfirst=False) -> "Node":
        """Creates an outer template with separators,
        adding brackets if requested,
        and capitalizing the first word if requested.
        """
        return join[
            self.left if brackets else '',
            sentence(
                capfirst=capfirst,
                add_period=False,
                sep=self.sep,
                sep2=self.sep2,
                last_sep=self.last_sep,
            )[children],
            self.right if brackets else '',
        ]


@dataclass
class PersonStyle:
    """A class providing additional data and helper functions
    to facilitate formatting of person names.
    """

    #: Plugin name of the style used for formatting person names.
    style: str = 'last'

    #: Plugin class instance used for formatting person names.
    #: Automatically initialised from :attr:`style`.
    style_plugin: "BaseNameStyle" = field(init=False)

    #: Whether or not to abbreviate first names.
    abbreviate: bool = True

    #: Separator between persons.
    sep: Union["BaseText", str] = ', '

    #: Separator between persons, if only two persons.
    sep2: Optional[Union["BaseText", str]] = ' and '

    #: Separator between persons, for last person if three or more persons.
    last_sep: Optional[Union["BaseText", str]] = ', and '

    #: Abbreviation text if three or more persons.
    other: Optional[Union["BaseText", str]] = field(
        default_factory=lambda: Text(' ', Tag('em', 'et al.')))

    def __post_init__(self):
        self.style_plugin = pybtex.plugin.find_plugin(
            'pybtex.style.names', name=self.style)()

    def names(self, role: str, full: bool) -> "Node":
        """Returns a template formatting the persons with correct separators
        and using the full person list if so requested.
        """
        return names(
            role=role,
            sep=self.sep,
            sep2=self.sep2,
            last_sep=self.last_sep,
            other=None if full else self.other,
        )

    def author_or_editor_or_title(self, full: bool) -> "Node":
        """Returns a template formatting the author, falling back on editor
        or title if author is not specified.
        """
        return author_or_editor_or_title(
            sep=self.sep,
            sep2=self.sep2,
            last_sep=self.last_sep,
            other=None if full else self.other,
        )


@dataclass
class GroupReferenceStyle(BaseReferenceStyle):
    """Composes a group of reference styles into a single consistent style."""

    #: List of style types.
    styles: List[BaseReferenceStyle] = field(default_factory=list)

    #: Dictionary from role names to styles.
    #: Automatically initialized from :attr:`styles`.
    role_style: Dict[str, BaseReferenceStyle] = field(default_factory=dict)

    def __post_init__(self):
        super().__post_init__()
        self.role_style.update(
            (role_name, style)
            for style in self.styles
            for role_name in style.role_names()
        )

    def role_names(self):
        return self.role_style.keys()

    def outer(self, role_name: str, children: List["BaseText"]) -> "Node":
        """Gets the outer template associated with *role_name*
        in one of the :attr:`styles`.
        """
        style = self.role_style[role_name]
        return style.outer(role_name, children)

    def inner(self, role_name: str) -> "Node":
        """Gets the inner template associated with *role_name*
        in one of the :attr:`styles`.
        """
        style = self.role_style[role_name]
        return style.inner(role_name)
