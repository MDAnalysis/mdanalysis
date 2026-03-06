from abc import ABC
from typing import TypeVar, Generic
from pybtex.richtext import BaseMultipartText, BaseText


ReferenceInfo = TypeVar('ReferenceInfo')
"""Generic type parameter for types that store reference information.
To be implemented by clients. See for instance
:class:`~sphinxcontrib.bibtex.domain.SphinxReferenceInfo`.
"""


class BaseReferenceText(BaseMultipartText, Generic[ReferenceInfo], ABC):
    """Generic rich text element for citation references.
    Instances store some extra reference info that can be used when formatting.
    This base class renders its children without further formatting.
    Implementations must create a derivation from this class which
    overrides the *render* method to create the desired output.
    See for instance
    :class:`~sphinxcontrib.bibtex.domain.SphinxReferenceText`.
    """

    def __init__(self, info: ReferenceInfo, *parts: "BaseText"):
        self.info = (info,)
        super().__init__(*parts)
