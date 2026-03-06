"""
    LaTeX Codec
    ~~~~~~~~~~~

    The :mod:`latexcodec.codec` module
    contains all classes and functions for LaTeX code
    translation. For practical use,
    you should only ever need to import the :mod:`latexcodec` module,
    which will automatically register the codec
    so it can be used by :meth:`str.encode`, :meth:`str.decode`,
    and any of the functions defined in the :mod:`codecs` module
    such as :func:`codecs.open` and so on.
    The other functions and classes
    are exposed in case someone would want to extend them.

    .. autofunction:: register

    .. autofunction:: find_latex

    .. autoclass:: LatexIncrementalEncoder
        :show-inheritance:
        :members:

    .. autoclass:: LatexIncrementalDecoder
        :show-inheritance:
        :members:

    .. autoclass:: LatexCodec
        :show-inheritance:
        :members:

    .. autoclass:: LatexUnicodeTable
        :members:
"""

# Copyright (c) 2003, 2008 David Eppstein
# Copyright (c) 2011-2020 Matthias C. M. Troffaes
#
# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation
# files (the "Software"), to deal in the Software without
# restriction, including without limitation the rights to use,
# copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following
# conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.

import codecs
import dataclasses
import importlib.resources as pkg_resources
import unicodedata
from codecs import CodecInfo
from typing import Any, Dict, Iterator, List, Optional, Tuple, Type, Union

from latexcodec import lexer


def register():
    """Register the :func:`find_latex` codec search function.

    .. seealso:: :func:`codecs.register`
    """
    codecs.register(find_latex)


# returns the codec search function
# this is used if latex_codec.py were to be placed in stdlib


def getregentry() -> Optional[CodecInfo]:
    """Encodings module API."""
    return find_latex("latex")


@dataclasses.dataclass
class UnicodeLatexTranslation:
    unicode: str
    latex: str
    encode: bool  #: Suitable for unicode -> latex.
    decode: bool  #: Suitable for latex -> unicode.
    text_mode: bool  #: Latex works in text mode.
    math_mode: bool  #: Latex works in math mode.


def load_unicode_latex_table() -> Iterator[UnicodeLatexTranslation]:
    with (pkg_resources.files("latexcodec") / "table.txt").open(
        "r", encoding="utf-8", errors="strict"
    ) as datafile:
        for line in datafile:
            marker, unicode_names, latex = line.rstrip("\r\n").split("\u0009")
            unicode = "".join(
                unicodedata.lookup(name) for name in unicode_names.split(",")
            )
            yield UnicodeLatexTranslation(
                unicode=unicode,
                latex=latex,
                encode=marker[1] in {"-", ">"},
                decode=marker[1] in {"-", "<"},
                text_mode=marker[0] in {"A", "T"},
                math_mode=marker[0] in {"A", "M"},
            )


class LatexUnicodeTable:
    """Tabulates a translation between LaTeX and unicode."""

    def __init__(self, lexer_):
        self.lexer: lexer.LatexIncrementalLexer = lexer_
        self.unicode_map: Dict[Tuple[lexer.Token, ...], str] = {}
        self.max_length: int = 0
        self.latex_map: Dict[str, Tuple[str, Tuple[lexer.Token, ...]]] = {}
        self.register_all()

    def register_all(self):
        """Register all symbols and their LaTeX equivalents
        (called by constructor).
        """
        # register special symbols
        self.register(
            UnicodeLatexTranslation(
                unicode="\n\n",
                latex=" \\par",
                encode=False,
                decode=True,
                text_mode=True,
                math_mode=False,
            )
        )
        self.register(
            UnicodeLatexTranslation(
                unicode="\n\n",
                latex="\\par",
                encode=False,
                decode=True,
                text_mode=True,
                math_mode=False,
            )
        )
        for trans in load_unicode_latex_table():
            self.register(trans)

    def register(self, trans: UnicodeLatexTranslation):
        """Register a correspondence between *unicode_text* and *latex_text*.

        :param UnicodeLatexTranslation trans: Description of translation.
        """
        if trans.math_mode and not trans.text_mode:
            # also register text version
            self.register(
                UnicodeLatexTranslation(
                    unicode=trans.unicode,
                    latex="$" + trans.latex + "$",
                    text_mode=True,
                    math_mode=False,
                    decode=trans.decode,
                    encode=trans.encode,
                )
            )
            self.register(
                UnicodeLatexTranslation(
                    unicode=trans.unicode,
                    latex=r"\(" + trans.latex + r"\)",
                    text_mode=True,
                    math_mode=False,
                    decode=trans.decode,
                    encode=trans.encode,
                )
            )
            # for the time being, we do not perform in-math substitutions
            return
        # tokenize, and register unicode translation
        self.lexer.reset()
        self.lexer.state = "M"
        tokens = tuple(self.lexer.get_tokens(trans.latex, final=True))
        if trans.decode:
            if tokens not in self.unicode_map:
                self.max_length = max(self.max_length, len(tokens))
                self.unicode_map[tokens] = trans.unicode
            # also register token variant with brackets, if appropriate
            # for instance, "\'{e}" for "\'e", "\c{c}" for "\c c", etc.
            # note: we do not remove brackets (they sometimes matter,
            # e.g. bibtex uses them to prevent lower case transformation)
            if (
                len(tokens) == 2
                and tokens[0].name.startswith("control")
                and tokens[1].name == "chars"
            ):
                self.register(
                    UnicodeLatexTranslation(
                        unicode=f"{{{trans.unicode}}}",
                        latex=f"{tokens[0].text}{{{tokens[1].text}}}",
                        decode=True,
                        encode=False,
                        math_mode=trans.math_mode,
                        text_mode=trans.text_mode,
                    )
                )
            if (
                len(tokens) == 4
                and tokens[0].text in {"$", r"\("}
                and tokens[1].name.startswith("control")
                and tokens[2].name == "chars"
                and tokens[3].text in {"$", r"\)"}
            ):
                # drop brackets in this case, since it is math mode
                self.register(
                    UnicodeLatexTranslation(
                        unicode=f"{trans.unicode}",
                        latex=f"{tokens[0].text}{tokens[1].text}"
                        f"{{{tokens[2].text}}}{tokens[3].text}",
                        decode=True,
                        encode=False,
                        math_mode=trans.math_mode,
                        text_mode=trans.text_mode,
                    )
                )
        if trans.encode and trans.unicode not in self.latex_map:
            assert len(trans.unicode) == 1
            self.latex_map[trans.unicode] = (trans.latex, tokens)


_LATEX_UNICODE_TABLE = LatexUnicodeTable(lexer.LatexIncrementalDecoder())

# incremental encoder does not need a buffer
# but decoder does


class LatexIncrementalEncoder(lexer.LatexIncrementalEncoder):
    """Translating incremental encoder for latex. Maintains a state to
    determine whether control spaces etc. need to be inserted.
    """

    emptytoken = lexer.Token("unknown", "")  #: The empty token.
    table = _LATEX_UNICODE_TABLE  #: Translation table.
    state: str

    def __init__(self, errors="strict"):
        super().__init__(errors=errors)
        self.reset()

    def reset(self):
        super(LatexIncrementalEncoder, self).reset()
        self.state = "M"

    def get_space_bytes(self, bytes_: str) -> Tuple[str, str]:
        """Inserts space bytes in space eating mode."""
        if self.state == "S":
            # in space eating mode
            # control space needed?
            if bytes_.startswith(" "):
                # replace by control space
                return "\\ ", bytes_[1:]
            else:
                # insert space (it is eaten, but needed for separation)
                return " ", bytes_
        else:
            return "", bytes_

    def _get_latex_chars_tokens_from_char(
        self, c: str
    ) -> Tuple[str, Tuple[lexer.Token, ...]]:
        # if ascii, try latex equivalents
        # (this covers \, #, &, and other special LaTeX characters)
        if ord(c) < 128:
            try:
                return self.table.latex_map[c]
            except KeyError:
                pass
        # next, try input encoding
        try:
            c.encode(self.inputenc, "strict")
        except UnicodeEncodeError:
            pass
        else:
            return c, (lexer.Token(name="chars", text=c),)
        # next, try latex equivalents of common unicode characters
        try:
            return self.table.latex_map[c]
        except KeyError:
            # translation failed
            if self.errors == "strict":
                raise UnicodeEncodeError(
                    "latex",  # codec
                    c,  # problematic input
                    0,
                    1,  # location of problematic character
                    "don't know how to translate {0} into latex".format(repr(c)),
                )
            elif self.errors == "ignore":
                return "", (self.emptytoken,)
            elif self.errors == "replace":
                # use the \\char command
                # this assumes
                # \usepackage[T1]{fontenc}
                # \usepackage[utf8]{inputenc}
                bytes_ = "{\\char" + str(ord(c)) + "}"
                return bytes_, (lexer.Token(name="chars", text=bytes_),)
            elif self.errors == "keep":
                return c, (lexer.Token(name="chars", text=c),)
            else:
                raise ValueError(
                    "latex codec does not support {0} errors".format(self.errors)
                )

    def get_latex_chars(self, unicode_: str, final: bool = False) -> Iterator[str]:
        if not isinstance(unicode_, str):
            raise TypeError(
                "expected unicode for encode input, but got {0} instead".format(
                    unicode_.__class__.__name__
                )
            )
        # convert character by character
        for pos, c in enumerate(unicode_):
            bytes_, tokens = self._get_latex_chars_tokens_from_char(c)
            space, bytes_ = self.get_space_bytes(bytes_)
            # update state
            if tokens and tokens[-1].name == "control_word":
                # we're eating spaces
                self.state = "S"
            elif tokens:
                self.state = "M"
            if space:
                yield space
            yield bytes_


class LatexIncrementalDecoder(lexer.LatexIncrementalDecoder):
    """Translating incremental decoder for LaTeX."""

    table = _LATEX_UNICODE_TABLE  #: Translation table.
    token_buffer: List[lexer.Token]  #: The token buffer of this decoder.

    def __init__(self, errors="strict"):
        lexer.LatexIncrementalDecoder.__init__(self, errors=errors)

    def reset(self):
        lexer.LatexIncrementalDecoder.reset(self)
        self.token_buffer = []

    # python codecs API does not support multibuffer incremental decoders

    def getstate(self) -> Any:
        raise NotImplementedError

    def setstate(self, state: Any) -> None:
        raise NotImplementedError

    def get_unicode_tokens(self, chars: str, final: bool = False) -> Iterator[str]:
        for token in self.get_tokens(chars, final=final):
            # at this point, token_buffer does not match anything
            self.token_buffer.append(token)
            # new token appended at the end, see if we have a match now
            # note: match is only possible at the *end* of the buffer
            # because all other positions have already been checked in
            # earlier iterations
            for i in range(len(self.token_buffer), 0, -1):
                last_tokens = tuple(self.token_buffer[-i:])  # last i tokens
                try:
                    unicode_text = self.table.unicode_map[last_tokens]
                except KeyError:
                    # no match: continue
                    continue
                else:
                    # match!! flush buffer, and translate last bit
                    # exclude last i tokens
                    for token2 in self.token_buffer[:-i]:
                        yield self.decode_token(token2)
                    yield unicode_text
                    self.token_buffer = []
                    break
            # flush tokens that can no longer match
            while len(self.token_buffer) >= self.table.max_length:
                yield self.decode_token(self.token_buffer.pop(0))
        # also flush the buffer at the end
        if final:
            for token in self.token_buffer:
                yield self.decode_token(token)
            self.token_buffer = []


class LatexCodec(codecs.Codec):
    IncrementalEncoder: Type[LatexIncrementalEncoder]
    IncrementalDecoder: Type[LatexIncrementalDecoder]

    def encode(
        self, unicode_: str, errors="strict"  # type: ignore
    ) -> Tuple[Union[bytes, str], int]:
        """Convert unicode string to LaTeX bytes."""
        encoder = self.IncrementalEncoder(errors=errors)
        return encoder.encode(unicode_, final=True), len(unicode_)

    def decode(self, bytes_: Union[bytes, str], errors="strict") -> Tuple[str, int]:
        """Convert LaTeX bytes to unicode string."""
        decoder = self.IncrementalDecoder(errors=errors)
        return decoder.decode(bytes_, final=True), len(bytes_)  # type: ignore


class UnicodeLatexIncrementalDecoder(LatexIncrementalDecoder):

    def decode(self, bytes_: str, final: bool = False) -> str:  # type: ignore
        return self.udecode(bytes_, final)


class UnicodeLatexIncrementalEncoder(LatexIncrementalEncoder):

    def encode(self, unicode_: str, final: bool = False) -> str:  # type: ignore
        return self.uencode(unicode_, final)


def find_latex(encoding: str) -> Optional[CodecInfo]:
    """Return a :class:`codecs.CodecInfo` instance for the requested
    LaTeX *encoding*, which must be equal to ``latex``,
    or to ``latex+<encoding>``
    where ``<encoding>`` describes another encoding.
    """
    if "_" in encoding:
        # Python 3.9 now normalizes "latex+latin1" to "latex_latin1"
        # https://bugs.python.org/issue37751
        encoding, _, inputenc_ = encoding.partition("_")
    else:
        encoding, _, inputenc_ = encoding.partition("+")
    if not inputenc_:
        inputenc_ = "ascii"
    if encoding == "latex":
        incremental_encoder = type(
            "incremental_encoder", (LatexIncrementalEncoder,), dict(inputenc=inputenc_)
        )
        incremental_decoder = type(
            "incremental_encoder", (LatexIncrementalDecoder,), dict(inputenc=inputenc_)
        )
    elif encoding == "ulatex":
        incremental_encoder = type(
            "incremental_encoder",
            (UnicodeLatexIncrementalEncoder,),
            dict(inputenc=inputenc_),
        )
        incremental_decoder = type(
            "incremental_encoder",
            (UnicodeLatexIncrementalDecoder,),
            dict(inputenc=inputenc_),
        )
    else:
        return None

    class Codec(LatexCodec):
        IncrementalEncoder = incremental_encoder
        IncrementalDecoder = incremental_decoder

    class StreamWriter(Codec, codecs.StreamWriter):
        pass

    class StreamReader(Codec, codecs.StreamReader):
        pass

    return codecs.CodecInfo(
        encode=Codec().encode,  # type: ignore
        decode=Codec().decode,  # type: ignore
        incrementalencoder=Codec.IncrementalEncoder,
        incrementaldecoder=Codec.IncrementalDecoder,
        streamreader=StreamReader,
        streamwriter=StreamWriter,
    )
