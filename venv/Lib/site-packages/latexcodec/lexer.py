"""
    LaTeX Lexer
    ~~~~~~~~~~~

    This module contains all classes for lexing LaTeX code, as well as
    general purpose base classes for incremental LaTeX decoders and
    encoders, which could be useful in case you are writing your own
    custom LaTeX codec.

    .. autoclass:: Token(name, text)

    .. autoclass:: LatexLexer
       :show-inheritance:
       :members:

    .. autoclass:: LatexIncrementalLexer
       :show-inheritance:
       :members:

    .. autoclass:: LatexIncrementalDecoder
       :show-inheritance:
       :members:

    .. autoclass:: LatexIncrementalEncoder
       :show-inheritance:
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
import re
import unicodedata
from abc import ABC, ABCMeta
from typing import Any, Iterator, NamedTuple, Sequence, Tuple


class Token(NamedTuple):
    name: str
    text: str


# implementation note: we derive from IncrementalDecoder because this
# class serves excellently as a base class for incremental decoders,
# but of course we don't decode yet until later


class MetaRegexpLexer(ABCMeta):
    """Metaclass for :class:`RegexpLexer`. Compiles tokens into a
    regular expression.
    """

    def __init__(cls, name, bases, dct):
        super().__init__(name, bases, dct)
        regexp_string = "|".join(
            "(?P<" + name + ">" + regexp + ")"
            for name, regexp in getattr(cls, "tokens", [])
        )
        cls.regexp = re.compile(regexp_string, re.DOTALL)


class RegexpLexer(codecs.IncrementalDecoder, metaclass=MetaRegexpLexer):
    """Abstract base class for regexp based lexers."""

    emptytoken = Token("unknown", "")  #: The empty token.
    tokens: Sequence[Tuple[str, str]] = ()  #: Sequence of token regexps.
    errors: str  #: How to respond to errors.
    raw_buffer: Token  #: The raw buffer of this lexer.
    regexp: Any  #: Compiled regular expression.

    def __init__(self, errors: str = "strict") -> None:
        """Initialize the codec."""
        super().__init__(errors=errors)
        self.errors = errors
        self.reset()

    def reset(self) -> None:
        """Reset state."""
        # buffer for storing last (possibly incomplete) token
        self.raw_buffer = self.emptytoken

    def getstate(self) -> Any:
        """Get state."""
        return self.raw_buffer.text, 0

    def setstate(self, state: Any) -> None:
        """Set state. The *state* must correspond to the return value
        of a previous :meth:`getstate` call.
        """
        self.raw_buffer = Token("unknown", state[0])

    def get_raw_tokens(self, chars: str, final: bool = False) -> Iterator[Token]:
        """Yield tokens without any further processing. Tokens are one of:

        - ``\\<word>``: a control word (i.e. a command)
        - ``\\<symbol>``: a control symbol (i.e. \\^ etc.)
        - ``#<n>``: a parameter
        - a series of byte characters
        """
        if self.raw_buffer.text:
            chars = self.raw_buffer.text + chars
        self.raw_buffer = self.emptytoken
        for match in self.regexp.finditer(chars):
            # yield the buffer token
            if self.raw_buffer.text:
                yield self.raw_buffer
            # fill buffer with next token
            assert match.lastgroup is not None
            self.raw_buffer = Token(match.lastgroup, match.group(0))
        if final:
            for token in self.flush_raw_tokens():
                yield token

    def flush_raw_tokens(self) -> Iterator[Token]:
        """Flush the raw token buffer."""
        if self.raw_buffer.text:
            yield self.raw_buffer
            self.raw_buffer = self.emptytoken


class LatexLexer(RegexpLexer, ABC):
    """A very simple lexer for tex/latex."""

    # implementation note: every token **must** be decodable by inputenc
    tokens = [
        # match newlines and percent first, to ensure comments match correctly
        ("control_symbol_x2", r"[\\][\\]|[\\]%"),
        # comment: for ease, and for speed, we handle it as a token
        ("comment", r"%[^\n]*"),
        # control tokens
        # in latex, some control tokens skip following whitespace
        # ('control-word' and 'control-symbol')
        # others do not ('control-symbol-x')
        # XXX TBT says no control symbols skip whitespace (except '\ ')
        # XXX but tests reveal otherwise?
        ("control_word", r"[\\][a-zA-Z]+"),
        ("control_symbol", r"[\\][~" r"'" r'"` =^!.]'),
        # TODO should only match ascii
        ("control_symbol_x", r"[\\][^a-zA-Z]"),
        # parameter tokens
        # also support a lone hash so we can lex things like '#a'
        ("parameter", r"\#[0-9]|\#"),
        # any remaining characters; for ease we also handle space and
        # newline as tokens
        # XXX TBT does not mention \t to be a space character as well
        # XXX but tests reveal otherwise?
        ("space", r" |\t"),
        ("newline", r"\n"),
        ("mathshift", r"[$][$]|[$]"),
        # note: some chars joined together to make it easier to detect
        # symbols that have a special function (i.e. --, ---, etc.)
        (
            "chars",
            r"---|--|-|[`][`]" r"|['][']" r"|[?][`]|[!][`]"
            # separate chars because brackets are optional
            # e.g. fran\\c cais = fran\\c{c}ais in latex
            # so only way to detect \\c acting on c only is this way
            r"|(?![ %#$\n\t\\]).",
        ),
        # trailing garbage which we cannot decode otherwise
        # (such as a lone '\' at the end of a buffer)
        # is never emitted, but used internally by the buffer
        ("unknown", r"."),
    ]
    """List of token names, and the regular expressions they match."""


class LatexIncrementalLexer(LatexLexer, ABC):
    """A very simple incremental lexer for tex/latex code. Roughly
    follows the state machine described in Tex By Topic, Chapter 2.

    The generated tokens satisfy:

    * no newline characters: paragraphs are separated by '\\par'
    * spaces following control tokens are compressed
    """

    partoken = Token("control_word", "\\par")
    spacetoken = Token("space", " ")
    replacetoken = Token("chars", "\ufffd")
    curlylefttoken = Token("chars", "{")
    curlyrighttoken = Token("chars", "}")
    state: str
    inline_math: bool

    def reset(self) -> None:
        super().reset()
        # three possible states:
        # newline (N), skipping spaces (S), and middle of line (M)
        self.state = "N"
        # inline math mode?
        self.inline_math = False

    def getstate(self) -> Any:
        # state 'M' is most common, so let that be zero
        return (
            self.raw_buffer,
            {"M": 0, "N": 1, "S": 2}[self.state] | (4 if self.inline_math else 0),
        )

    def setstate(self, state: Any):
        self.raw_buffer = state[0]
        self.state = {0: "M", 1: "N", 2: "S"}[state[1] & 3]
        self.inline_math = bool(state[1] & 4)

    def get_tokens(self, chars: str, final: bool = False) -> Iterator[Token]:
        """Yield tokens while maintaining a state. Also skip
        whitespace after control words and (some) control symbols.
        Replaces newlines by spaces and \\par commands depending on
        the context.
        """
        # current position relative to the start of chars in the sequence
        # of bytes that have been decoded
        pos = -len(self.raw_buffer.text)
        for token in self.get_raw_tokens(chars, final=final):
            pos = pos + len(token.text)
            assert pos >= 0  # first token includes at least self.raw_buffer
            if token.name == "newline":
                if self.state == "N":
                    # if state was 'N', generate new paragraph
                    yield self.partoken
                elif self.state == "S":
                    # switch to 'N' state, do not generate a space
                    self.state = "N"
                elif self.state == "M":
                    # switch to 'N' state, generate a space
                    self.state = "N"
                    yield self.spacetoken
                else:
                    raise AssertionError("unknown tex state {0!r}".format(self.state))
            elif token.name == "space":
                if self.state == "N":
                    # remain in 'N' state, no space token generated
                    pass
                elif self.state == "S":
                    # remain in 'S' state, no space token generated
                    pass
                elif self.state == "M":
                    # in M mode, generate the space,
                    # but switch to space skip mode
                    self.state = "S"
                    yield token
                else:
                    raise AssertionError("unknown state {0!r}".format(self.state))
            elif token.name == "mathshift":
                self.inline_math = not self.inline_math
                self.state = "M"
                yield token
            elif token.name == "parameter":
                self.state = "M"
                yield token
            elif token.name == "control_word":
                # go to space skip mode
                self.state = "S"
                yield token
            elif token.name == "control_symbol":
                # go to space skip mode
                self.state = "S"
                yield token
            elif token.name == "control_symbol_x" or token.name == "control_symbol_x2":
                # don't skip following space, so go to M mode
                self.state = "M"
                yield token
            elif token.name == "comment":
                # no token is generated
                # note: comment does not include the newline
                self.state = "S"
            elif token.name == "chars":
                self.state = "M"
                yield token
            elif token.name == "unknown":
                if self.errors == "strict":
                    # current position within chars
                    # this is the position right after the unknown token
                    raise UnicodeDecodeError(
                        "latex",  # codec
                        chars.encode("utf8"),  # problematic input
                        pos - len(token.text),  # start of problematic token
                        pos,  # end of it
                        "unknown token {0!r}".format(token.text),
                    )
                elif self.errors == "ignore":
                    # do nothing
                    pass
                elif self.errors == "replace":
                    yield self.replacetoken
                else:
                    raise NotImplementedError(
                        "error mode {0!r} not supported".format(self.errors)
                    )
            else:
                raise AssertionError("unknown token name {0!r}".format(token.name))


class LatexIncrementalDecoder(LatexIncrementalLexer):
    """Simple incremental decoder. Transforms lexed LaTeX tokens into
    unicode.

    To customize decoding, subclass and override
    :meth:`get_unicode_tokens`.
    """

    inputenc = "ascii"
    """Input encoding. **Must** extend ascii."""

    def __init__(self, errors: str = "strict") -> None:
        super(LatexIncrementalDecoder, self).__init__(errors)
        self.decoder = codecs.getincrementaldecoder(self.inputenc)(errors)

    def decode_token(self, token: Token) -> str:
        """Returns the decoded token text.

        .. note::

           Control words get an extra space added at the back to make
           sure separation from the next token, so that decoded token
           sequences can be joined together.

           For example, the tokens ``'\\hello'`` and ``'world'``
           will correctly result in ``'\\hello world'`` (remember
           that LaTeX eats space following control words). If no space
           were added, this would wrongfully result in
           ``'\\helloworld'``.

        """
        text = token.text
        return text if token.name != "control_word" else text + " "

    def get_unicode_tokens(self, chars: str, final: bool = False) -> Iterator[str]:
        """Decode every token. Override to
        process the tokens in some other way (for example, for token
        translation).
        """
        for token in self.get_tokens(chars, final=final):
            yield self.decode_token(token)

    def udecode(self, bytes_: str, final: bool = False) -> str:
        """Decode LaTeX *bytes_* into a unicode string.

        This implementation calls :meth:`get_unicode_tokens` and joins
        the resulting unicode strings together.
        """
        return "".join(self.get_unicode_tokens(bytes_, final=final))

    def decode(self, bytes_: bytes, final: bool = False) -> str:
        """Decode LaTeX *bytes_* into a unicode string. Implementation uses
        :meth:`udecode`.
        """
        try:
            chars = self.decoder.decode(bytes_, final=final)
        except UnicodeDecodeError as e:
            # API requires that the encode method raises a ValueError
            # in this case
            raise ValueError(e)
        return self.udecode(chars, final)


class LatexIncrementalEncoder(codecs.IncrementalEncoder):
    """Simple incremental encoder for LaTeX. Transforms unicode into
    :class:`bytes`.

    To customize decoding, subclass and override
    :meth:`get_latex_bytes`.
    """

    inputenc = "ascii"
    """Input encoding. **Must** extend ascii."""

    buffer: str

    def __init__(self, errors: str = "strict") -> None:
        """Initialize the codec."""
        super().__init__(errors=errors)
        self.errors = errors
        self.reset()

    def reset(self) -> None:
        """Reset state."""
        # buffer for storing last (possibly incomplete) token
        self.buffer = ""

    def getstate(self) -> Any:
        """Get state."""
        return self.buffer

    def setstate(self, state: Any) -> None:
        """Set state. The *state* must correspond to the return value
        of a previous :meth:`getstate` call.
        """
        self.buffer = state

    def get_unicode_tokens(self, unicode_: str, final: bool = False) -> Iterator[str]:
        """Split unicode into tokens so that every token starts with a
        non-combining character.
        """
        if not isinstance(unicode_, str):
            raise TypeError(
                "expected unicode for encode input, but got {0} instead".format(
                    unicode_.__class__.__name__
                )
            )
        for c in unicode_:
            if not unicodedata.combining(c):
                for token in self.flush_unicode_tokens():
                    yield token
            self.buffer += c
        if final:
            for token in self.flush_unicode_tokens():
                yield token

    def flush_unicode_tokens(self) -> Iterator[str]:
        """Flush the buffer."""
        if self.buffer:
            yield self.buffer
            self.buffer = ""

    def get_latex_chars(self, unicode_: str, final: bool = False) -> Iterator[str]:
        """Encode every character. Override to
        process the unicode in some other way (for example, for character
        translation).
        """
        for token in self.get_unicode_tokens(unicode_, final=final):
            yield token

    def uencode(self, unicode_: str, final: bool = False) -> str:
        """Encode the *unicode_* string into LaTeX :class:`bytes`.

        This implementation calls :meth:`get_latex_chars` and joins
        the resulting :class:`bytes` together.
        """
        return "".join(self.get_latex_chars(unicode_, final=final))

    def encode(self, unicode_: str, final: bool = False) -> bytes:
        """Encode the *unicode_* string into LaTeX :class:`bytes`.

        This implementation calls :meth:`get_latex_chars` and joins
        the resulting :class:`bytes` together.
        """
        chars = self.uencode(unicode_, final)
        try:
            return chars.encode(self.inputenc, self.errors)
        except UnicodeEncodeError as e:
            # API requires that the encode method raises a ValueError
            # in this case
            raise ValueError(e)
