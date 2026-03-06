# Copyright (c) 2006-2021  Andrey Golovizin
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

r"""(simple but) rich text formatting tools

Usage:

>>> t = Text('this ', 'is a ', Tag('em', 'very'), Text(' rich', ' text'))
>>> print(t.render_as('latex'))
this is a \emph{very} rich text
>>> print(str(t))
this is a very rich text
>>> t = t.capitalize().add_period()
>>> print(t.render_as('latex'))
This is a \emph{very} rich text.
>>> print(str(t))
This is a very rich text.
>>> print(Symbol('ndash').render_as('latex'))
--
>>> t = Text('Some ', Tag('em', Text('nested ', Tag('tt', 'Text', Text(' objects')))), '.')
>>> print(t.render_as('latex'))
Some \emph{nested \texttt{Text objects}}.
>>> print(str(t))
Some nested Text objects.
>>> t = t.upper()
>>> print(t.render_as('latex'))
SOME \emph{NESTED \texttt{TEXT OBJECTS}}.
>>> print(str(t))
SOME NESTED TEXT OBJECTS.

>>> t = Text(', ').join(['one', 'two', Tag('em', 'three')])
>>> print(t.render_as('latex'))
one, two, \emph{three}
>>> print(str(t))
one, two, three
>>> t = Text(Symbol('nbsp')).join(['one', 'two', Tag('em', 'three')])
>>> print(t.render_as('latex'))
one~two~\emph{three}
>>> print(str(t))
one<nbsp>two<nbsp>three
"""
from __future__ import absolute_import, unicode_literals

import itertools
import warnings
from abc import ABCMeta, abstractmethod

from pybtex import textutils
from pybtex.utils import collect_iterable, deprecated


# workaround for doctests in Python 2/3
def str_repr(string):
    """
    >>> print(str_repr('test'))
    'test'
    >>> print(str_repr(u'test'))
    'test'
    """

    result = repr(string)
    if result.startswith('u'):
        return result[1:]
    else:
        return result


def ensure_text(value):
    if isinstance(value, str):
        return String(value)
    elif isinstance(value, BaseText):
        return value
    else:
        bad_type = type(value).__name__
        raise ValueError('parts must be strings or BaseText instances, not ' + bad_type)


class BaseText(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def __str__(self):
        raise NotImplementedError

    @abstractmethod
    def __eq__(self, other):
        raise NotImplementedError

    def __ne__(self, other):
        return not self == other

    @abstractmethod
    def __len__(self):
        raise NotImplementedError

    @abstractmethod
    def __contains__(self, item):
        raise NotImplementedError

    @abstractmethod
    def __getitem__(self, key):
        raise NotImplementedError

    def __add__(self, other):
        """
        Concatenate this Text with another Text or string.

        >>> Text('Longcat is ') + Tag('em', 'long')
        Text('Longcat is ', Tag('em', 'long'))
        """

        return Text(self, other)

    def append(self, text):
        """
        Append text to the end of this text.

        Normally, this is the same as concatenating texts with +,
        but for tags and similar objects the appended text is placed _inside_ the tag.

        >>> text = Tag('em', 'Look here')
        >>> print((text +  '!').render_as('html'))
        <em>Look here</em>!
        >>> print(text.append('!').render_as('html'))
        <em>Look here!</em>
        """

        return self + text

    def join(self, parts):
        """Join a list using this text (like string.join)

        >>> letters = ['a', 'b', 'c']
        >>> print(str(String('-').join(letters)))
        a-b-c
        >>> print(str(String('-').join(iter(letters))))
        a-b-c
        """

        if not parts:
            return Text()
        joined = []
        for part in parts:
            if joined:
                joined.append(self)
            joined.append(part)
        return Text(*joined)

    @abstractmethod
    def split(self, sep=None, keep_empty_parts=None):
        raise NotImplementedError

    @abstractmethod
    def startswith(self, prefix):
        """
        Return True if string starts with the prefix,
        otherwise return False.

        prefix can also be a tuple of suffixes to look for.
        """

        raise NotImplementedError

    @abstractmethod
    def endswith(self, suffix):
        """
        Return True if the string ends with the specified suffix,
        otherwise return False.

        suffix can also be a tuple of suffixes to look for.
        """

        raise NotImplementedError

    @abstractmethod
    def isalpha(self):
        raise NotImplementedError

    def add_period(self, period='.'):
        """
        Add a period to the end of text, if the last character is not ".", "!" or "?".

        >>> text = Text("That's all, folks")
        >>> print(str(text.add_period()))
        That's all, folks.

        >>> text = Text("That's all, folks!")
        >>> print(str(text.add_period()))
        That's all, folks!

        """

        if self and not textutils.is_terminated(self):
            return self.append(period)
        else:
            return self

    def abbreviate(self):
        def abbreviate_word(word):
            if word.isalpha():
                return word[0].add_period()
            else:
                return word

        parts = self.split(textutils.delimiter_re)
        return String('').join(abbreviate_word(part) for part in parts)

    def capfirst(self):
        """
        Capitalize the first letter of the text.

        >>> Text(Tag('em', 'long Cat')).capfirst()
        Text(Tag('em', 'Long Cat'))

        """
        return self[:1].upper() + self[1:]

    def capitalize(self):
        """
        Capitalize the first letter of the text and lowercase the rest.

        >>> Text(Tag('em', 'LONG CAT')).capitalize()
        Text(Tag('em', 'Long cat'))

        """
        return self[:1].upper() + self[1:].lower()

    @abstractmethod
    def lower(self):
        raise NotImplementedError

    @abstractmethod
    def upper(self):
        raise NotImplementedError

    @abstractmethod
    def render(self, backend):
        raise NotImplementedError

    def render_as(self, backend_name):
        r"""
        Render this :py:class:`Text` into markup.
        This is a wrapper method that loads a formatting backend plugin
        and calls :py:meth:`Text.render`.

        >>> text = Text('Longcat is ', Tag('em', 'looooooong'), '!')
        >>> print(text.render_as('html'))
        Longcat is <em>looooooong</em>!
        >>> print(text.render_as('latex'))
        Longcat is \emph{looooooong}!
        >>> print(text.render_as('text'))
        Longcat is looooooong!

        :param backend_name: The name of the output backend (like ``"latex"`` or
            ``"html"``).

        """
        from pybtex.plugin import find_plugin
        backend_cls = find_plugin('pybtex.backends', backend_name)
        return self.render(backend_cls())

    def _unpack(self):
        """
        For Text object, iterate over all text parts.
        Else, yield the object itself.

        Used for unpacking Text objects passed as children to another Text object.
        """

        yield self

    def _typeinfo(self):
        """

        Return the type of this object and its parameters
        (not including the actual text content).

        Used for:

        - merging similar tags together (<em>A</em><em>B</em> -> <em>AB</em>),
        - creating similar text objects with different text content.

        """

        return None, ()


class BaseMultipartText(BaseText):
    info = ()

    def __init__(self, *parts):
        """Create a text object consisting of one or more parts.

        Empty parts are ignored:

        >>> Text() == Text('') == Text('', '', '')
        True
        >>> Text('Word', '') == Text('Word')
        True

        Text() objects are unpacked and their children are included directly:

        >>> Text(Text('Multi', ' '), Tag('em', 'part'), Text(' ', Text('text!')))
        Text('Multi ', Tag('em', 'part'), ' text!')
        >>> Tag('strong', Text('Multi', ' '), Tag('em', 'part'), Text(' ', 'text!'))
        Tag('strong', 'Multi ', Tag('em', 'part'), ' text!')

        Similar objects are merged together:

        >>> Text('Multi', Tag('em', 'part'), Text(Tag('em', ' ', 'text!')))
        Text('Multi', Tag('em', 'part text!'))
        >>> Text('Please ', HRef('/', 'click'), HRef('/', ' here'), '.')
        Text('Please ', HRef('/', 'click here'), '.')
        """

        parts = (ensure_text(part) for part in parts)
        nonempty_parts = (part for part in parts if part)
        unpacked_parts = itertools.chain(*[part._unpack() for part in nonempty_parts])
        merged_parts = self._merge_similar(unpacked_parts)
        self.parts = list(merged_parts)
        self.length = sum(len(part) for part in self.parts)

    def __str__(self):
        return ''.join(str(part) for part in self.parts)

    def __eq__(self, other):
        """
        Rich text objects support equality comparison:

        >>> Text('Cat') == Text('cat')
        False
        >>> Text('Cat') == Text('Cat')
        True

        """
        return (
            isinstance(other, BaseText) and
            self._typeinfo() == other._typeinfo() and
            self.parts == other.parts
        )

    def __len__(self):
        """
        ``len(text)`` returns the number of characters in the text, ignoring
        the markup:

        >>> len(Text('Long cat'))
        8
        >>> len(Text(Tag('em', 'Long'), ' cat'))
        8
        >>> len(Text(HRef('http://example.com/', 'Long'), ' cat'))
        8

        """
        return self.length

    def __contains__(self, item):
        """
        ``value in text`` returns ``True`` if any part of the ``text``
        contains the substring ``value``:

        >>> 'Long cat' in Text('Long cat!')
        True

        Substrings splitted across multiple text parts are not matched:

        >>> 'Long cat' in Text(Tag('em', 'Long'), 'cat!')
        False

        """
        if not isinstance(item, str):
            raise TypeError(item)
        return not item or any(part.__contains__(item) for part in self.parts)

    def __getitem__(self, key):
        """
        Slicing and extracting characters works like with regular strings,
        formatting is preserved.

        >>> Text('Longcat is ', Tag('em', 'looooooong!'))[:15]
        Text('Longcat is ', Tag('em', 'looo'))
        >>> Text('Longcat is ', Tag('em', 'looooooong!'))[-1]
        Text(Tag('em', '!'))
        """

        if isinstance(key, int):
            start = key
            end = None
        elif isinstance(key, slice):
            start, end, step = key.indices(len(self))
            if step != 1:
                raise NotImplementedError
        else:
            raise TypeError(key, type(key))

        if start < 0:
            start = len(self) + start
        if end is None:
            end = start + 1
        if end < 0:
            end = len(self) + end
        return self._slice_end(len(self) - start)._slice_beginning(end - start)

    def _slice_beginning(self, slice_length):
        """
        Return a text consistng of the first slice_length characters
        of this text (with formatting preserved).
        """

        parts = []
        length = 0
        for part in self.parts:
            if length + len(part) > slice_length:
                parts.append(part[:slice_length - length])
                break
            else:
                parts.append(part)
                length += len(part)
        return self._create_similar(parts)

    def _slice_end(self, slice_length):
        """
        Return a text consistng of the last slice_length characters
        of this text (with formatting preserved).
        """

        parts = []
        length = 0
        for part in reversed(self.parts):
            if length + len(part) > slice_length:
                parts.append(part[len(part) - (slice_length - length):])
                break
            else:
                parts.append(part)
                length += len(part)
        return self._create_similar(reversed(parts))

    def append(self, text):
        """
        Append text to the end of this text.

        For Tags, HRefs, etc. the appended text is placed *inside* the tag.

        >>> text = Tag('strong', 'Chuck Norris')
        >>> print((text +  ' wins!').render_as('html'))
        <strong>Chuck Norris</strong> wins!
        >>> print(text.append(' wins!').render_as('html'))
        <strong>Chuck Norris wins!</strong>
        """

        return self._create_similar(self.parts + [text])

    @collect_iterable
    def split(self, sep=None, keep_empty_parts=None):
        """
        >>> Text('a + b').split()
        [Text('a'), Text('+'), Text('b')]

        >>> Text('a, b').split(', ')
        [Text('a'), Text('b')]
        """

        if keep_empty_parts is None:
            keep_empty_parts = sep is not None

        tail = [''] if keep_empty_parts else []
        for part in self.parts:
            split_part = part.split(sep, keep_empty_parts=True)
            if not split_part:
                continue
            for item in split_part[:-1]:
                if tail:
                    yield self._create_similar(tail + [item])
                    tail = []
                else:
                    if item or keep_empty_parts:
                        yield self._create_similar([item])
            tail.append(split_part[-1])
        if tail:
            tail_text = self._create_similar(tail)
            if tail_text or keep_empty_parts:
                yield tail_text

    def startswith(self, prefix):
        """
        Return True if the text starts with the given prefix.

        >>> Text('Longcat!').startswith('Longcat')
        True

        Prefixes split across multiple parts are not matched:

        >>> Text(Tag('em', 'Long'), 'cat!').startswith('Longcat')
        False

        """

        if not self.parts:
            return False
        else:
            return self.parts[0].startswith(prefix)

    def endswith(self, suffix):
        """
        Return True if the text ends with the given suffix.

        >>> Text('Longcat!').endswith('cat!')
        True

        Suffixes split across multiple parts are not matched:

        >>> Text('Long', Tag('em', 'cat'), '!').endswith('cat!')
        False

        """

        if not self.parts:
            return False
        else:
            return self.parts[-1].endswith(suffix)

    def isalpha(self):
        """
        Return True if all characters in the string are alphabetic and there is
        at least one character, False otherwise.
        """
        return bool(self) and all(part.isalpha() for part in self.parts)

    def lower(self):
        """
        Convert rich text to lowercase.

        >>> Text(Tag('em', 'Long cat')).lower()
        Text(Tag('em', 'long cat'))
        """

        return self._create_similar(part.lower() for part in self.parts)

    def upper(self):
        """
        Convert rich text to uppsercase.

        >>> Text(Tag('em', 'Long cat')).upper()
        Text(Tag('em', 'LONG CAT'))
        """
        return self._create_similar(part.upper() for part in self.parts)

    def render(self, backend):
        """
        Render this :py:class:`Text` into markup.

        :param backend: The formatting backend (an instance of
            :py:class:`pybtex.backends.BaseBackend`).
        """

        rendered_list = [part.render(backend) for part in self.parts]
        assert all(isinstance(item, backend.RenderType)
                   for item in rendered_list)
        return backend.render_sequence(rendered_list)

    def _typeinfo(self):
        """Return the type and the parameters used to create this text object.

        >>> text = Tag('strong', 'Heavy rain!')
        >>> text._typeinfo() == (Tag, ('strong',))
        True

        """

        return type(self), self.info

    def _create_similar(self, parts):
        """
        Create a new text object of the same type with the same parameters,
        with different text content.

        >>> text = Tag('strong', 'Bananas!')
        >>> text._create_similar(['Apples!'])
        Tag('strong', 'Apples!')
        """

        cls, cls_args = self._typeinfo()
        args = list(cls_args) + list(parts)
        return cls(*args)

    def _merge_similar(self, parts):
        """Merge adjacent text objects with the same type and parameters together.

        >>> text = Text()
        >>> parts = [Tag('em', 'Breaking'), Tag('em', ' '), Tag('em', 'news!')]
        >>> list(text._merge_similar(parts))
        [Tag('em', 'Breaking news!')]
        """

        groups = itertools.groupby(parts, lambda value: value._typeinfo())
        for typeinfo, group in groups:
            cls, info = typeinfo
            group = list(group)
            if cls and len(group) > 1:
                group_parts = itertools.chain(*(text.parts for text in group))
                args = list(info) + list(group_parts)
                yield cls(*args)
            else:
                for text in group:
                    yield text

    @deprecated('0.19', 'use __unicode__() instead')
    def plaintext(self):
        return str(self)

    @deprecated('0.19')
    def enumerate(self):
        for n, child in enumerate(self.parts):
            try:
                for p in child.enumerate():
                    yield p
            except AttributeError:
                yield self, n

    @deprecated('0.19')
    def reversed(self):
        for n, child in reversed(list(enumerate(self.parts))):
            try:
                for p in child.reversed():
                    yield p
            except AttributeError:
                yield self, n

    @deprecated('0.19', 'use slicing instead')
    def get_beginning(self):
        try:
            l, i = next(self.enumerate())
        except StopIteration:
            pass
        else:
            return l.parts[i]

    @deprecated('0.19', 'use slicing instead')
    def get_end(self):
        try:
            l, i = next(self.reversed())
        except StopIteration:
            pass
        else:
            return l.parts[i]

    @deprecated('0.19', 'use slicing instead')
    def apply_to_start(self, f):
        return self.map(f, lambda index, length: index == 0)

    @deprecated('0.19', 'use slicing instead')
    def apply_to_end(self, f):
        return self.map(f, lambda index, length: index == length - 1)

    @deprecated('0.19')
    def map(self, f, condition=None):
        if condition is None:
            condition = lambda index, length: True

        def iter_map_with_condition():
            length = len(self)
            for index, child in enumerate(self.parts):
                if hasattr(child, 'map'):
                    yield child.map(f, condition) if condition(index, length) else child
                else:
                    yield f(child) if condition(index, length) else child
        return self._create_similar(iter_map_with_condition())


class String(BaseText):
    """
    A :py:class:`String` is a wrapper for a plain Python string.

    >>> from pybtex.richtext import String
    >>> print(String('Crime & Punishment').render_as('text'))
    Crime & Punishment
    >>> print(String('Crime & Punishment').render_as('html'))
    Crime &amp; Punishment

    :py:class:`String` supports the same methods as :py:class:`Text`.

    """

    def __init__(self, *parts):
        """
        All arguments must be plain unicode strings.
        Arguments are concatenated together.

        >>> print(str(String('November', ', ', 'December', '.')))
        November, December.
        """

        self.value = ''.join(parts)

    def __repr__(self):
        return str_repr(self.value)

    def __str__(self):
        return str(self.value)

    def __eq__(self, other):
        """
        Compare two :py:class:`.String` objects.


        """
        return type(other) == type(self) and self.value == other.value

    def __len__(self):
        return self.value.__len__()

    def __contains__(self, item):
        return self.value.__contains__(item)

    def __getitem__(self, index):
        return String(self.value.__getitem__(index))

    def __add__(self, other):
        return BaseText.__add__(self, other)

    def split(self, sep=None, keep_empty_parts=None):
        if keep_empty_parts is None:
            keep_empty_parts = sep is not None

        if sep is None:
            from .textutils import whitespace_re
            parts = whitespace_re.split(self.value)
        elif isinstance(sep, str):
            parts = self.value.split(sep)
        else:
            try:
                split_method = sep.split
            except AttributeError:
                raise TypeError('sep must be None, string or compiled regular expression')
            else:
                parts = split_method(self.value)
        return [String(part) for part in parts if part or keep_empty_parts]

    def startswith(self, prefix):
        """
        Return True if string starts with the prefix,
        otherwise return False.

        prefix can also be a tuple of suffixes to look for.
        """
        return self.value.startswith(prefix)

    def endswith(self, suffix):
        """
        Return True if the string ends with the specified suffix,
        otherwise return False.

        suffix can also be a tuple of suffixes to look for.
        return self.value.endswith(text)
        """
        return self.value.endswith(suffix)

    def isalpha(self):
        return self.value.isalpha()

    def lower(self):
        return String(self.value.lower())

    def upper(self):
        return String(self.value.upper())

    @property
    def parts(self):
        return [str(self)]

    def _typeinfo(self):
        return String, ()

    def render(self, backend):
        return backend.format_str(self.value)


class Text(BaseMultipartText):
    """
    The :py:class:`Text` class is the top level container that may contain
    :py:class:`String`, :py:class:`Tag` or :py:class:`HRef` objects.

    """

    def __repr__(self):
        return 'Text({})'.format(', '.join(repr(part) for part in self.parts))

    def _unpack(self):
        for part in self.parts:
            yield part

    @classmethod
    def from_latex(cls, latex):
        import codecs
        import latexcodec  # noqa
        from pybtex.markup import LaTeXParser

        return LaTeXParser(codecs.decode(latex, 'ulatex')).parse()


class Tag(BaseMultipartText):
    r"""
    A :py:class:`Tag` represents something like an HTML tag
    or a LaTeX formatting command:

    >>> from pybtex.richtext import Tag
    >>> tag = Tag('em', 'The TeXbook')
    >>> print(tag.render_as('html'))
    <em>The TeXbook</em>
    >>> print(tag.render_as('latex'))
    \emph{The TeXbook}

    :py:class:`Tag` supports the same methods as :py:class:`Text`.
    """

    def __check_name(self, name):
        depr_map = {}
        depr_map[u'emph'] = u'em'
        if name in depr_map:
            msg = u"The tag '%s' is deprecated" % name
            msg += u", use '%s' instead." % depr_map[name]
            warnings.warn(msg, DeprecationWarning, stacklevel=3)
            return depr_map[name]
        return name

    def __init__(self, name, *args):
        if not isinstance(name, (str, Text)):
            raise ValueError(
                "name must be str or Text (got %s)" % name.__class__.__name__)
        self.name = self.__check_name(str(name))
        self.info = self.name,
        super(Tag, self).__init__(*args)

    def __repr__(self):
        if self.parts:
            reprparts = ', '.join(repr(part) for part in self.parts)
            return 'Tag({}, {})'.format(str_repr(self.name), reprparts)
        else:
            return 'Tag({})'.format(str_repr(self.name))

    def render(self, backend):
        text = super(Tag, self).render(backend)
        return backend.format_tag(self.name, text)


class HRef(BaseMultipartText):
    """
    A :py:class:`HRef` represends a hyperlink:

    >>> from pybtex.richtext import Tag
    >>> href = HRef('http://ctan.org/', 'CTAN')
    >>> print(href.render_as('html'))
    <a href="http://ctan.org/">CTAN</a>
    >>> print(href.render_as('latex'))
    \\href{http://ctan.org/}{CTAN}

    >>> href = HRef(String('http://ctan.org/'), String('http://ctan.org/'))
    >>> print(href.render_as('latex'))
    \\url{http://ctan.org/}

    :py:class:`HRef` supports the same methods as :py:class:`Text`.

    """

    def __init__(self, url, *args, external=False):
        if not isinstance(url, (str, BaseText)):
            raise ValueError(
                "url must be str or Text (got %s)" % url.__class__.__name__)
        self.url = str(url)
        self.info = self.url,
        self.external = external
        super(HRef, self).__init__(*args)

    def __repr__(self):
        reprparts = ', '.join(repr(part) for part in self.parts)
        return 'HRef({}, {})'.format(str_repr(self.url), reprparts)

    def render(self, backend):
        text = super(HRef, self).render(backend)
        return backend.format_href(self.url, text, self.external)


class Protected(BaseMultipartText):
    r"""
    A :py:class:`Protected` represents a "protected" piece of text.

    - :py:meth:`Protected.lower`, :py:meth:`Protected.upper`,
      :py:meth:`Protected.capitalize`, and :py:meth:`Protected.capitalize()`
      are no-ops and just return the :py:class:`Protected` object itself.
    - :py:meth:`Protected.split` never splits the text. It always returns a
      one-element list containing the :py:class:`Protected` object itself.
    - In LaTeX output, :py:class:`Protected` is {surrounded by braces}.  HTML
      and plain text backends just output the text as-is.

    >>> from pybtex.richtext import Protected
    >>> text = Protected('The CTAN archive')
    >>> text.lower()
    Protected('The CTAN archive')
    >>> text.split()
    [Protected('The CTAN archive')]
    >>> print(text.render_as('latex'))
    {The CTAN archive}
    >>> print(text.render_as('html'))
    <span class="bibtex-protected">The CTAN archive</span>

    .. versionadded:: 0.20

    """

    def __init__(self, *args):
        super(Protected, self).__init__(*args)

    def __repr__(self):
        reprparts = ', '.join(repr(part) for part in self.parts)
        return 'Protected({})'.format(reprparts)

    def capfirst(self):
        return self

    def capitalize(self):
        return self

    def lower(self):
        return self

    def upper(self):
        return self

    def split(self, sep=None, keep_empty_parts=None):
        return [self]

    def render(self, backend):
        text = super(Protected, self).render(backend)
        return backend.format_protected(text)


class Symbol(BaseText):
    """A special symbol. This class is rarely used and may be removed in
    future versions.

    Examples of special symbols are non-breaking spaces and dashes.

    :py:class:`Symbol` supports the same methods as :py:class:`Text`.
    """

    def __init__(self, name):
        self.name = name
        self.info = self.name,

    def __len__(self):
        return 1

    def __repr__(self):
        return "Symbol(%s)" % str_repr(self.name)

    def __str__(self):
        # XXX
        return u'<%s>' % self.name

    def __eq__(self, other):
        return self.name == other.name

    def __contains__(self, item):
        return False

    def __getitem__(self, index):
        # mimic the behavior of a 1-element string
        try:
            result = 'a'[index]
        except IndexError:
            raise IndexError('richtext.Symbol index out of range')
        else:
            return self if result else String()

    def split(self, sep=None, keep_empty_parts=None):
        return [self]

    def startswith(self, text):
        return False

    def endswith(self, text):
        return False

    def isalpha(self):
        return False

    def render(self, backend):
        return backend.symbols[self.name]

    def upper(self):
        return self

    def lower(self):
        return self


nbsp = Symbol('nbsp')
