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

from __future__ import absolute_import, unicode_literals

import re

from pybtex.bibtex.exceptions import BibTeXError
from pybtex.utils import pairwise

whitespace_re = re.compile(r'(\s)')
purify_special_char_re = re.compile(r'^\\[A-Za-z]+')


def wrap(string, width=79, subsequent_indent='  '):
    r"""
    Wrap long string into multiple lines by inserting line breaks.

    The string is broken at whitespace characters so that each line is as long
    as possible, but no longer than ``width`` characters.

    If there are no possible break points in the first ``width`` characters, a
    longer line will be produced, with the line break inserted at the first
    possible whitespace characters after ``width``.

    After each line break, the subsequent line is indented with
    ``subsequent_indent`` (two spaces by default).

    The lines are not allowed to be shorter than ``len(subsequent_indent) + 1``
    (3 characters by default), so that each line contains at least one
    non-whitespace character after the indent.

    >>> print(wrap('', width=3))
    <BLANKLINE>
    >>> print(wrap('0123456789 12345', width=10))
    0123456789
      12345
    >>> print(wrap('01234 6789 12345', width=10))
    01234 6789
      12345
    >>> print(wrap('01234 6789 12345', width=11))
    01234 6789
      12345
    >>> print(wrap('01234 6789 12345', width=9))
    01234
      6789
      12345
    >>> print(wrap(' a b c', width=3))
     a b
      c
    >>> print(wrap('aa bb c', width=3))
    aa bb
      c

    """

    min_width = len(subsequent_indent)

    def find_break(string):
        for prev_match, match in pairwise(whitespace_re.finditer(string)):
            if (match is None or match.start() > width) and prev_match.start() > min_width:
                return prev_match.start()

    def iter_lines(string):
        while len(string) > width:
            break_pos = find_break(string)
            if not break_pos:
                yield string
                return
            yield string[:break_pos]
            string = subsequent_indent + string[break_pos + 1:]
        if string:
            yield string

    return '\n'.join(line.rstrip() for line in iter_lines(string))


class BibTeXString(object):
    def __init__(self, chars, level=0, max_level=100):
        if level > max_level:
            raise BibTeXError('too many nested braces')

        self.level = level
        self.is_closed = False
        self.contents = list(self.find_closing_brace(iter(chars)))

    def __iter__(self):
        return self.traverse()

    def find_closing_brace(self, chars):
        for char in chars:
            if char == '{':
                yield BibTeXString(chars, self.level + 1)
            elif char == '}' and self.level > 0:
                self.is_closed = True
                return
            else:
                yield char

    def is_special_char(self):
        return self.level == 1 and self.contents and self.contents[0] == '\\'

    def traverse(self, open=None, f=lambda char, string: char, close=None):
        if open is not None and self.level > 0:
            yield open(self)

        for child in self.contents:
            if hasattr(child, 'traverse'):
                if child.is_special_char():
                    if open is not None:
                        yield open(child)
                    yield f(child.inner_string(), child)
                    if close is not None:
                        yield close(child)
                else:
                    for result in child.traverse(open, f, close):
                        yield result
            else:
                yield f(child, self)

        if close is not None and self.level > 0 and self.is_closed:
            yield close(self)

    def __str__(self):
        return ''.join(self.traverse(open=lambda string: '{', close=lambda string: '}'))

    def inner_string(self):
        return ''.join(str(child) for child in self.contents)


def change_case(string, mode):
    r"""
    >>> print(change_case('aBcD', 'l'))
    abcd
    >>> print(change_case('aBcD', 'u'))
    ABCD
    >>> print(change_case('ABcD', 't'))
    Abcd
    >>> print(change_case(r'The {\TeX book \noop}', 'u'))
    THE {\TeX BOOK \noop}
    >>> print(change_case(r'And Now: BOOO!!!', 't'))
    And now: Booo!!!
    >>> print(change_case(r'And {Now: BOOO!!!}', 't'))
    And {Now: BOOO!!!}
    >>> print(change_case(r'And {Now: {BOOO}!!!}', 'l'))
    and {Now: {BOOO}!!!}
    >>> print(change_case(r'And {\Now: BOOO!!!}', 't'))
    And {\Now: booo!!!}
    >>> print(change_case(r'And {\Now: {BOOO}!!!}', 'l'))
    and {\Now: {booo}!!!}
    >>> print(change_case(r'{\TeX\ and databases\Dash\TeX DBI}', 't'))
    {\TeX\ and databases\Dash\TeX DBI}
    """

    def title(char, state):
        if state == 'start':
            return char
        else:
            return char.lower()

    lower = lambda char, state: char.lower()
    upper = lambda char, state: char.upper()

    convert = {'l': lower, 'u': upper, 't': title}[mode]

    def convert_special_char(special_char, state):
        # FIXME BibTeX treats some accented and foreign characterss specially
        def convert_words(words):
            for word in words:
                if word.startswith('\\'):
                    yield word
                else:
                    yield convert(word, state)

        return ' '.join(convert_words(special_char.split(' ')))

    def change_case_iter(string, mode):
        state = 'start'
        for char, brace_level in scan_bibtex_string(string):
            if brace_level == 0:
                yield convert(char, state)
                if char == ':':
                    state = 'after colon'
                elif char.isspace() and state == 'after colon':
                    state = 'start'
                else:
                    state = 'normal'
            else:
                if brace_level == 1 and char.startswith('\\'):
                    yield convert_special_char(char, state)
                else:
                    yield char

    return ''.join(change_case_iter(string, mode))


def bibtex_substring(string, start, length):
    r"""
    Return a substring of the given length, starting from the given position.

    start and length are 1-based. If start is < 0, it is counted from the end
    of the string. If start is 0, an empty string is returned.

    >>> print(bibtex_substring('abcdef', 1, 3))
    abc
    >>> print(bibtex_substring('abcdef', 2, 3))
    bcd
    >>> print(bibtex_substring('abcdef', 2, 1000))
    bcdef
    >>> print(bibtex_substring('abcdef', 0, 1000))
    <BLANKLINE>
    >>> print(bibtex_substring('abcdef', -1, 1))
    f
    >>> print(bibtex_substring('abcdef', -1, 2))
    ef
    >>> print(bibtex_substring('abcdef', -2, 3))
    cde
    >>> print(bibtex_substring('abcdef', -2, 1000))
    abcde
    """

    if start > 0:
        start0 = start - 1
        end0 = start0 + length
    elif start < 0:
        end0 = len(string) + start + 1
        start0 = end0 - length
    else: # start == 0:
        return u''
    return string[start0:end0]


def bibtex_len(string):
    r"""Return the number of characters in the string.

    Braces are ignored. "Special characters" are ignored. A "special character"
    is a substring at brace level 1, if the first character after the opening
    brace is a backslash, like in "de la Vall{\'e}e Poussin".

    >>> print(bibtex_len(r"de la Vall{\'e}e Poussin"))
    20
    >>> print(bibtex_len(r"de la Vall{e}e Poussin"))
    20
    >>> print(bibtex_len(r"de la Vallee Poussin"))
    20
    >>> print(bibtex_len(r'\ABC 123'))
    8
    >>> print(bibtex_len(r'{\abc}'))
    1
    >>> print(bibtex_len(r'{\abc'))
    1
    >>> print(bibtex_len(r'}\abc'))
    4
    >>> print(bibtex_len(r'\abc}'))
    4
    >>> print(bibtex_len(r'\abc{'))
    4
    >>> print(bibtex_len(r'level 0 {1 {2}}'))
    11
    >>> print(bibtex_len(r'level 0 {\1 {2}}'))
    9
    >>> print(bibtex_len(r'level 0 {1 {\2}}'))
    12
    """
    length = 0
    for char, brace_level in scan_bibtex_string(string):
        if char not in '{}':
            length += 1
    return length


def bibtex_width(string):
    r"""
    Determine the width of the given string, in relative units.

    >>> bibtex_width('')
    0
    >>> bibtex_width('abc')
    1500
    >>> bibtex_width('ab{c}')
    2500
    >>> bibtex_width(r"ab{\'c}")
    1500
    >>> bibtex_width(r"ab{\'c{}}")
    1500
    >>> bibtex_width(r"ab{\'c{}")
    1500
    >>> bibtex_width(r"ab{\'c{d}}")
    2056
    """

    from pybtex.charwidths import charwidths
    width = 0
    for token, brace_level in scan_bibtex_string(string):
        if brace_level == 1 and token.startswith('\\'):
            for char in token[2:]:
                if char not in '{}':
                    width += charwidths.get(char, 0)
            width -= 1000  # two braces
        else:
            width += charwidths.get(token, 0)
    return width


def bibtex_prefix(string, num_chars):
    r"""Return the firxt num_char characters of the string.

    Braces and "special characters" are ignored, as in bibtex_len.  If the
    resulting prefix ends at brace level > 0, missing closing braces are
    appended.

    >>> print(bibtex_prefix('abc', 1))
    a
    >>> print(bibtex_prefix('abc', 5))
    abc
    >>> print(bibtex_prefix('ab{c}d', 3))
    ab{c}
    >>> print(bibtex_prefix('ab{cd}', 3))
    ab{c}
    >>> print(bibtex_prefix('ab{cd', 3))
    ab{c}
    >>> print(bibtex_prefix(r'ab{\cd}', 3))
    ab{\cd}
    >>> print(bibtex_prefix(r'ab{\cd', 3))
    ab{\cd}

    """
    def prefix():
        length = 0
        for char, brace_level in scan_bibtex_string(string):
            yield char
            if char not in '{}':
                length += 1
            if length >= num_chars:
                break
        for i in range(brace_level):
            yield '}'
    return ''.join(prefix())


def bibtex_purify(string):
    r"""Strip special characters from the string.

    >>> print(bibtex_purify('Abc 1234'))
    Abc 1234
    >>> print(bibtex_purify('Abc  1234'))
    Abc  1234
    >>> print(bibtex_purify('Abc-Def'))
    Abc Def
    >>> print(bibtex_purify('Abc-~-Def'))
    Abc   Def
    >>> print(bibtex_purify('{XXX YYY}'))
    XXX YYY
    >>> print(bibtex_purify('{XXX {YYY}}'))
    XXX YYY
    >>> print(bibtex_purify(r'XXX {\YYY} XXX'))
    XXX  XXX
    >>> print(bibtex_purify(r'{XXX {\YYY} XXX}'))
    XXX YYY XXX
    >>> print(bibtex_purify(r'\\abc def'))
    abc def
    >>> print(bibtex_purify('a@#$@#$b@#$@#$c'))
    abc
    >>> print(bibtex_purify(r'{\noopsort{1973b}}1973'))
    1973b1973
    >>> print(bibtex_purify(r'{sort{1973b}}1973'))
    sort1973b1973
    >>> print(bibtex_purify(r'{sort{\abc1973b}}1973'))
    sortabc1973b1973
    >>> print(bibtex_purify(r'{\noopsort{1973a}}{\switchargs{--90}{1968}}'))
    1973a901968
    """

    # FIXME BibTeX treats some accented and foreign characterss specially
    def purify_iter(string):
        for token, brace_level in scan_bibtex_string(string):
            if brace_level == 1 and token.startswith('\\'):
                for char in purify_special_char_re.sub('', token):
                    if char.isalnum():
                        yield char
            else:
                if token.isalnum():
                    yield token
                elif token.isspace() or token in '-~':
                    yield ' '

    return ''.join(purify_iter(string))


def scan_bibtex_string(string):
    """ Yield (char, brace_level) tuples.

    "Special characters", as in bibtex_len, are treated as a single character

    """
    return BibTeXString(string).traverse(
        open=lambda string: ('{', string.level),
        f=lambda char, string: (char, string.level),
        close=lambda string: ('}', string.level - 1),
    )


def split_name_list(string):
    r"""
    Split a list of names, separated by ' and '.

    >>> split_name_list('Johnson and Peterson')
    ['Johnson', 'Peterson']
    >>> split_name_list('Johnson AND Peterson')
    ['Johnson', 'Peterson']
    >>> split_name_list('Johnson AnD Peterson')
    ['Johnson', 'Peterson']
    >>> split_name_list('Armand and Peterson')
    ['Armand', 'Peterson']
    >>> split_name_list('Armand and anderssen')
    ['Armand', 'anderssen']
    >>> split_name_list('{Armand and Anderssen}')
    ['{Armand and Anderssen}']
    >>> split_name_list('What a Strange{ }and Bizzare Name! and Peterson')
    ['What a Strange{ }and Bizzare Name!', 'Peterson']
    >>> split_name_list('What a Strange and{ }Bizzare Name! and Peterson')
    ['What a Strange and{ }Bizzare Name!', 'Peterson']
    """
    return split_tex_string(string, ' [Aa][Nn][Dd] ')


def _find_closing_brace(string):
    r"""
    >>> _find_closing_brace('')
    ('', '')
    >>> _find_closing_brace('no braces')
    ('no braces', '')
    >>> _find_closing_brace('brace at the end}')
    ('brace at the end}', '')
    >>> _find_closing_brace('two closing braces}}')
    ('two closing braces}', '}')
    >>> _find_closing_brace('two closing} braces} and some text')
    ('two closing}', ' braces} and some text')
    >>> _find_closing_brace('more {nested{}}{braces}} and the rest}')
    ('more {nested{}}{braces}}', ' and the rest}')
    """
    up_to_brace = []
    brace_level = 1
    while brace_level >= 1:
        next_brace = BRACE_RE.search(string)
        if not next_brace:
            break

        up_to_brace.append(string[:next_brace.end()])
        string = string[next_brace.end():]

        if next_brace.group() == '{':
            brace_level += 1
        elif next_brace.group() == '}':
            brace_level -= 1
        else:
            raise ValueError(next_brace.group())

    if not up_to_brace:
        up_to_brace, string = [string], ''
    return ''.join(up_to_brace), string


# "\ " is a "control space" in TeX, i. e. "a space that is not to be ignored"
#     -- The TeXbook, Chapter 3: Controlling TeX, p 8
# ~ is a space character, according to BibTeX
# \~ is not a space character
BIBTEX_SPACE_RE = re.compile(r'(?:\\ |\s|(?<!\\)~)+')
BRACE_RE = re.compile(r'{|}')


def split_tex_string(string, sep=None, strip=True, filter_empty=False):
    r"""Split a string using the given separator (regexp).

    Everything at brace level > 0 is ignored.

    >>> split_tex_string('')
    []
    >>> split_tex_string('     ')
    []
    >>> split_tex_string('.a.b.c.', r'\.')
    ['', 'a', 'b', 'c', '']
    >>> split_tex_string('.a.b.c.{d.}.', r'\.')
    ['', 'a', 'b', 'c', '{d.}', '']
    >>> split_tex_string('Matsui      Fuuka')
    ['Matsui', 'Fuuka']
    >>> split_tex_string('{Matsui      Fuuka}')
    ['{Matsui      Fuuka}']
    >>> split_tex_string(r'Matsui\ Fuuka')
    ['Matsui', 'Fuuka']
    >>> split_tex_string(r'{Matsui\ Fuuka}')
    ['{Matsui\\ Fuuka}']
    >>> split_tex_string('a')
    ['a']
    >>> split_tex_string('on a')
    ['on', 'a']
    >>> split_tex_string(r'Qui\~{n}onero-Candela, J.')
    ['Qui\\~{n}onero-Candela,', 'J.']
    """

    if sep is None:
        sep = BIBTEX_SPACE_RE
        filter_empty = True

    sep = re.compile(sep)

    result = []
    word_parts = []

    while True:
        head, brace, string = string.partition('{')

        if head:
            head_parts = sep.split(head)
            for word in head_parts[:-1]:
                result.append(''.join(word_parts + [word]))
                word_parts = []
            word_parts.append(head_parts[-1])

        if brace:
            word_parts.append(brace)
            up_to_closing_brace, string = _find_closing_brace(string)
            word_parts.append(up_to_closing_brace)
        else:
            break

    if word_parts:
        result.append(''.join(word_parts))

    if strip:
        result = [part.strip() for part in result]
    if filter_empty:
        result = [part for part in result if part]
    return result


def bibtex_first_letter(string):
    r""" Return the first letter or special character of the string.

    >>> print(bibtex_first_letter('Andrew Blake'))
    A
    >>> print(bibtex_first_letter('{Andrew} Blake'))
    A
    >>> print(bibtex_first_letter('1Andrew'))
    A
    >>> print(bibtex_first_letter(r'{\TeX} markup'))
    {\TeX}
    >>> print(bibtex_first_letter(''))
    <BLANKLINE>
    >>> print(bibtex_first_letter('123 123 123 {}'))
    <BLANKLINE>
    >>> print(bibtex_first_letter(r'\LaTeX Project Team'))
    L

    """

    for char in BibTeXString(string):
        if char.startswith('\\') and char != '\\':
            return '{{{0}}}'.format(char)
        elif char.isalpha():
            return char
    return ''


def bibtex_abbreviate(string, delimiter=None, separator='-'):
    r"""
    Abbreviate string.

    >>> print(bibtex_abbreviate('Andrew Blake'))
    A
    >>> print(bibtex_abbreviate('Jean-Pierre'))
    J.-P
    >>> print(bibtex_abbreviate('Jean--Pierre'))
    J.-P
    
    """

    def _bibtex_abbreviate():
        for token in split_tex_string(string, sep=separator):
            letter = bibtex_first_letter(token)
            if letter:
                yield letter

    if delimiter is None:
        delimiter = '.-'
    return delimiter.join(_bibtex_abbreviate())
