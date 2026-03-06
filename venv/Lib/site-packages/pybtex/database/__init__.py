# vim: fileencoding=utf-8
# Copyright (c) 2006-2021  Andrey Golovizin
# Copyright (c)      2022  Julian Rüth
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

from __future__ import unicode_literals
from __future__ import print_function
import re
from typing import Iterable, Tuple

try:
    from collections.abc import Mapping
except ImportError:
    from collections import Mapping

import textwrap

from pybtex.exceptions import PybtexError
from pybtex.utils import (
    deprecated,
    OrderedCaseInsensitiveDict, CaseInsensitiveDefaultDict, CaseInsensitiveSet
)
from pybtex.richtext import Text
from pybtex.bibtex.utils import split_tex_string, scan_bibtex_string
from pybtex.errors import report_error
from pybtex.plugin import find_plugin


# for python2 compatibility
def indent(text, prefix):
    if hasattr(textwrap, "indent"):
        return textwrap.indent(text, prefix)
    else:
        return ''.join(prefix + line for line in text.splitlines(True))


class BibliographyDataError(PybtexError):
    pass


class InvalidNameString(PybtexError):
    def __init__(self, name_string):
        message = 'Too many commas in {}'.format(repr(name_string))
        super(InvalidNameString, self).__init__(message)


class BibliographyData(object):
    def __init__(self, entries=None, preamble=None, wanted_entries=None, min_crossrefs=2):
        """
        A :py:class:`.BibliographyData` object contains a dictionary of bibliography
        entries referenced by their keys.
        Each entry represented by an :py:class:`.Entry` object.

        Additionally, :py:class:`.BibliographyData` may contain a LaTeX
        preamble defined by ``@PREAMBLE`` commands in the BibTeX file.
        """

        self.entries = OrderedCaseInsensitiveDict()
        '''A dictionary of bibliography entries referenced by their keys.

        The dictionary is case insensitive:

        >>> bib_data = parse_string("""
        ...     @ARTICLE{gnats,
        ...         author = {L[eslie] A. Aamport},
        ...         title = {The Gnats and Gnus Document Preparation System},
        ...     }
        ... """, 'bibtex')
        >>> bib_data.entries['gnats'] == bib_data.entries['GNATS']
        True

        '''

        self.crossref_count = CaseInsensitiveDefaultDict(int)
        self.min_crossrefs = min_crossrefs
        self._preamble = []
        if wanted_entries is not None:
            self.wanted_entries = CaseInsensitiveSet(wanted_entries)
            self.citations = CaseInsensitiveSet(wanted_entries)
        else:
            self.wanted_entries = None
            self.citations = CaseInsensitiveSet()
        if entries:
            if isinstance(entries, Mapping):
                entries = entries.items()
            for (key, entry) in entries:
                self.add_entry(key, entry)
        if preamble:
            self._preamble.extend(preamble)

    def __eq__(self, other):
        if not isinstance(other, BibliographyData):
            return super(BibliographyData, self) == other
        return (
            self.entries == other.entries
            and self._preamble == other._preamble
        )

    def __repr__(self):

        repr_entry = repr(self.entries)
        keys = self.entries.keys()

        for key in keys:
            ind = repr_entry.index(key) - 2  # find first instance
            repr_entry = repr_entry[:ind] + "\n" + repr_entry[ind:]

        repr_entry = indent(repr_entry, prefix="    ")
        repr_entry = repr_entry[4:]  # drop 1st indent

        return (
            "BibliographyData(\n"
            "  entries={0},\n\n"
            "  preamble={1})".format(repr_entry, repr(self._preamble))
        )

    def add_to_preamble(self, *values):
        self._preamble.extend(values)

    @property
    def preamble_list(self):
        r'''
        LaTeX preamble as list of strings

        >>> bib_data = parse_string(r"""
        ...     @PREAMBLE{"\newcommand{\noopsort}[1]{}"}
        ...     @PREAMBLE{"\newcommand{\nooptilde}[1]{}"}
        ... """, 'bibtex')
        >>> print(bib_data.preamble_list)
        ['\\newcommand{\\noopsort}[1]{}', '\\newcommand{\\nooptilde}[1]{}']

        .. versionadded:: 0.19
            Earlier versions used :py:meth:`.get_preamble()`, which is now deprecated.
        '''
        return self._preamble

    @property
    def preamble(self):
        r'''
        LaTeX preamble.

        >>> bib_data = parse_string(r"""
        ...     @PREAMBLE{"\newcommand{\noopsort}[1]{}"}
        ... """, 'bibtex')
        >>> print(bib_data.preamble)
        \newcommand{\noopsort}[1]{}

        .. versionadded:: 0.19
            Earlier versions used :py:meth:`.get_preamble()`, which is now deprecated.
        '''
        return ''.join(self._preamble)

    @deprecated('0.19', 'use BibliographyData.preamble instead')
    def get_preamble(self):
        """
        .. deprecated:: 0.19
            Use :py:attr:`.preamble` instead.
        """
        return self.preamble

    def want_entry(self, key):
        return (
            self.wanted_entries is None
            or key in self.wanted_entries
            or '*' in self.wanted_entries
        )

    def get_canonical_key(self, key):
        if key in self.citations:
            return self.citations.get_canonical_key(key)
        else:
            return key

    def add_entry(self, key, entry):
        if not self.want_entry(key):
            return
        if key in self.entries:
            report_error(BibliographyDataError('repeated bibliography entry: %s' % key))
            return
        entry.key = self.get_canonical_key(key)
        self.entries[entry.key] = entry
        try:
            crossref = entry.fields['crossref']
        except KeyError:
            pass
        else:
            if self.wanted_entries is not None:
                self.wanted_entries.add(crossref)

    def add_entries(self, entries: Iterable[Tuple[str, "Entry"]]) -> None:
        for key, entry in entries:
            self.add_entry(key, entry)

    def _get_crossreferenced_citations(self, citations, min_crossrefs):
        r"""
        Get cititations not cited explicitly but referenced by other citations.

        >>> from pybtex.database import Entry
        >>> data = BibliographyData({
        ...     'main_article': Entry('article', {'crossref': 'xrefd_arcicle'}),
        ...     'xrefd_arcicle': Entry('article'),
        ... })
        >>> list(data._get_crossreferenced_citations([], min_crossrefs=1))
        []
        >>> list(data._get_crossreferenced_citations(['main_article'], min_crossrefs=1))
        ['xrefd_arcicle']
        >>> list(data._get_crossreferenced_citations(['Main_article'], min_crossrefs=1))
        ['xrefd_arcicle']
        >>> list(data._get_crossreferenced_citations(['main_article'], min_crossrefs=2))
        []
        >>> list(data._get_crossreferenced_citations(['xrefd_arcicle'], min_crossrefs=1))
        []

        >>> data2 = BibliographyData(data.entries, wanted_entries=data.entries.keys())
        >>> list(data2._get_crossreferenced_citations([], min_crossrefs=1))
        []
        >>> list(data2._get_crossreferenced_citations(['main_article'], min_crossrefs=1))
        ['xrefd_arcicle']
        >>> list(data2._get_crossreferenced_citations(['Main_article'], min_crossrefs=1))
        ['xrefd_arcicle']
        >>> list(data2._get_crossreferenced_citations(['main_article'], min_crossrefs=2))
        []
        >>> list(data2._get_crossreferenced_citations(['xrefd_arcicle'], min_crossrefs=1))
        []
        >>> list(data2._get_crossreferenced_citations(['xrefd_arcicle'], min_crossrefs=1))
        []

        """

        crossref_count = CaseInsensitiveDefaultDict(int)
        citation_set = CaseInsensitiveSet(citations)
        for citation in citations:
            try:
                entry = self.entries[citation]
                crossref = entry.fields['crossref']
            except KeyError:
                continue
            try:
                crossref_entry = self.entries[crossref]
            except KeyError:
                report_error(BibliographyDataError(
                    'bad cross-reference: entry "{key}" refers to '
                    'entry "{crossref}" which does not exist.'.format(
                        key=citation, crossref=crossref,
                    )
                ))
                continue

            canonical_crossref = crossref_entry.key
            crossref_count[canonical_crossref] += 1
            if crossref_count[canonical_crossref] >= min_crossrefs and canonical_crossref not in citation_set:
                citation_set.add(canonical_crossref)
                yield canonical_crossref

    def _expand_wildcard_citations(self, citations):
        r"""
        Expand wildcard citations (\citation{*} in .aux file).

        >>> from pybtex.database import Entry
        >>> data = BibliographyData((
        ...     ('uno', Entry('article')),
        ...     ('dos', Entry('article')),
        ...     ('tres', Entry('article')),
        ...     ('cuatro', Entry('article')),
        ... ))
        >>> list(data._expand_wildcard_citations([]))
        []
        >>> list(data._expand_wildcard_citations(['*']))
        ['uno', 'dos', 'tres', 'cuatro']
        >>> list(data._expand_wildcard_citations(['uno', '*']))
        ['uno', 'dos', 'tres', 'cuatro']
        >>> list(data._expand_wildcard_citations(['dos', '*']))
        ['dos', 'uno', 'tres', 'cuatro']
        >>> list(data._expand_wildcard_citations(['*', 'uno']))
        ['uno', 'dos', 'tres', 'cuatro']
        >>> list(data._expand_wildcard_citations(['*', 'DOS']))
        ['uno', 'dos', 'tres', 'cuatro']

        """

        citation_set = CaseInsensitiveSet()
        for citation in citations:
            if citation == '*':
                for key in self.entries:
                    if key not in citation_set:
                        citation_set.add(key)
                        yield key
            else:
                if citation not in citation_set:
                    citation_set.add(citation)
                    yield citation

    def add_extra_citations(self, citations, min_crossrefs):
        expanded_citations = list(self._expand_wildcard_citations(citations))
        crossrefs = list(self._get_crossreferenced_citations(expanded_citations, min_crossrefs))
        return expanded_citations + crossrefs

    def to_string(self, bib_format, **kwargs):
        """
        Return the data as a unicode string in the given format.

        :param bib_format: Data format ("bibtex", "yaml", etc.).

        .. versionadded:: 0.19
        """
        writer = find_plugin('pybtex.database.output', bib_format)(**kwargs)
        return writer.to_string(self)

    @classmethod
    def from_string(cls, value, bib_format, **kwargs):
        """
        Return the data from a unicode string in the given format.

        :param bib_format: Data format ("bibtex", "yaml", etc.).

        .. versionadded:: 0.22.2
        """
        return parse_string(value, bib_format, **kwargs)

    def to_bytes(self, bib_format, **kwargs):
        """
        Return the data as a byte string in the given format.

        :param bib_format: Data format ("bibtex", "yaml", etc.).

        .. versionadded:: 0.19
        """
        writer = find_plugin('pybtex.database.output', bib_format)(**kwargs)
        return writer.to_bytes(self)

    def to_file(self, file, bib_format=None, **kwargs):
        """
        Save the data to a file.

        :param file: A file name or a file-like object.
        :param bib_format: Data format ("bibtex", "yaml", etc.).
            If not specified, Pybtex will try to guess by the file name.

        .. versionadded:: 0.19
        """
        if isinstance(file, str):
            filename = file
        else:
            filename = getattr(file, 'name', None)
        writer = find_plugin('pybtex.database.output', bib_format, filename=filename)(**kwargs)
        return writer.write_file(self, file)

    def lower(self):
        '''
        Return another :py:class:`.BibliographyData` with all identifiers converted to lowercase.

        >>> data = parse_string("""
        ...     @BOOK{Obrazy,
        ...         title = "Obrazy z Rus",
        ...         author = "Karel Havlíček Borovský",
        ...     }
        ...     @BOOK{Elegie,
        ...         title = "Tirolské elegie",
        ...         author = "Karel Havlíček Borovský",
        ...     }
        ... """, 'bibtex')
        >>> data_lower = data.lower()
        >>> list(data_lower.entries.keys())
        ['obrazy', 'elegie']
        >>> for entry in data_lower.entries.values():
        ...     entry.key
        ...     list(entry.persons.keys())
        ...     list(entry.fields.keys())
        'obrazy'
        ['author']
        ['title']
        'elegie'
        ['author']
        ['title']

        '''

        entries_lower = ((key.lower(), entry.lower()) for key, entry in self.entries.items())
        return type(self)(
            entries=entries_lower,
            preamble=self._preamble,
            wanted_entries=self.wanted_entries,
            min_crossrefs=self.min_crossrefs,
        )


class RichFieldProxyDict(Mapping):
    def __init__(self, fields):
        self._fields = fields

    def __contains__(self):
        return self._fields.__contains__()

    def __iter__(self):
        return self._fields.__iter__()

    def __len__(self):
        return self._fields.__len__()

    def __getitem__(self, key):
        return Text.from_latex(self._fields[key])


class Entry(object):
    """A bibliography entry."""

    type = None
    """Entry type (``'book'``, ``'article'``, etc.)."""

    key = None
    """Entry key (for example, ``'fukushima1980neocognitron'``)."""

    fields = None
    """A dictionary of entry fields.
    The dictionary is ordered and case-insensitive."""

    persons = None
    """
    A dictionary of entry persons, by their roles.

    The most often used roles are ``'author'`` and ``'editor'``.
    """

    """A reference to the containing :py:class:`.BibliographyData` object. Used to resolve crossrefs."""

    def __init__(self, type_, fields=None, persons=None):
        if fields is None:
            fields = {}
        if persons is None:
            persons = {}
        self.type = type_.lower()
        self.original_type = type_
        self.fields = OrderedCaseInsensitiveDict(fields)
        self.persons = OrderedCaseInsensitiveDict(persons)

    def __eq__(self, other):
        if not isinstance(other, Entry):
            return super(Entry, self) == other
        return (
            self.type == other.type
            and self.fields == other.fields
            and self.persons == other.persons
        )

    def __repr__(self):
        r"""
        Return a string representation of this entry.

        >>> from pybtex.database import Entry
        >>> Entry("article")
        Entry('article')

        >>> Entry("article", fields={"title": "An article about journals", "journal": "Royal Society of Journals"}, persons={"editor": Person("Donald E. Knuth")})
        Entry('article',
          fields=[
            ('title', 'An article about journals'),
            ('journal', 'Royal Society of Journals')],
          persons={'editor': Person('Knuth, Donald E.')})

        """
        parameters = [repr(self.type)]

        if self.fields:
            parameters.append("fields=[\n    {}]".format(",\n    ".join(
                ["({!r}, {!r})".format(key, value) for (key, value) in self.fields.items()])))

        if self.persons:
            parameters.append("persons={!r}".format(dict(self.persons)))

        return "Entry({})".format(",\n  ".join(parameters))

    def add_person(self, person, role):
        self.persons.setdefault(role, []).append(person)

    def lower(self):
        return type(self)(
            self.type,
            fields=self.fields.lower(),
            persons=self.persons.lower(),
        )

    def _find_person_field(self, role):
        persons = self.persons[role]
        return ' and '.join(str(person) for person in persons)

    def _find_crossref_field(self, name, bib_data):
        if bib_data is None or 'crossref' not in self.fields:
            raise KeyError(name)
        referenced_entry = bib_data.entries[self.fields['crossref']]
        return referenced_entry._find_field(name, bib_data)

    def _find_field(self, name, bib_data=None):
        """
        Find the field with the given ``name`` according to this rules:

        - If the given field ``name`` in in ``self.fields``, just return
          self.fields[name].

        - Otherwise, if ``name`` is ``"authors"`` or ``"editors"`` (or any other
          person role), return the list of names as a string, separated by
          ``" and "``.

        - Otherwise, if this entry has a ``crossreff`` field, look up for the
          cross-referenced entry and try to find its field with the given
          ``name``.
        """
        try:
            return self.fields[name]
        except KeyError:
            try:
                return self._find_person_field(name)
            except KeyError:
                return self._find_crossref_field(name, bib_data)

    def to_string(self, bib_format, **kwargs):
        """
        Return the data as a unicode string in the given format.

        :param bib_format: Data format ("bibtex", "yaml", etc.).

        """
        writer = find_plugin('pybtex.database.output', bib_format)(**kwargs)
        return writer.to_string(BibliographyData(entries={self.key: self}))

    @classmethod
    def from_string(cls, value, bib_format, entry_number=0, **kwargs):
        """
        Return the data from a unicode string in the given format.

        :param bib_format: Data format ("bibtex", "yaml", etc.).
        :param entry_number: entry number if the string has more than one.

        .. versionadded:: 0.22.2
        """
        # get bibliography
        bibdata = BibliographyData.from_string(value, bib_format, **kwargs)
        # grab specific instance
        key = tuple(bibdata.entries.keys())[entry_number]
        return bibdata.entries[key]


class Person(object):
    """A person or some other person-like entity.

    >>> knuth = Person('Donald E. Knuth')
    >>> knuth.first_names
    ['Donald']
    >>> knuth.middle_names
    ['E.']
    >>> knuth.last_names
    ['Knuth']

    """

    first_names = None
    """
    A list of first names.

    .. versionadded:: 0.19
        Earlier versions used :py:meth:`.first`, which is now deprecated.
    """

    middle_names = None
    """
    A list of middle names.

    .. versionadded:: 0.19
        Earlier versions used :py:meth:`.middle`, which is now deprecated.
    """

    prelast_names = None
    """
    A list of pre-last (aka von) name parts.

    .. versionadded:: 0.19
        Earlier versions used :py:meth:`.middle`, which is now deprecated.
    """

    last_names = None
    """
    A list of last names.

    .. versionadded:: 0.19
        Earlier versions used :py:meth:`.last`, which is now deprecated.
    """

    lineage_names = None
    """
    A list of linage (aka Jr) name parts.

    .. versionadded:: 0.19
        Earlier versions used :py:meth:`.lineage`, which is now deprecated.
    """

    valid_roles = ['author', 'editor']
    style1_re = re.compile(r'^(.+),\s*(.+)$')
    style2_re = re.compile(r'^(.+),\s*(.+),\s*(.+)$')

    def __init__(self, string="", first="", middle="", prelast="", last="", lineage=""):
        """
        :param string: The full name string.
            It will be parsed and split into separate first, last, middle,
            pre-last and lineage name parst.

            Supported name formats are:

            - von Last, First
            - von Last, Jr, First
            - First von Last

            (see BibTeX manual for explanation)

        """

        self.first_names = []
        self.middle_names = []
        self.prelast_names = []
        self.last_names = []
        self.lineage_names = []

        string = string.strip()
        if string:
            self._parse_string(string)
        self.first_names.extend(split_tex_string(first))
        self.middle_names.extend(split_tex_string(middle))
        self.prelast_names.extend(split_tex_string(prelast))
        self.last_names.extend(split_tex_string(last))
        self.lineage_names.extend(split_tex_string(lineage))

    @property
    def bibtex_first_names(self):
        """A list of first and middle names together.
        (BibTeX treats all middle names as first.)

        .. versionadded:: 0.19
            Earlier versions used :py:meth:`Person.bibtex_first`, which is now deprecated.


        >>> knuth = Person('Donald E. Knuth')
        >>> knuth.bibtex_first_names
        ['Donald', 'E.']
        """
        return self.first_names + self.middle_names

    def _parse_string(self, name):
        """Extract various parts of the name from a string.

        >>> p = Person('Avinash K. Dixit')
        >>> print(p.first_names)
        ['Avinash']
        >>> print(p.middle_names)
        ['K.']
        >>> print(p.prelast_names)
        []
        >>> print(p.last_names)
        ['Dixit']
        >>> print(p.lineage_names)
        []
        >>> print(str(p))
        Dixit, Avinash K.
        >>> p == Person(str(p))
        True
        >>> p = Person('Dixit, Jr, Avinash K. ')
        >>> print(p.first_names)
        ['Avinash']
        >>> print(p.middle_names)
        ['K.']
        >>> print(p.prelast_names)
        []
        >>> print(p.last_names)
        ['Dixit']
        >>> print(p.lineage_names)
        ['Jr']
        >>> print(str(p))
        Dixit, Jr, Avinash K.
        >>> p == Person(str(p))
        True

        >>> p = Person('abc')
        >>> print(p.first_names, p.middle_names, p.prelast_names, p.last_names, p.lineage_names)
        [] [] [] ['abc'] []
        >>> p = Person('Viktorov, Michail~Markovitch')
        >>> print(p.first_names, p.middle_names, p.prelast_names, p.last_names, p.lineage_names)
        ['Michail'] ['Markovitch'] [] ['Viktorov'] []
        """
        def process_first_middle(parts):
            try:
                self.first_names.append(parts[0])
                self.middle_names.extend(parts[1:])
            except IndexError:
                pass

        def process_von_last(parts):
            # von cannot be the last name in the list
            von_last = parts[:-1]
            definitely_not_von = parts[-1:]

            if von_last:
                von, last = rsplit_at(von_last, is_von_name)
                self.prelast_names.extend(von)
                self.last_names.extend(last)
            self.last_names.extend(definitely_not_von)

        def find_pos(lst, pred):
            for i, item in enumerate(lst):
                if pred(item):
                    return i
            return i + 1

        def split_at(lst, pred):
            """Split the given list into two parts.

            The second part starts with the first item for which the given
            predicate is True.
            """
            pos = find_pos(lst, pred)
            return lst[:pos], lst[pos:]

        def rsplit_at(lst, pred):
            rpos = find_pos(reversed(lst), pred)
            pos = len(lst) - rpos
            return lst[:pos], lst[pos:]

        def is_von_name(string):
            if string[0].isupper():
                return False
            if string[0].islower():
                return True
            else:
                for char, brace_level in scan_bibtex_string(string):
                    if brace_level == 0 and char.isalpha():
                        return char.islower()
                    elif brace_level == 1 and char.startswith('\\'):
                        return special_char_islower(char)
            return False

        def special_char_islower(special_char):
            control_sequence = True
            for char in special_char[1:]:  # skip the backslash
                if control_sequence:
                    if not char.isalpha():
                        control_sequence = False
                else:
                    if char.isalpha():
                        return char.islower()
            return False

        parts = split_tex_string(name, ',')
        if len(parts) > 3:
            report_error(InvalidNameString(name))
            last_parts = parts[2:]
            parts = parts[:2] + [' '.join(last_parts)]

        if len(parts) == 3:  # von Last, Jr, First
            process_von_last(split_tex_string(parts[0]))
            self.lineage_names.extend(split_tex_string(parts[1]))
            process_first_middle(split_tex_string(parts[2]))
        elif len(parts) == 2:  # von Last, First
            process_von_last(split_tex_string(parts[0]))
            process_first_middle(split_tex_string(parts[1]))
        elif len(parts) == 1:  # First von Last
            parts = split_tex_string(name)
            first_middle, von_last = split_at(parts, is_von_name)
            if not von_last and first_middle:
                last = first_middle.pop()
                von_last.append(last)
            process_first_middle(first_middle)
            process_von_last(von_last)
        else:
            # should hot really happen
            raise ValueError(name)

    def __eq__(self, other):
        if not isinstance(other, Person):
            return super(Person, self) == other
        return (
            self.first_names == other.first_names
            and self.middle_names == other.middle_names
            and self.prelast_names == other.prelast_names
            and self.last_names == other.last_names
            and self.lineage_names == other.lineage_names
        )

    def __str__(self):
        # von Last, Jr, First
        von_last = ' '.join(self.prelast_names + self.last_names)
        jr = ' '.join(self.lineage_names)
        first = ' '.join(self.first_names + self.middle_names)
        return ', '.join(part for part in (von_last, jr, first) if part)

    def __repr__(self):
        return 'Person({0})'.format(repr(str(self)))

    def get_part_as_text(self, type):
        names = getattr(self, type + '_names')
        return ' '.join(names)

    def get_part(self, type, abbr=False):
        """Get a list of name parts by `type`.

        >>> knuth = Person('Donald E. Knuth')
        >>> knuth.get_part('first')
        ['Donald']
        >>> knuth.get_part('last')
        ['Knuth']
        """

        names = getattr(self, type + '_names')
        if abbr:
            import warnings
            warnings.warn('Person.get_part(abbr=True) is deprecated since 0.19: use pybtex.textutils.abbreviate()', stacklevel=2)
            from pybtex.textutils import abbreviate
            names = [abbreviate(name) for name in names]
        return names

    @property
    def rich_first_names(self):
        """
        A list of first names converted to :ref:`rich text <rich-text>`.

        .. versionadded:: 0.20
        """

        return [Text.from_latex(name) for name in self.first_names]

    @property
    def rich_middle_names(self):
        """
        A list of middle names converted to :ref:`rich text <rich-text>`.

        .. versionadded:: 0.20
        """
        return [Text.from_latex(name) for name in self.middle_names]

    @property
    def rich_prelast_names(self):
        """
        A list of pre-last (aka von) name parts converted to :ref:`rich text <rich-text>`.

        .. versionadded:: 0.20
        """
        return [Text.from_latex(name) for name in self.prelast_names]

    @property
    def rich_last_names(self):
        """
        A list of last names converted to :ref:`rich text <rich-text>`.

        .. versionadded:: 0.20
        """
        return [Text.from_latex(name) for name in self.last_names]

    @property
    def rich_lineage_names(self):
        """
        A list of lineage (aka Jr) name parts converted to :ref:`rich text <rich-text>`.

        .. versionadded:: 0.20
        """
        return [Text.from_latex(name) for name in self.lineage_names]

    @deprecated('0.19', 'use Person.first_names instead')
    def first(self, abbr=False):
        """
        .. deprecated:: 0.19
            Use :py:attr:`.first_names` instead.
        """
        return self.get_part('first', abbr)

    @deprecated('0.19', 'use Person.middle_names instead')
    def middle(self, abbr=False):
        """
        .. deprecated:: 0.19
            Use :py:attr:`.middle_names` instead.
        """
        return self.get_part('middle', abbr)

    @deprecated('0.19', 'use Person.prelast_names instead')
    def prelast(self, abbr=False):
        """
        .. deprecated:: 0.19
            Use :py:attr:`.prelast_names` instead.
        """
        return self.get_part('prelast', abbr)

    @deprecated('0.19', 'use Person.last_names instead')
    def last(self, abbr=False):
        """
        .. deprecated:: 0.19
            Use :py:attr:`.last_names` instead.
        """
        return self.get_part('last', abbr)

    @deprecated('0.19', 'use Person.lineage_names instead')
    def lineage(self, abbr=False):
        """
        .. deprecated:: 0.19
            Use :py:attr:`.lineage_names` instead.
        """
        return self.get_part('lineage', abbr)

    @deprecated('0.19', 'use Person.bibtex_first_names instead')
    def bibtex_first(self):
        """
        .. deprecated:: 0.19
            Use :py:attr:`.bibtex_first_names` instead.
        """
        return self.bibtex_first_names


def parse_file(file, bib_format=None, **kwargs):
    """
    Read bibliography data from file and return a :py:class:`.BibliographyData` object.

    :param file: A file name or a file-like object.
    :param bib_format: Data format ("bibtex", "yaml", etc.).
        If not specified, Pybtex will try to guess by the file name.

    .. versionadded:: 0.19
    """

    if isinstance(file, str):
        filename = file
    else:
        filename = getattr(file, 'name', None)

    parser = find_plugin('pybtex.database.input', bib_format, filename=filename)(**kwargs)
    return parser.parse_file(file)


def parse_string(value, bib_format, **kwargs):
    """
    Parse a Unicode string containing bibliography data and return a :py:class:`.BibliographyData` object.

    :param value: Unicode string.
    :param bib_format: Data format ("bibtex", "yaml", etc.).

    .. versionadded:: 0.19
    """

    parser = find_plugin('pybtex.database.input', bib_format)(**kwargs)
    return parser.parse_string(value)


def parse_bytes(value, bib_format, **kwargs):
    """
    Parse a byte string containing bibliography data and return a :py:class:`.BibliographyData` object.

    :param value: Byte string.
    :param bib_format: Data format (for example, "bibtexml").

    .. versionadded:: 0.19
    """

    parser = find_plugin('pybtex.database.input', bib_format)(**kwargs)
    return parser.parse_bytes(value)
