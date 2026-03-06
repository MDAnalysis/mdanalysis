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

from __future__ import unicode_literals

import io
from xml.sax.saxutils import XMLGenerator
from xml.sax.xmlreader import AttributesImpl

from pybtex.database.output import BaseWriter


doctype = b"""<!DOCTYPE bibtex:file PUBLIC
    "-//BibTeXML//DTD XML for BibTeX v1.0//EN"
        "bibtexml.dtd" >
"""


class _PrettyXMLWriter(object):
    def __init__(self, output, encoding='UTF-8', namespace=('bibtex', 'http://bibtexml.sf.net/'), header=True):
        self.prefix, self.uri = namespace
        self.generator = XMLGenerator(output, encoding=encoding)
        if header:
            self.generator.startDocument()
        self.generator.startPrefixMapping(self.prefix, self.uri)
        self.stack = []

    def write(self, data):
        self.generator.characters(data)

    def newline(self):
        self.write('\n')

    def indent_line(self):
        self.write(' ' * (len(self.stack) * 4))

    def start(self, tag, attrs=None, newline=True):
        if attrs is None:
            attrs = {}
        else:
            attrs = {(None, key): value for key, value in attrs.items()}
        self.indent_line()
        self.stack.append(tag)
        self.generator.startElementNS((self.uri, tag), tag, AttributesImpl(attrs))
        if newline:
            self.newline()

    def end(self, indent=True):
        tag = self.stack.pop()
        if indent:
            self.indent_line()
        self.generator.endElementNS((self.uri, tag), tag)
        self.newline()

    def element(self, tag, data):
        self.start(tag, newline=False)
        self.write(data)
        self.end(indent=False)

    def close(self):
        self.generator.endDocument()


class Writer(BaseWriter):
    """Outputs BibTeXML markup"""

    def write_stream(self, bib_data, stream):
        xml_writer = _PrettyXMLWriter(stream, self.encoding)
        self._write(bib_data, xml_writer)

    def to_string(self, bib_data):
        """
        Return a unicode XML string without encoding declaration.

        >>> from pybtex.database import BibliographyData
        >>> data = BibliographyData()
        >>> unicode_xml = Writer().to_string(data)
        >>> isinstance(unicode_xml, str)
        True
        >>> print(unicode_xml)
        <bibtex:file xmlns:bibtex="http://bibtexml.sf.net/">
        <BLANKLINE>
        </bibtex:file>
        """

        output = io.BytesIO()
        writer = _PrettyXMLWriter(output, encoding='UTF-8', header=None)
        self._write(bib_data, writer)
        return output.getvalue().decode('UTF-8').strip()

    def _write(self, bib_data, writer):
        def write_persons(persons, role):
            if persons:
                writer.start(role)
                for person in persons:
                    writer.start('person')
                    for type in ('first', 'middle', 'prelast', 'last', 'lineage'):
                        name = person.get_part_as_text(type)
                        if name:
                            writer.element(type, name)
                    writer.end()
                writer.end()

        writer.start('file')
        writer.newline()
        for key, entry in bib_data.entries.items():
            writer.start('entry', dict(id=key))
            writer.start(entry.original_type)
            for field_name, field_value in entry.fields.items():
                writer.element(field_name, field_value)
            for role, persons in entry.persons.items():
                write_persons(persons, role)
            writer.end()
            writer.end()
            writer.newline()
        writer.end()
        writer.close()
