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


"""
LaTeX output backend.

>>> from pybtex.richtext import Tag, HRef
>>> latex = Backend()
>>> print(Tag('em', '').render(latex))
<BLANKLINE>
>>> print(Tag('em', 'Non-', 'empty').render(latex))
\\emph{Non-empty}
>>> print(HRef('/', '').render(latex))
<BLANKLINE>
>>> print(HRef('/', 'Non-', 'empty').render(latex))
\\href{/}{Non-empty}
>>> print(HRef('http://example.org/', 'http://example.org/').render(latex))
\\url{http://example.org/}
>>> print(Tag('sup', 'hello').render(latex))
\\textsuperscript{hello}
"""
from __future__ import unicode_literals

import codecs

import latexcodec  # noqa
from pybtex.backends import BaseBackend


class Backend(BaseBackend):
    default_suffix = '.bbl'
    symbols = {
        'ndash': u'--',
        'newblock': u'\n\\newblock ',
        'nbsp': u'~'
    }
    tags = {
        u'em': u'emph',
        u'strong': None,
        u'i': u'textit',
        u'b': u'textbf',
        u'tt': u'texttt',
        u'sup': u'textsuperscript',
        u'sub': u'textsubscript',
    }

    def __init__(self, encoding=None):
        super(Backend, self).__init__(encoding)
        self.latex_encoding = 'ulatex+' + self.encoding

    def format_str(self, str_):
        return codecs.encode(str_, self.latex_encoding)

    def format_tag(self, tag_name, text):
        tag = self.tags.get(tag_name)
        if tag is None:
            return u'{%s}' % text if text else u''
        else:
            return r'\%s{%s}' % (tag, text) if text else u''

    def format_href(self, url, text, external=False):
        if not text:
            return ''
        opts = '[pdfnewwindow]' if external else ''
        if text == self.format_str(url):
            return u'\\url{%s}' % url
        else:
            return r'\href%s{%s}{%s}' % (opts, url, text)

    def format_protected(self, text):
        """
        >>> from pybtex.richtext import Protected
        >>> print(Protected('CTAN').render_as('latex'))
        {CTAN}
        """

        return '{%s}' % text

    def write_prologue(self):
        if self.formatted_bibliography.preamble:
            self.output(self.formatted_bibliography.preamble + u'\n')

        longest_label = self.formatted_bibliography.get_longest_label()
        self.output(u'\\begin{thebibliography}{%s}' % longest_label)

    def write_epilogue(self):
        self.output(u'\n\n\\end{thebibliography}\n')

    def write_entry(self, key, label, text):
        self.output(u'\n\n\\bibitem[%s]{%s}\n' % (label, key))
        self.output(text)
