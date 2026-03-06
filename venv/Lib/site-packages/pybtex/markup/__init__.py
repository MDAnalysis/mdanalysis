from __future__ import unicode_literals
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


from pybtex.scanner import Scanner, Literal, PybtexSyntaxError
from pybtex.richtext import String, Text, Protected


class LaTeXParser(Scanner):
    LBRACE = Literal(u'{')
    RBRACE = Literal(u'}')

    def parse(self, level=0):
        """
        >>> LaTeXParser('abc').parse()
        Text('abc')

        >>> LaTeXParser('abc{def}').parse()
        Text('abc', Protected('def'))

        >>> LaTeXParser('abc{def {xyz}} !').parse()
        Text('abc', Protected('def ', Protected('xyz')), ' !')
        """

        return Text(*self.iter_string_parts(level=level))

    def iter_string_parts(self, level=0):
        while True:
            token = self.skip_to([self.LBRACE, self.RBRACE])
            if not token:
                remainder = self.get_remainder()
                if remainder:
                    yield String(remainder)
                if level != 0:
                    raise PybtexSyntaxError('unbalanced braces', self)
                break
            elif token.pattern is self.LBRACE:
                yield String(token.value[:-1])
                yield Protected(*self.iter_string_parts(level=level + 1))
            else:  # brace.pattern is self.RBRACE
                yield String(token.value[:-1])
                if level == 0:
                    raise PybtexSyntaxError('unbalanced braces', self)
                break
