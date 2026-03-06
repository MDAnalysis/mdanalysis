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

"""convert bibliography database from one format to another
"""
from __future__ import unicode_literals
from pybtex.exceptions import PybtexError
from pybtex import database


class ConvertError(PybtexError):
    pass


def convert(
    from_filename, to_filename,
    from_format=None, to_format=None,
    input_encoding=None, output_encoding=None,
    parser_options=None,
    preserve_case=True,
    **kwargs
):
    if parser_options is None:
        parser_options = {}

    if from_filename == to_filename:
        raise ConvertError('input and output file can not be the same')

    bib_data = database.parse_file(
        from_filename,
        bib_format=from_format, encoding=input_encoding,
        **parser_options
    )
    if not preserve_case:
        bib_data = bib_data.lower()
    bib_data.to_file(to_filename, bib_format=to_format, encoding=output_encoding)
