# Copyright (c) 2006-2021  Andrey Golovizin
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

"""BibTeX unnamed stack language interpreter and related stuff
"""
from __future__ import unicode_literals


from os import path


from pybtex import Engine


class BibTeXEngine(Engine):
    """
    The Python fomatting engine.

    See :py:class:`pybtex.Engine` for inherited methods.
    """

    def format_from_files(
        self,
        bib_files_or_filenames,
        style,
        citations=['*'],
        bib_format=None,
        bib_encoding=None,
        output_encoding=None,
        bst_encoding=None,
        min_crossrefs=2,
        output_filename=None,
        add_output_suffix=False,
        **kwargs
    ):
        """
        Read the bigliography data from the given files and produce a formated
        bibliography.

        :param bib_files_or_filenames: A list of file names or file objects.
        :param style: The name of the formatting style.
        :param citations: A list of citation keys.
        :param bib_format: The name of the bibliography format. The default
            format is ``bibtex``.
        :param bib_encoding: Encoding of bibliography files.
        :param output_encoding: Encoding that will be used by the output backend.
        :param bst_encoding: Encoding of the ``.bst`` file.
        :param min_crossrefs: Include cross-referenced entries after this many
            crossrefs. See BibTeX manual for details.
        :param output_filename: If ``None``, the result will be returned as a
            string. Else, the result will be written to the specified file.
        :param add_output_suffix: Append a ``.bbl`` suffix to the output file name.
        """

        from io import StringIO
        import pybtex.io
        from pybtex.bibtex import bst
        from pybtex.bibtex.interpreter import Interpreter

        if bib_format is None:
            from pybtex.database.input.bibtex import Parser as bib_format
        bst_filename = style + path.extsep + 'bst'
        bst_script = bst.parse_file(bst_filename, bst_encoding)
        interpreter = Interpreter(bib_format, bib_encoding)
        bbl_data = interpreter.run(bst_script, citations, bib_files_or_filenames, min_crossrefs=min_crossrefs)

        if add_output_suffix:
            output_filename = output_filename + '.bbl'
        if output_filename:
            output_file = pybtex.io.open_unicode(output_filename, 'w', encoding=output_encoding)
        else:
            output_file = StringIO()
        with output_file:
            output_file.write(bbl_data)
            if isinstance(output_file, StringIO):
                return output_file.getvalue()


def make_bibliography(*args, **kwargs):
    """A convenience function that calls :py:meth:`.BibTeXEngine.make_bibliography`."""
    return BibTeXEngine().make_bibliography(*args, **kwargs)


def format_from_file(*args, **kwargs):
    """A convenience function that calls :py:meth:`.BibTeXEngine.format_from_file`."""
    return BibTeXEngine().format_from_file(*args, **kwargs)


def format_from_files(*args, **kwargs):
    """A convenience function that calls :py:meth:`.BibTeXEngine.format_from_files`."""
    return BibTeXEngine().format_from_files(*args, **kwargs)


def format_from_string(*args, **kwargs):
    """A convenience function that calls :py:meth:`.BibTeXEngine.format_from_string`."""
    return BibTeXEngine().format_from_string(*args, **kwargs)


def format_from_strings(*args, **kwargs):
    """A convenience function that calls :py:meth:`.BibTeXEngine.format_from_strings`."""
    return BibTeXEngine().format_from_strings(*args, **kwargs)
