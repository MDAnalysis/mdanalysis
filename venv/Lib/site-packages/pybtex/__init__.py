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

"""bibliography processor
"""

from __future__ import absolute_import
from __future__ import unicode_literals
from os import path


__version__ = "0.25.1"


class Engine(object):
    def make_bibliography(self, aux_filename, style=None, output_encoding=None, bib_format=None, **kwargs):
        """
        Read the given ``.aux`` file and produce a formatted bibliography
        using :py:meth:`~.Engine.format_from_files`.

        :param style: If not ``None``, use this style instead of specified in the ``.aux`` file.
        """

        from pybtex import auxfile
        if bib_format is None:
            from pybtex.database.input.bibtex import Parser as bib_format

        aux_data = auxfile.parse_file(aux_filename, output_encoding)
        if style is None:
            style = aux_data.style
        base_filename = path.splitext(aux_filename)[0]
        bib_filenames = [filename + bib_format.default_suffix for filename in aux_data.data]
        return self.format_from_files(
            bib_filenames,
            style=aux_data.style,
            citations=aux_data.citations,
            output_encoding=output_encoding,
            output_filename=base_filename,
            add_output_suffix=True,
            **kwargs
        )

    def format_from_string(self, bib_string, *args, **kwargs):
        """
        Parse the bigliography data from the given string and produce a formated
        bibliography using :py:meth:`~.Engine.format_from_files`.

        This is a convenience method that calls
        :py:meth:`~.Engine.format_from_strings` with a single string.
        """
        return self.format_from_strings([bib_string], *args, **kwargs)

    def format_from_strings(self, bib_strings, *args, **kwargs):
        """
        Parse the bigliography data from the given strings and produce a formated
        bibliography.

        This is a convenience method that wraps each string into a StringIO,
        then calls :py:meth:`~.Engine.format_from_files`.
        """
        from io import StringIO
        inputs = [StringIO(bib_string) for bib_string in bib_strings]
        return self.format_from_files(inputs, *args, **kwargs)

    def format_from_file(self, filename, *args, **kwargs):
        """
        Read the bigliography data from the given file and produce a formated
        bibliography.

        This is a convenience method that calls :py:meth:`~.Engine.format_from_files`
        with a single file. All extra arguments are passed to
        :py:meth:`~.Engine.format_from_files`.
        """
        return self.format_from_files([filename], *args, **kwargs)

    def format_from_files(*args, **kwargs):
        """
        Read the bigliography data from the given files and produce a formated
        bibliography.

        This is an abstract method overridden by both
        :py:class:`pybtex.PybtexEngine` and :py:class:`pybtex.bibtex.BibTeXEngine`.
        """
        raise NotImplementedError


class PybtexEngine(Engine):
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
        output_backend=None,
        output_encoding=None,
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
        :param output_backend: Which output backend to use. The default is ``latex``.
        :param output_encoding: Encoding that will be used by the output backend.
        :param bst_encoding: Encoding of the ``.bst`` file.
        :param min_crossrefs: Include cross-referenced entries after this many
            crossrefs. See BibTeX manual for details.
        :param output_filename: If ``None``, the result will be returned as a
            string. Else, the result will be written to the specified file.
        :param add_output_suffix: Append default suffix to the output file
            name (``.bbl`` for LaTeX, ``.html`` for HTML, etc.).
        """

        from pybtex.plugin import find_plugin

        bib_parser = find_plugin('pybtex.database.input', bib_format)
        bib_data = bib_parser(
            encoding=bib_encoding,
            wanted_entries=citations,
            min_crossrefs=min_crossrefs,
        ).parse_files(bib_files_or_filenames)

        style_cls = find_plugin('pybtex.style.formatting', style)
        style = style_cls(
            label_style=kwargs.get('label_style'),
            name_style=kwargs.get('name_style'),
            sorting_style=kwargs.get('sorting_style'),
            abbreviate_names=kwargs.get('abbreviate_names'),
            min_crossrefs=min_crossrefs,
        )
        formatted_bibliography = style.format_bibliography(bib_data, citations)

        output_backend = find_plugin('pybtex.backends', output_backend)
        if add_output_suffix:
            output_filename = output_filename + output_backend.default_suffix
        if not output_filename:
            import io
            output_filename = io.StringIO()
        return output_backend(output_encoding).write_to_file(formatted_bibliography, output_filename)


def make_bibliography(*args, **kwargs):
    """A convenience function that calls :py:meth:`.PybtexEngine.make_bibliography`."""
    return PybtexEngine().make_bibliography(*args, **kwargs)


def format_from_file(*args, **kwargs):
    """A convenience function that calls :py:meth:`.PybtexEngine.format_from_file`."""
    return PybtexEngine().format_from_file(*args, **kwargs)


def format_from_files(*args, **kwargs):
    """A convenience function that calls :py:meth:`.PybtexEngine.format_from_files`."""
    return PybtexEngine().format_from_files(*args, **kwargs)


def format_from_string(*args, **kwargs):
    """A convenience function that calls :py:meth:`.PybtexEngine.format_from_string`."""
    return PybtexEngine().format_from_string(*args, **kwargs)


def format_from_strings(*args, **kwargs):
    """A convenience function that calls :py:meth:`.PybtexEngine.format_from_strings`."""
    return PybtexEngine().format_from_strings(*args, **kwargs)
