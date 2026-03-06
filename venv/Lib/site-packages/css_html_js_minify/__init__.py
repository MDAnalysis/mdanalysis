#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""CSS-HTML-JS-Minify.

Minifier for the Web.
"""


from .minify import (process_single_html_file, process_single_js_file,
                     process_single_css_file, html_minify, js_minify,
                     css_minify)


__version__ = '2.5.5'
__all__ = ('__version__', 'process_single_html_file', 'process_single_js_file',
           'process_single_css_file', 'html_minify', 'js_minify', 'css_minify',
           'minify')
