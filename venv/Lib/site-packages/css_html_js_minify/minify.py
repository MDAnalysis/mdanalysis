#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""CSS-HTML-JS-Minify.

Minifier for the Web.
"""


import atexit
import os
import sys
import gzip
import logging as log

from argparse import ArgumentParser
from datetime import datetime
from functools import partial
from hashlib import sha1
from multiprocessing import Pool, cpu_count
from subprocess import getoutput
from time import sleep

from .css_minifier import css_minify
from .html_minifier import html_minify
from .js_minifier import js_minify


__all__ = ('process_multiple_files', 'prefixer_extensioner',
           'process_single_css_file', 'process_single_html_file',
           'process_single_js_file', 'make_arguments_parser', 'main')


##############################################################################


def process_multiple_files(file_path, watch=False, wrap=False, timestamp=False,
                           comments=False, sort=False, overwrite=False,
                           zipy=False, prefix='', add_hash=False):
    """Process multiple CSS, JS, HTML files with multiprocessing."""
    print(f"Process {os.getpid()} is Compressing {file_path}.")
    if watch:
        previous = int(os.stat(file_path).st_mtime)
        print(f"Process {os.getpid()} is Watching {file_path}.")
        while True:
            actual = int(os.stat(file_path).st_mtime)
            if previous == actual:
                sleep(60)
            else:
                previous = actual
                print(f"Modification detected on {file_path}.")
                if file_path.endswith(".css"):
                    process_single_css_file(
                        file_path, wrap=wrap, timestamp=timestamp,
                        comments=comments, sort=sort, overwrite=overwrite,
                        zipy=zipy, prefix=prefix, add_hash=add_hash)
                elif file_path.endswith(".js"):
                    process_single_js_file(
                        file_path, timestamp=timestamp,
                        overwrite=overwrite, zipy=zipy)
                else:
                    process_single_html_file(
                        file_path, comments=comments,
                        overwrite=overwrite, prefix=prefix, add_hash=add_hash)
    else:
        if file_path.endswith(".css"):
            process_single_css_file(
                file_path, wrap=wrap, timestamp=timestamp,
                comments=comments, sort=sort, overwrite=overwrite, zipy=zipy,
                prefix=prefix, add_hash=add_hash)
        elif file_path.endswith(".js"):
            process_single_js_file(
                file_path, timestamp=timestamp,
                overwrite=overwrite, zipy=zipy)
        else:
            process_single_html_file(
                file_path, comments=comments,
                overwrite=overwrite, prefix=prefix, add_hash=add_hash)


def prefixer_extensioner(file_path, old, new,
                         file_content=None, prefix='', add_hash=False):
    """Take a file path and safely preppend a prefix and change extension.

    This is needed because filepath.replace('.foo', '.bar') sometimes may
    replace '/folder.foo/file.foo' into '/folder.bar/file.bar' wrong!.
    >>> prefixer_extensioner('/tmp/test.js', '.js', '.min.js')
    '/tmp/test.min.js'
    """
    print(f"Prepending '{new.upper()}' Prefix to {file_path}.")
    extension = os.path.splitext(file_path)[1].lower().replace(old, new)
    filenames = os.path.splitext(os.path.basename(file_path))[0]
    filenames = prefix + filenames if prefix else filenames
    if add_hash and file_content:  # http://stackoverflow.com/a/25568916
        filenames += "-" + sha1(file_content.encode("utf-8")).hexdigest()[:11]
        print(f"Appending SHA1 HEX-Digest Hash to '{file_path}'.")
    dir_names = os.path.dirname(file_path)
    file_path = os.path.join(dir_names, filenames + extension)
    return file_path


def process_single_css_file(css_file_path, wrap=False, timestamp=False,
                            comments=False, sort=False, overwrite=False,
                            zipy=False, prefix='', add_hash=False,
                            output_path=None):
    """Process a single CSS file."""
    print(f"Processing CSS file: {css_file_path}.")
    with open(css_file_path, encoding="utf-8") as css_file:
        original_css = css_file.read()
    print(f"INPUT: Reading CSS file {css_file_path}.")
    minified_css = css_minify(original_css, wrap=wrap,
                              comments=comments, sort=sort)
    if timestamp:
        taim = "/* {0} */ ".format(datetime.now().isoformat()[:-7].lower())
        minified_css = taim + minified_css
    if output_path is None:
        min_css_file_path = prefixer_extensioner(
            css_file_path, ".css", ".css" if overwrite else ".min.css",
            original_css, prefix=prefix, add_hash=add_hash)
        if zipy:
            gz_file_path = prefixer_extensioner(
                css_file_path, ".css",
                ".css.gz" if overwrite else ".min.css.gz", original_css,
                prefix=prefix, add_hash=add_hash)
            print(f"OUTPUT: Writing ZIP CSS {gz_file_path}.")
    else:
        min_css_file_path = gz_file_path = output_path
    if not zipy or output_path is None:
        # if specific output path is requested,write write only one output file
        with open(min_css_file_path, "w", encoding="utf-8") as output_file:
            output_file.write(minified_css)
    if zipy:
        with gzip.open(gz_file_path, "wt", encoding="utf-8") as output_gz:
            output_gz.write(minified_css)
    print(f"OUTPUT: Writing CSS Minified {min_css_file_path}.")
    return min_css_file_path


def process_single_html_file(html_file_path, comments=False, overwrite=False,
                             prefix='', add_hash=False, output_path=None):
    """Process a single HTML file."""
    print(f"Processing HTML file: {html_file_path}.")
    with open(html_file_path, encoding="utf-8") as html_file:
        minified_html = html_minify(html_file.read(), comments=comments)
    print(f"INPUT: Reading HTML file {html_file_path}.")
    if output_path is None:
        html_file_path = prefixer_extensioner(
            html_file_path, ".html" if overwrite else ".htm", ".html",
            prefix=prefix, add_hash=add_hash)
    else:
        html_file_path = output_path
    with open(html_file_path, "w", encoding="utf-8") as output_file:
        output_file.write(minified_html)
    print(f"OUTPUT: Writing HTML Minified {html_file_path}.")
    return html_file_path


def process_single_js_file(js_file_path, timestamp=False, overwrite=False,
                           zipy=False, output_path=None):
    """Process a single JS file."""
    print(f"Processing JS file: {js_file_path}.")
    with open(js_file_path, encoding="utf-8") as js_file:
        original_js = js_file.read()
    print(f"INPUT: Reading JS file {js_file_path}.")
    minified_js = js_minify(original_js)
    if timestamp:
        taim = "/* {} */ ".format(datetime.now().isoformat()[:-7].lower())
        minified_js = taim + minified_js
    if output_path is None:
        min_js_file_path = prefixer_extensioner(
            js_file_path, ".js", ".js" if overwrite else ".min.js",
            original_js)
        if zipy:
            gz_file_path = prefixer_extensioner(
                js_file_path, ".js", ".js.gz" if overwrite else ".min.js.gz",
                original_js)
            print(f"OUTPUT: Writing ZIP JS {gz_file_path}.")
    else:
        min_js_file_path = gz_file_path = output_path
    if not zipy or output_path is None:
        # if specific output path is requested,write write only one output file
        with open(min_js_file_path, "w", encoding="utf-8") as output_file:
            output_file.write(minified_js)
    if zipy:
        with gzip.open(gz_file_path, "wt", encoding="utf-8") as output_gz:
            output_gz.write(minified_js)
    print(f"OUTPUT: Writing JS Minified {min_js_file_path}.")
    return min_js_file_path


def make_arguments_parser():
    """Build and return a command line agument parser."""
    parser = ArgumentParser(description=__doc__, epilog="""CSS-HTML-JS-Minify:
    Takes a file or folder full path string and process all CSS/HTML/JS found.
    If argument is not file/folder will fail. Check Updates works on Python3.
    Std-In to Std-Out is deprecated since it may fail with unicode characters.
    SHA1 HEX-Digest 11 Chars Hash on Filenames is used for Server Cache.
    CSS Properties are Alpha-Sorted, to help spot cloned ones, Selectors not.
    Watch works for whole folders, with minimum of ~60 Secs between runs.""")
    # parser.add_argument('--version', action='version', version=css_html_js_minify.__version__)
    parser.add_argument('fullpath', metavar='fullpath', type=str,
                        help='Full path to local file or folder.')
    parser.add_argument('--wrap', action='store_true',
                        help="Wrap output to ~80 chars per line, CSS only.")
    parser.add_argument('--prefix', type=str,
                        help="Prefix string to prepend on output filenames.")
    parser.add_argument('--timestamp', action='store_true',
                        help="Add a Time Stamp on all CSS/JS output files.")
    parser.add_argument('--quiet', action='store_true', help="Quiet, Silent.")
    parser.add_argument('--hash', action='store_true',
                        help="Add SHA1 HEX-Digest 11chars Hash to Filenames.")
    parser.add_argument('--zipy', action='store_true',
                        help="GZIP Minified files as '*.gz', CSS/JS only.")
    parser.add_argument('--sort', action='store_true',
                        help="Alphabetically Sort CSS Properties, CSS only.")
    parser.add_argument('--comments', action='store_true',
                        help="Keep comments, CSS/HTML only (Not Recommended)")
    parser.add_argument('--overwrite', action='store_true',
                        help="Force overwrite all in-place (Not Recommended)")
    parser.add_argument('--after', type=str,
                        help="Command to execute after run (Experimental).")
    parser.add_argument('--before', type=str,
                        help="Command to execute before run (Experimental).")
    parser.add_argument('--watch', action='store_true', help="Watch changes.")
    parser.add_argument('--multiple', action='store_true',
                        help="Allow Multiple instances (Not Recommended).")
    return parser.parse_args()


def walk2list(folder: str, target: tuple, omit: tuple=(),
              showhidden: bool=False, topdown: bool=True,
              onerror: object=None, followlinks: bool=False) -> tuple:
    """Perform full walk, gather full path of all files."""
    oswalk = os.walk(folder, topdown=topdown,
                     onerror=onerror, followlinks=followlinks)

    return [os.path.abspath(os.path.join(r, f))
            for r, d, fs in oswalk
            for f in fs if not f.startswith(() if showhidden else ".") and
            not f.endswith(omit) and f.endswith(target)]


def main():
    """Main Loop."""
    args = make_arguments_parser()
    if os.path.isfile(args.fullpath) and args.fullpath.endswith(".css"):
        print("Target is a CSS File.")  # Work based on if argument is
        list_of_files = str(args.fullpath)  # file or folder, folder is slower.
        process_single_css_file(
            args.fullpath, wrap=args.wrap, timestamp=args.timestamp,
            comments=args.comments, sort=args.sort, overwrite=args.overwrite,
            zipy=args.zipy, prefix=args.prefix, add_hash=args.hash)
    elif os.path.isfile(args.fullpath) and args.fullpath.endswith(
            ".html" if args.overwrite else ".htm"):
        print("Target is HTML File.")
        list_of_files = str(args.fullpath)
        process_single_html_file(
            args.fullpath, comments=args.comments,
            overwrite=args.overwrite, prefix=args.prefix, add_hash=args.hash)
    elif os.path.isfile(args.fullpath) and args.fullpath.endswith(".js"):
        print("Target is a JS File.")
        list_of_files = str(args.fullpath)
        process_single_js_file(
            args.fullpath, timestamp=args.timestamp,
            overwrite=args.overwrite, zipy=args.zipy)
    elif os.path.isdir(args.fullpath):
        print("Target is a Folder with CSS, HTML, JS files !.")
        print("Processing a whole Folder may take some time...")
        list_of_files = walk2list(
            args.fullpath,
            (".css", ".js", ".html" if args.overwrite else ".htm"),
            (".min.css", ".min.js", ".htm" if args.overwrite else ".html"))
        print(f'Total Maximum CPUs used: ~{cpu_count()} Cores.')
        pool = Pool(cpu_count())  # Multiprocessing Async
        pool.map_async(partial(
                process_multiple_files, watch=args.watch,
                wrap=args.wrap, timestamp=args.timestamp,
                comments=args.comments, sort=args.sort,
                overwrite=args.overwrite, zipy=args.zipy,
                prefix=args.prefix, add_hash=args.hash),
            list_of_files)
        pool.close()
        pool.join()
    else:
        print("File or folder not found,or cant be read,or I/O Error.")
        sys.exit(1)
    if args.after and getoutput:
        print(getoutput(str(args.after)))
    print(f'\n {"-" * 80} \n Files Processed: {list_of_files}.')
    print(f'''Number of Files Processed:
          {len(list_of_files) if isinstance(list_of_files, tuple) else 1}.''')
