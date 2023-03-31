# Copyright (c) 2013 Michael Dowling <mtdowling@gmail.com>
# Copyright (c) 2017 Jared Dillard <jared.dillard@gmail.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

import os
import queue
import xml.etree.ElementTree as ET
from multiprocessing import Manager

from sphinx.util.logging import getLogger

__version__ = "2.5.0"

logger = getLogger(__name__)


def setup(app):
    """Setup connects events to the sitemap builder"""
    app.add_config_value("site_url", default=None, rebuild="")
    app.add_config_value(
        "sitemap_url_scheme", default="{lang}{version}{link}", rebuild=""
    )
    app.add_config_value("sitemap_locales", default=None, rebuild="")

    app.add_config_value("sitemap_filename", default="sitemap.xml", rebuild="")

    try:
        app.add_config_value("html_baseurl", default=None, rebuild="")
    except BaseException:
        pass

    app.connect("builder-inited", record_builder_type)
    app.connect("html-page-context", add_html_link)
    app.connect("build-finished", create_sitemap)

    return {
        "parallel_read_safe": True,
        "parallel_write_safe": True,
        "version": __version__,
    }


def get_locales(app, exception):
    # Manually configured list of locales
    sitemap_locales = app.builder.config.sitemap_locales
    if sitemap_locales:
        # special value to add nothing -> use primary language only
        if sitemap_locales == [None]:
            return []

        # otherwise, add each locale
        return [locale for locale in sitemap_locales]

    # Or autodetect locales
    locales = []
    for locale_dir in app.builder.config.locale_dirs:
        locale_dir = os.path.join(app.confdir, locale_dir)
        if os.path.isdir(locale_dir):
            for locale in os.listdir(locale_dir):
                if os.path.isdir(os.path.join(locale_dir, locale)):
                    locales.append(locale)
    return locales


def record_builder_type(app):
    # builder isn't initialized in the setup so we do it here
    builder = getattr(app, "builder", None)
    if builder is None:
        return
    builder.env.is_directory_builder = type(builder).__name__ == "DirectoryHTMLBuilder"
    builder.env.sitemap_links = Manager().Queue()


def hreflang_formatter(lang):
    """
    sitemap hreflang should follow correct format.
        Use hyphen instead of underscore in language and country value.
    ref: https://en.wikipedia.org/wiki/Hreflang#Common_Mistakes
    source: https://github.com/readthedocs/readthedocs.org/pull/5638
    """
    if "_" in lang:
        return lang.replace("_", "-")
    return lang


def add_html_link(app, pagename, templatename, context, doctree):
    """As each page is built, collect page names for the sitemap"""
    env = app.builder.env
    if app.builder.config.html_file_suffix is None:
        file_suffix = ".html"
    else:
        file_suffix = app.builder.config.html_file_suffix

    # Support DirectoryHTMLBuilder path structure
    # where generated links between pages omit the index.html
    if env.is_directory_builder:
        if pagename == "index":
            sitemap_link = ""
        elif pagename.endswith("/index"):
            sitemap_link = pagename[:-6] + "/"
        else:
            sitemap_link = pagename + "/"
    else:
        sitemap_link = pagename + file_suffix

    env.sitemap_links.put(sitemap_link)


def create_sitemap(app, exception):
    """Generates the sitemap.xml from the collected HTML page links"""
    site_url = app.builder.config.site_url or app.builder.config.html_baseurl
    if site_url:
        site_url.rstrip("/") + "/"
    else:
        logger.warning(
            "sphinx-sitemap: html_baseurl is required in conf.py." "Sitemap not built.",
            type="sitemap",
            subtype="configuration",
        )
        return

    env = app.builder.env
    if env.sitemap_links.empty():
        logger.info(
            "sphinx-sitemap: No pages generated for %s" % app.config.sitemap_filename,
            type="sitemap",
            subtype="information",
        )
        return

    ET.register_namespace("xhtml", "http://www.w3.org/1999/xhtml")

    root = ET.Element("urlset", xmlns="http://www.sitemaps.org/schemas/sitemap/0.9")

    locales = get_locales(app, exception)

    if app.builder.config.version:
        version = app.builder.config.version + "/"
    else:
        version = ""

    while True:
        try:
            link = env.sitemap_links.get_nowait()
        except queue.Empty:
            break

        url = ET.SubElement(root, "url")
        scheme = app.config.sitemap_url_scheme
        if app.builder.config.language:
            lang = app.builder.config.language + "/"
        else:
            lang = ""

        ET.SubElement(url, "loc").text = site_url + scheme.format(
            lang=lang, version=version, link=link
        )

        for lang in locales:
            lang = lang + "/"
            ET.SubElement(
                url,
                "{http://www.w3.org/1999/xhtml}link",
                rel="alternate",
                hreflang=hreflang_formatter(lang.rstrip("/")),
                href=site_url + scheme.format(lang=lang, version=version, link=link),
            )

    filename = app.outdir + "/" + app.config.sitemap_filename
    ET.ElementTree(root).write(
        filename, xml_declaration=True, encoding="utf-8", method="xml"
    )

    logger.info(
        "sphinx-sitemap: %s was generated for URL %s in %s"
        % (app.config.sitemap_filename, site_url, filename),
        type="sitemap",
        subtype="information",
    )
