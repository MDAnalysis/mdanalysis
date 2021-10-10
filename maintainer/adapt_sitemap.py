#!/usr/bin/env python
#
#
# Adjust path in sitemap.xml

from xml.etree import ElementTree
import argparse

# defaults for MDAnalysis, see https://github.com/MDAnalysis/mdanalysis/pull/1890
# and https://github.com/MDAnalysis/MDAnalysis.github.io/issues/78

RELEASE_URL = "https://www.mdanalysis.org/docs/"
DEVELOP_URL = "https://www.mdanalysis.org/mdanalysis/"

# change if sitemaps.org updates their schema
NAMESPACE = {"sitemaps": "http://www.sitemaps.org/schemas/sitemap/0.9"}

def replace_loc(tree, search, replace, namespace=NAMESPACE):
    root = tree.getroot()
    urls = root.findall("sitemaps:url", namespace)
    if len(urls) == 0:
        raise ValueError("No sitemaps:url element found: check if the namespace in the XML file "
                         "is still xmlns='{0[sitemaps]}'".format(namespace))
    for url in urls:
        loc = url.find("sitemaps:loc", namespace)
        try:
            loc.text = loc.text.replace(search, replace)
        except AttributError:
            raise ValueError("No sitemaps:loc element found: check if the namespace in the XML file "
                             "is still xmlns='{0[sitemaps]}'".format(namespace))
    return tree


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Change top level loc in sitemap.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('sitemap', metavar="FILE",
                        help="path to sitemap.xml file, will be changed in place")
    parser.add_argument('--output', '-o', metavar="FILE",
                        default="sitemap_release.xml",
                        help="write altered XML to output FILE")
    parser.add_argument("--search", "-s", metavar="URL",
                        default=DEVELOP_URL,
                        help="search this URL in the loc elements")
    parser.add_argument("--replace", "-r", metavar="URL",
                        default=RELEASE_URL,
                        help="replace the searched URL with this URL in the loc elements")
    args = parser.parse_args()


    with open(args.sitemap) as xmlfile:
        tree = ElementTree.parse(xmlfile)

    tree = replace_loc(tree, args.search, args.replace)

    with open(args.output, "wb") as xmlfile:
        tree.write(xmlfile, encoding="utf-8", xml_declaration=True,
                   default_namespace=NAMESPACE['sitemaps'])

    print("adapt_sitemap.py: Created output file {} with change in loc:".format(args.output))
    print("adapt_sitemap.py: {0} --> {1}".format(args.search, args.replace))
