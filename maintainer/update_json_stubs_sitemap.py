# Adapted from mdtraj/devtools/travis-ci/update_versions_json.py,
# by Robert T McGibbon, under the LGPLv2.1 license
# Similarly by the terms of LGPL this is also in the public domain
# Lily Wang, 2020
#
# This script is called by deploy_docs_via_travis.sh to:
#  1. update versions.json
#  2. Write some redirect stubs
#  3. Write a sitemap.xml file for the root directory
#

import json
import os
import shutil
import xml.etree.ElementTree as ET
import errno
import glob
import textwrap
import shutil

try:
    from urllib.request import Request, urlopen
except ImportError:
    from urllib2 import Request, urlopen

URL = os.environ['URL']
VERSION = os.environ['VERSION']

if "http" not in URL:
    raise ValueError("URL should have the transfer protocol (HTTP/S). "
                     f"Given: $URL={URL}")

try:
    int(VERSION[0])
except ValueError:
    raise ValueError("$VERSION should start with a number. "
                     f"Given: $VERSION={VERSION}") from None


def get_web_file(filename, callback, default):
    url = os.path.join(URL, filename)
    try:
        page = Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        data = urlopen(page).read().decode()
    except Exception as e:
        print(e)
        try:
            with open(filename, 'r') as f:
                return callback(f)
        except IOError as e:
            print(e)
            return default
    else:
        return callback(data)


def write_redirect(file, version='', outfile=None):
    if outfile is None:
        outfile = file
    url = os.path.join(URL, version, file)
    REDIRECT = textwrap.dedent(f"""
    <!DOCTYPE html>
    <meta charset="utf-8">
    <title>Redirecting to {url}</title>
    <meta http-equiv="refresh" content="0; URL={url}">
    <link rel="canonical" href="{url}">
    """)
    with open(outfile, 'w') as f:
        f.write(REDIRECT)
    print(f"Wrote redirect from {url} to {outfile}")


# ========= WRITE JSON =========
# Update $root/versions.json with links to the right version
versions = get_web_file('versions.json', json.loads, [])
existing = [item['version'] for item in versions]
already_exists = VERSION in existing
latest = 'dev' not in VERSION

if not already_exists:
    if latest:
        for ver in versions:
            ver['latest'] = False

    versions.append({
        'version': VERSION,
        'display': VERSION,
        'url': os.path.join(URL, VERSION),
        'latest': latest
    })

for ver in versions[::-1]:
    if ver['latest']:
        latest_version = ver['version']
        break
else:
    try:
        latest_version = versions[-1]['version']
    except IndexError:
        latest_version = None

for ver in versions[::-1]:
    if 'dev' in ver['version']:
        dev_version = ver['version']
        break
else:
    try:
        dev_version = versions[-1]['version']
    except IndexError:
        dev_version = None

versions.sort(key=lambda x: x["version"])

# ========= WRITE HTML STUBS =========
# Add HTML files to redirect:
# index.html -> latest release
# latest/index.html -> latest release
# dev/index.html -> dev docs

if latest:
    shutil.copytree(VERSION, "stable")
    print(f"Copied {VERSION} to stable/")
    html_files = glob.glob(f'stable/**/*.html', recursive=True)
    for file in html_files:
        outfile = file[7:]  # strip "stable/"
        dirname = os.path.dirname(outfile)
        if dirname and not os.path.exists(dirname):
            try:
                os.makedirs(dirname)
            except OSError as exc:
                if exc.errno != errno.EEXIST:
                    raise

        write_redirect(file, '', outfile)



if latest_version:
    write_redirect('index.html', "stable")
    write_redirect('index.html', latest_version, 'latest/index.html')
    for ver in versions:
        if ver["version"] == "stable":
            ver["url"] = os.path.join(URL, "stable")
            break
    else:
        versions.append({
            "version": "stable",
            "display": "stable",
            "url": os.path.join(URL, "stable"),
            "latest": False
        })


if dev_version:
    if dev_version == VERSION:
        shutil.copytree(VERSION, "dev")
        print(f"Copied {VERSION} to dev/")

    for ver in versions:
        if ver["version"] == "dev":
            ver["url"] = os.path.join(URL, "dev")
            break
    else:
        versions.append({
            "version": "dev",
            "display": "dev",
            "url": os.path.join(URL, "dev"),
            "latest": False
        })


with open("versions.json", 'w') as f:
    json.dump(versions, f, indent=2)

# ========= WRITE SUPER SITEMAP.XML =========
# make one big sitemap.xml
ET.register_namespace('xhtml', "http://www.w3.org/1999/xhtml")

# so we could make 1 big sitemap as commented
# below, but they must be max 50 MB / 50k URL.
# Yes, this is 100+ releases, but who knows when
# that'll happen and who'll look at this then?
# bigroot = ET.Element("urlset")
# bigroot.set("xmlns", "http://www.sitemaps.org/schemas/sitemap/0.9")
# for ver in versions:
#     tree = get_web_file(ver['version']+'/sitemap.xml', ET.fromstring,
#                         ET.fromstring(''))
#     root = tree.getroot()
#     bigroot.extend(root.getchildren())
# ET.ElementTree(bigroot).write('sitemap.xml',
#                               xml_declaration=True,
#                               encoding='utf-8',
#                               method="xml")

# so instead we make a sitemap of sitemaps.
bigroot = ET.Element("sitemapindex")
bigroot.set("xmlns", "http://www.sitemaps.org/schemas/sitemap/0.9")
for ver in versions:
    path = os.path.join(URL, '{}/sitemap.xml'.format(ver['version']))
    sitemap = ET.SubElement(bigroot, 'sitemap')
    ET.SubElement(sitemap, 'loc').text = path

ET.ElementTree(bigroot).write('sitemap_index.xml',
                              xml_declaration=True,
                              encoding='utf-8',
                              method="xml")
