# Adapted from mdtraj/devtools/travis-ci/update_versions_json.py by Robert T McGibbon
# by Lily Wang, 2020
# 
# This script is called by deploy_docs_via_travis.sh to update versions.json
# and re-write the following redirect stubs:
# index.html (-> latest release)
# latest.html (-> latest release)
# dev.html (-> dev docs)

import json
import os

try:
    from urllib.request import Request, urlopen
except ImportError:
    from urllib2 import Request, urlopen

# ========= WRITE JSON =========
URL = "https://mdanalysis.github.io/mdanalysis/"

VERSION = os.environ['VERSION']
url = os.path.join(URL, 'versions.json')

try:
    page = Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    data = urlopen(page).read().decode()
except Exception as e:
    print(e)
    try:
        with open('versions.json', 'r') as f:
            versions = json.loads(f)
    except IOError as e:
        print(e)
        versions = []
else:
    versions = json.loads(data)

existing = [item['version'] for item in versions]
already_exists = VERSION in existing

if not already_exists:
    latest = 'dev' not in VERSION
    if latest:
        for ver in versions:
            ver['latest'] = False

    versions.append({
        'version': VERSION,
        'display': VERSION,
        'url': os.path.join(URL, VERSION),
        'latest': latest
        })

with open("versions.json", 'w') as f:
    json.dump(versions, f, indent=2)

# ========= WRITE HTML STUBS =========
REDIRECT = """
<!DOCTYPE html>
<meta charset="utf-8">
<title>Redirecting to {url}</title>
<meta http-equiv="refresh" content="0; URL={url}">
<link rel="canonical" href="{url}">
"""

for ver in versions[::-1]:
    if ver['latest']:
        latest_url = ver['url']
        break
else:
    try:
        latest_url = versions[-1]['url']
    except IndexError:
        latest_url = URL

for ver in versions[::-1]:
    if 'dev' in ver['version']:
        dev_url = ver['url']
        break
else:
    try:
        dev_url = versions[-1]['url']
    except IndexError:
        dev_url = URL

with open('index.html', 'w') as f:
    f.write(REDIRECT.format(url=latest_url))

with open('latest.html', 'w') as f:
    f.write(REDIRECT.format(url=latest_url))

with open('dev.html', 'w') as f:
    f.write(REDIRECT.format(url=dev_url))