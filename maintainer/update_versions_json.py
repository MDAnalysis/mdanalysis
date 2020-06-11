import json
import os

try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen

URL = "https://www.mdanalysis.org/mdanalysis/"

VERSION = os.environ['VERSION']
url = os.path.join(URL, 'versions.json')
try:
    data = urlopen(url).read().decode()
    versions = json.loads(data)
except:
    try:
        with open('versions.json', 'r') as f:
            versions = json.loads(f)
    except:
        versions = []

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

