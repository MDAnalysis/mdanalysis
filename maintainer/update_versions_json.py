import json
import MDAnalysis as mda

try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen

URL = "mdadocs.minium.com.au"

data = urlopen(URL + '/versions.json').read().decode()
versions = json.loads(data)

existing = [item['version'] for item in versions]
already_exists = mda.__version__ in existing

if not already_exists:
    for ver in versions:
        ver['latest'] = False

    versions.append({
        'version': mda.__version__,
        'display': mda.__version__,
        'url': os.path.join(URL, mda.__version),
        'latest': True
        })

print(json.dumps(versions, indent=2))
