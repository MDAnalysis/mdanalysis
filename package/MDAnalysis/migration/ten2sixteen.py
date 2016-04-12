import sys
from lib2to3.main import main

if __name__ == '__main__':
    sys.exit(main("fixes"))

def ten2sixteen(*files):
    """Convert MDAnalysis 0.10.0 scripts to 0.16.0

    Usage:
       ten2sixteen('myscript.py')
    """
    main('MDAnalysis.migration.fixes', ['-w'] + list(files))
