import sys
from lib2to3.main import main

if __name__ == '__main__':
    sys.exit(main("fixes"))

def ten2eleven(*files):
    """Convert MDAnalysis 0.10.0 scripts to 0.11.0

    Usage:
       ten2eleven('myscript.py')
    """
    main('MDAnalysis.migration.fixes', ['-w'] + list(files))
