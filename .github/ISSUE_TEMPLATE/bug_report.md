---
name: Bug report
about: Create a report to help us improve

---

## Expected behavior ##

<!-- A clear and concise description of what you want to do and what you think should happen. (Code to reproduce the behavior can be added below). -->


## Actual behavior ##

<!-- What happened instead. Add as much detail as you can. Include (copy and paste) stack traces and any output. -->


## Code to reproduce the behavior ##

<!-- Show us how to reproduce the failure. If you can, use trajectory files from the test data. Use the code snipped below as a starting point. -->

``` python
import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PSF, DCD,  GRO, PDB, TPR, XTC, TRR,  PRMncdf, NCDF

u = mda.Universe(PSF, DCD)

....

```

## Current version of MDAnalysis ##

- Which version are you using? (run `python -c "import MDAnalysis as mda; print(mda.__version__)"`)
- Which version of Python (`python -V`)?
- Which operating system?
