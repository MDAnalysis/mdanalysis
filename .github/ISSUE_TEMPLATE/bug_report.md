---
name: Bug report
about: Create a report to help us improve

---

**Expected behavior**

A clear and concise description of what you want to do and what you think should happen. (Code to reproduce the behavior can be added below).


**Actual behavior**

What happened instead. Add as much detail as you can. Include (copy and paste) stack traces and any output.


**Code to reproduce the behavior**

Show us how to reproduce the failiure. If you can, use trajectory files from the test data.

``` python
import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PSF, DCD,  GRO, PDB, TPR, XTC, TRR,  PRMncdf, NCDF

u = mda.Universe(PSF, DCD)

....

```

**Currently version of MDAnalysis**

- Which version are you using? (run `python -c "import MDAnalysis as mda; print(mda.__version__)"`)
- Which version of Python (`python -V`)?
- Which operating system?
