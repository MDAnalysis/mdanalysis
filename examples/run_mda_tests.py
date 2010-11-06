#!/usr/bin/env python
import MDAnalysis.tests
import sys

try:
  nproc = sys.argv[1]
except IndexError:
  nproc = 1

MDAnalysis.tests.test(label='full',verbose=3, 
  	extra_argv=['--exe', '--processes=%d' % int(nproc)])

