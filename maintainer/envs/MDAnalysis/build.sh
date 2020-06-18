#!/bin/bash

MDA_USE_OPENMP=FALSE pip install package/ --no-deps
pip install testsuite/  --no-deps
