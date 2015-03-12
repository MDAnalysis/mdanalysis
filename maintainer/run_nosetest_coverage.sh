#!/bin/bash
# Coverage test of MDAnalysis testsuite with the coverage plugin of nose.
# author: Oliver Beckstein 2014
# licence: Public Domain

# Install MDAnalysis and MDAnalysisTests ('python setup.py develop' is good enough).

QUALITYDIR="$(cd $(dirname $0) && pwd)"

TESTSUITE="$QUALITYDIR/../testsuite"
HTMLDIR="$QUALITYDIR/htmlcov"
LOGFILE="$QUALITYDIR/testing.log"

echo "-- testing everything under $TESTSUITE"
echo "-- HTML coverage report in $HTMLDIR"
echo "-- all output in $LOGFILE"

# clean everything because otherwise coverage does not seem to update all files
rm -rf "$HTMLDIR"

# For some reason I had to run nosetests on separate file to make it
# generate the html reports; it didn't work if I just said
# "MDAnalysisTests" or "."; running 'coverage html' on the existing
# .coverage file works but uses long pathnames for filenames, which is
# ugly.

(cd "$TESTSUITE";
rm -f .coverage .noseids testing.log
nosetests-2.7 -v --with-id \
   --with-coverage --cover-erase --cover-html-dir="$HTMLDIR" --cover-html --cover-package=MDAnalysis \
   MDAnalysisTests/test_*.py  \
   2>&1 | tee "$LOGFILE";
)


