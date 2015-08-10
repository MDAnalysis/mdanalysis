#!/bin/bash
# -*- coding: utf-8 -*-
# Written by Sébastien Buchoux, 2012-2013
# Placed into the Public Domain

LOGFILE="$PWD/unittests_$(date +'%Y%m%d').txt"
NPROC=6
usage="$0 [options]

Build MDAnalysis and run the unit tests.

Options
-h                this help
-i                print information about python and paths
-o FILE           also write output to FILE,
                  set to /dev/null if you do not want it.
                  [default: ${LOGFILE}]
-n NPROC          run in parallel with NPROC processes [$NPROC]
"

function get_root_dir () {
    # Get root dir
    local cur_dir root_dir
    cur_dir="`pwd`"
    root_dir="`pwd`/`dirname $0`/.."
    root_dir="`cd \"$root_dir\"; pwd; cd \"$cur_dir\";`" # Trick to get a clean string
    echo "$root_dir"
}

function report_python_env () {
   # Report Env
    python -c "\
from __future__ import print_function;import sys;\
print('Python executable: %s' % sys.executable);
print('Python version: %s' % sys.version);
#print('Architecture: %s' % sys.arch);"
    echo ""
}


while getopts hio:n: OPT; do
    case "$OPT" in
	h)  echo "$usage"; exit 0;;
	i)  echo "rootdir: $(get_root_dir)";
            report_python_env;
	    exit 0;;        
        o)  LOGFILE=${OPTARG};;
	n)  NPROC=${OPTARG};;
	*)  echo "Unknown option $OPT"; exit 2;;
    esac
done




root_dir=$(get_root_dir)

# Build core
echo "Building MDAnalysis..."
cd "$root_dir/package"
rm -rf build
python setup.py build
cd "$cur_dir"

# Build testsuite
echo
echo "Building MDAnalysis Testsuite..."
cd "$root_dir/testsuite"
rm -rf build
python setup.py build
cd "$cur_dir"

# Get the lib dirs
lib_dir="`ls \"$root_dir/package/build\" | grep \"lib\"`"
package_lib="$root_dir/package/build/$lib_dir"

testsuite_lib_dir="`ls \"$root_dir/testsuite/build\" | grep \"lib\"`"
testsuite_lib="$root_dir/testsuite/build/$testsuite_lib_dir"

# Run the tests
NOSETESTS=$root_dir/testsuite/MDAnalysisTests/mda_nosetests 
echo
echo "Running the tests... (output to ${LOGFILE})"
echo "(Using MDAnalysis test script $NOSETESTS.)"
echo "Tests from $tessuite_lib"

cd "$package_lib"
$NOSETESTS -v --processes=$NPROC --process-timeout=120 \
	    --with-memleak "$testsuite_lib/MDAnalysisTests" 2>&1 | tee $LOGFILE
cd "$cur_dir"
