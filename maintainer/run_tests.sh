#!/bin/bash
# Written by Sébastien Buchoux, 2012-2013
# Placed into the Public Domain

# Get root dir
cur_dir="`pwd`"
root_dir="`pwd`/`dirname $0`/.."
root_dir="`cd \"$root_dir\"; pwd; cd \"$cur_dir\";`" # Trick to get a clean string

# Report Env
python -c "from __future__ import print_function;import sys;\
print('Python executable: %s' % sys.executable);
print('Python version: %s' % sys.version);
print('Architecture: %s' % sys.arch);"
echo ""



# Build core
echo "Building MDAnalysis..."
cd "$root_dir/package"
python setup.py build
cd "$cur_dir"

# Build testsuite
echo
echo "Building MDAnalysis Testsuite..."
cd "$root_dir/testsuite"
python setup.py build
cd "$cur_dir"

# Get the lib dirs
lib_dir="`ls "$root_dir/package/build" | grep \"lib\"`"
package_lib="$root_dir/package/build/$lib_dir"
testsuite_lib="$root_dir/testsuite/build/lib"

# Run the tests
echo
echo "Running the tests..."

cd "$package_lib"
nosetests -v "$testsuite_lib/MDAnalysisTests"
cd "$cur_dir"
