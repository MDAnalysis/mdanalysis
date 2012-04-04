#!/bin/bash
# Oliver Beckstein 2012
# Placed into the Public Domain

: ${NETRC:=$HOME/.netrc}

usage="$(basename $0) [options] RELEASE

Upload the MDAnalysis package and tests MDAnalysisTests to google code.
Requires the googlecode_upload.py script.

Options

-h       help
-n       dry run (only check if everything needed is present)

Username and password are read from NETRC=$NETRC.
"

DRYRUN=""
while getopts nh OPT; do
    case $OPT in
	h) echo "$usage";
           exit 0;;
        n) echo "Ok, only doing a dry-run (no upload)";
           DRYRUN="echo";;
       \?) exit 1;;
    esac
done
shift $((OPTIND - 1))

RELEASE=$1
test -n "$RELEASE" || { echo "Need a release, see -h for usage"; exit 1; }

PACKAGE=./package/dist/MDAnalysis-${RELEASE}.tar.gz 
TESTSUITE=./testsuite/dist/MDAnalysisTests-${RELEASE}.tar.gz 

for tarball in $PACKAGE $TESTSUITE; do
    test -e "$tarball" || { echo "missing tarball ${tarball}: generate with 'python sdist'"; exit 2; }
done

test -e "$NETRC"  || { echo "need google code username and password in $NETRC"; exit 2; }

GUSER=$(awk '/^machine *code.google.com/ {sub("@gmail\.com", "", $4); printf $4}' $NETRC)
GPASSWORD=$(awk '/^machine *code.google.com/ {printf $6}' $NETRC)

UPLOAD="$DRYRUN googlecode_upload.py -p mdanalysis -u $GUSER -w $GPASSWORD -l Featured,Type-Source,OpSys-Linux,OpSys-OSX"

$UPLOAD -s "release $RELEASE (fixes and enhancements)"  $PACKAGE
$UPLOAD -s "tests and test data $RELEASE (for MDAnalysis ${RELEASE})" $TESTSUITE
