#!/bin/bash
# Written by Oliver Beckstein, 2012
# Placed into the Public Domain

usage="usage: $(basename $0) RELEASE

Change the RELEASE in all relevant files in MDAnalysis. RELEASE should be something
like '0.7.6' or 0.7.6-devel'.

Run from the top directory of the git checkout.
"

FILES="package/setup.py package/MDAnalysis/__init__.py testsuite/setup.py testsuite/MDAnalysisTests/__init__.py"

die () {
    echo "ERROR: $1"
    exit ${2:-1}
}

findcommand() {
    for name in $*; do
	path=$(which $name)
	if [ -n "$path" ]; then
	    echo $path
	    return 0
	fi
    done
    die "None of the commands $* found." 2
}


while getopts h OPT; do
    case $OPT in
	h) echo "$usage";
	   exit 0
	   ;;
	\?) exit 2;;
    esac
done

shift $((OPTIND - 1))
RELEASE=$1

test -n "$RELEASE" || die "Required argument missing. See -h for help." 2

# find a sed with -i
SED=$(findcommand gsed sed)

echo "Setting RELEASE/__version__ in MDAnalysis to $RELEASE"

git grep -E -l 'RELEASE.*[0-9]+\.[0-9]+\.[0-9]+(-devel)?' $FILES  \
   | xargs -I FILE $SED -i~  '/RELEASE/s/[0-9]\+\.[0-9]\+\.[0-9]\+\(-devel\)\?/'${RELEASE}'/' FILE
git grep -E -l '__version__ =.*[0-9]+\.[0-9]+\.[0-9]+(-devel)?'  $FILES \
   | xargs -I FILE $SED -i~ '/__version__/s/[0-9]\+\.[0-9]\+\.[0-9]\+\(-devel\)\?/'${RELEASE}'/' FILE

