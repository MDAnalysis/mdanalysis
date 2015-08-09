#!/bin/bash
# Written by Oliver Beckstein, 2012
# Placed into the Public Domain

usage="usage: $(basename $0) RELEASE

Change the RELEASE in all relevant files in MDAnalysis. RELEASE should be something
like '0.7.6' or 0.7.6-dev'.

Run from the top directory of the git checkout.
"

FILES="package/setup.py package/MDAnalysis/version.py testsuite/setup.py testsuite/MDAnalysisTests/__init__.py"

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

sed_eregex () {
  # does this sed understand extended regular expressions?
  # - FreeBSD (Mac OS X) sed -E
  # - GNU sed --regexp-extended (undocumented: also -E ...)
  local SED=$1
  if [ "good" = "$(echo 'bad' | $SED -E 's/(abc)?bad|foo/good/')" ]; then
     echo "$SED -E"
     return 0
  elif [ "good" = "$(echo 'bad' | $SED --regexp-extended 's/(abc)?bad|foo/good/')" ]; then     
     echo "$SED --regexp-extended"
     return 0
  elif [ "good" = "$(echo 'bad' | $SED 's/(abc)?bad|foo/good/')" ]; then
     echo "$SED"
     return 0
  fi
  echo "false"
  return 1
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

# find a sed with -i and -E
for cmd in gsed sed; do 
   SED=$(sed_eregex $(findcommand $cmd))
   [ "$SED" != "false" ] && break
done
[ "$SED" = "false" ] && { echo "ERROR: cannot find suitable sed."; exit 1; } 
# should check for -i but we just hope for the best ...
# modern(ish) seds have -i
echo "Using sed = $SED"

echo "Setting RELEASE/__version__ in MDAnalysis to $RELEASE"

git grep -E -l 'RELEASE.*[0-9]+\.[0-9]+\.[0-9]+(\.[0-9]+)?(dev|rc*[0-9])?' $FILES  \
   | xargs -I FILE $SED '/RELEASE/s/[0-9]+\.[0-9]+\.[0-9]+(\.[0-9]+)?(dev|rc*[0-9])?/'${RELEASE}'/' -i.bak FILE
git grep -E -l '__version__ =.*[0-9]+\.[0-9]+\.[0-9]+(\.[0-9]+)?(dev|rc*[0-9])?'  $FILES \
   | xargs -I FILE $SED '/__version__/s/[0-9]+\.[0-9]+\.[0-9]+(\.[0-9]+)?(dev|rc*[0-9])?/'${RELEASE}'/' -i.bak FILE
git status

