#!/bin/bash
# Oliver Beckstein 2014, Placed into the Public Domain
# run setup.py develop for package and tests
usage="$(basename $0) [options]

Run 'python setup.py' in the package and testsuite directory to install
the library and the unit tests.

By default does a standard 'install'. Use -d for a 'develop' installation
that uses the sources.

Options

-h       help
-d       develop install
"

die () {
    local msg="$1" err="${2:-1}"
    echo 1>&2 "EE ERROR ($err): $msg"
    exit $err
}

setup_py () {
    local command="${1}" dir="$2"
    pushd "$dir" || die "failed to cd to '$dir'" $?
    echo ">>     cd $dir"
    echo ">>     python setup.py ${command}"
    python setup.py ${command}
    popd
}

command="install"
while getopts hd OPT; do
    case $OPT in 
	h) echo "$usage"; exit 0;;
	d) command="develop";;
    esac
done

setup_py $command package
setup_py $command testsuite


