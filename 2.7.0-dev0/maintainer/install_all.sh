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
-d       develop install (MODE = develop)
-u       user install (setup.py MODE --user)
"

die () {
    local msg="$1" err="${2:-1}"
    echo 1>&2 "EE ERROR ($err): $msg"
    exit $err
}

setup_py () {
    local command="${1}" dir="$2" options="$3"
    pushd "$dir" || die "failed to cd to '$dir'" $?
    echo ">>     cd $dir"
    echo ">>     python setup.py ${command} ${options}"
    python setup.py ${command} ${options}
    popd
}

command="install"
options=""
while getopts hdu OPT; do
    case $OPT in 
	h) echo "$usage"; exit 0;;
	d) command="develop";;
	u) options="${options} --user";;
	'?') die "Unknown option. Try -h for help." 2;;
    esac
done

setup_py $command package "${options}"
setup_py $command testsuite "${options}"


