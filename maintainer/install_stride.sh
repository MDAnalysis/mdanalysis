#!/usr/bin/env bash
#
# Written by Lily Wang, 2020
#
# This script is placed into the Public Domain using the CC0 1.0 Public Domain
# Dedication https://creativecommons.org/publicdomain/zero/1.0/
#
# NOTE: IF YOU RUN THIS SCRIPT YOU MUST BE ABLE TO COMPLY WITH THE DSSP
#       LICENSE BELOW.
#
# All rights reserved, whether the	whole or part	of the program is
# concerned. Permission to use, copy, and modify	this software and its
# documentation is	granted	for academic use, provided that:
#
# i. this copyright notice appears in all copies of the	software and
# related documentation;
#
# ii. the reference given below (Frishman and Argos, 1995) must be
# cited in any publication of scientific results based in part or
# completely	on the use of the program;
#
# iii. bugs will be reported to the authors.
#
# The use of the  software in  commercial	 activities  is	 not  allowed
# without a prior written commercial license agreement.
#
# IF YOU USE STRIDE IN PUBLISHED WORK, PLEASE CITE:
#
# Frishman, D and Argos, P. (1995) Knowledge-based secondary structure
# assignment. Proteins: structure, function and genetics, 23,  
# 566-579.
#
# The source code is available at http://webclu.bio.wzw.tum.de/stride/
#
# arguments
#    OSNAME : linux | darwin
#    PREFIX : "path/to/installdir"
#
# PREFIX can be relative to CWD or an absolute path; it will be created
# and STRIDE unpacked under PREFIX/stride.
#
#


set -o errexit -o nounset

SCRIPT=$(basename $0)


DOWNLOAD_URL='http://webclu.bio.wzw.tum.de/stride/stride.tar.gz'
TARFILE='stride.tar.gz'


# path to dir with executables in the current HOLE distribution
STRIDE_EXE_DIR=stride

function die () {
    local msg="$1" err=$2
    echo "[${SCRIPT}] ERROR: $msg [$err]"
    exit $err
}

function is_native_executable () {
    local filename="$1" OSNAME="$2"
    file "${filename}" | grep -qi ${OFORMAT}
    return $?
}

OSNAME="$1"
case "$OSNAME" in
    Linux|linux)
	OSNAME=Linux
	OFORMAT=Linux
	;;
    OSX|osx|Darwin|darwin)
	OSNAME=Darwin
	OFORMAT="Mach-O"
	;;
    *)
	die "OS '${OSNAME}' not supported." 10;;
esac


PREFIX="$2"
test -n "$PREFIX" || die "missing argument PREFIX (installation directory)" 10
PREFIX="${PREFIX}/${STRIDE_EXE_DIR}"

#------------------------------------------------------------
# start installation
#------------------------------------------------------------

mkdir -p "$PREFIX" && cd "$PREFIX" || die "Failed to create and cd to $PREFIX" 1
if ! test -f ${TARFILE}; then
    echo "Downloading ${TARFILE} from ${DOWNLOAD_URL}..."
    # fixing curl on travis/anaconda, see http://stackoverflow.com/a/31060428/334357
    export CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt
    curl -L ${DOWNLOAD_URL} -o ${TARFILE} || \
    die "Failed to download ${TARFILE} from ${DOWNLOAD_URL}" 1
fi

# only install the executables in STRIDE_EXE_DIR 
tar xvf ${TARFILE} || die "Failed 'tar xvf ${TARFILE}'" 1
STRIDE_EXE="${PWD}/stride"

echo "Making"
make

test -f "${STRIDE_EXE}" || die "stride executable ${STRIDE_EXE} not installed" 2
is_native_executable ${STRIDE_EXE} ${OFORMAT} || \
    { file ${STRIDE_EXE}; \
      die "${STRIDE_EXE} will not run on ${OSNAME} (object format ${OFORMAT})" 3; }

echo "stride was installed into ${PWD}"
echo ">>> export PATH=\${PATH}:${PWD}"

# repeat this line in .travis.yml
export PATH=${PATH}:${PWD}
