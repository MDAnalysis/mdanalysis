#!/bin/bash
#
# Written by Oliver Beckstein, 2016
#
# This script is placed into the Public Domain using the CC0 1.0 Public Domain
# Dedication https://creativecommons.org/publicdomain/zero/1.0/
#
#
# NOTE: IF YOU RUN THIS SCRIPT YOU MUST BE ABLE TO COMPLY WITH THE HOLE
#       NOT-FOR-PROFIT LICENSE.
# 
#
# Install HOLE from http://www.holeprogram.org/
#
# arguments
#    OSNAME : linux | darwin
#    PREFIX : "path/to/installdir"
#
# PREFIX can be relative to CWD or an absolute path; it will be created
# and HOLE unpacked under PREFIX (the tar file contains "hole2/...")
#
#
# HOLE v2.2004 is used under the terms of the 'HOLE END USER LICENCE AGREEMENT
# NOT-FOR-PROFIT VERSION' as provided in the installation tarball as file
# 'doc/Licence-not-for-profit.asciidoc' and see copy at
# https://github.com/MDAnalysis/mdanalysis/files/372246/Licence-not-for-profit.txt
# (recent as of 2016-07-19).
#


set -o errexit -o nounset

SCRIPT=$(basename $0)

# We still use the 2.2004 versions (academic only) because
# for the v2.2005 (Apache license) there are no pre-compiled binaries
# available at http://www.holeprogram.org/ as of 2018-02-22

# could switch to http://www.holeprogram.org/downloads/2.2.004/hole2-NotForProfit-2.2.004-Linux-i686.tar.gz
DOWNLOAD_URL_LINUX='https://www.holeprogram.org/downloads/2.2.004/hole2-NotForProfit-2.2.004-Linux-x86_64.tar.gz'
TARFILE_LINUX='hole2-NotForProfit-2.2.004-Linux-x86_64.tar.gz'

# could switch to http://www.holeprogram.org/downloads/2.2.004/hole2-NotForProfit-2.2.004-Darwin-i386.tar.gz
DOWNLOAD_URL_DARWIN='https://www.holeprogram.org/downloads/2.2.004/hole2-NotForProfit-2.2.004-Darwin-i386.tar.gz'
TARFILE_DARWIN=hole2-NotForProfit-2.2.004-Darwin-i386.tar.gz

# path to dir with executables in the current HOLE distribution
HOLE_EXE_DIR=hole2/exe

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
	DOWNLOAD_URL="${DOWNLOAD_URL_LINUX}"
	TARFILE=${TARFILE_LINUX}
    CURLARGS=""
	;;
    OSX|osx|Darwin|darwin)
	OSNAME=Darwin
	OFORMAT="Mach-O"
	DOWNLOAD_URL="${DOWNLOAD_URL_DARWIN}"
	TARFILE=${TARFILE_DARWIN}
    security export -p -t certs -k `security list-keychains -d system|cut -d '"' -f 2` -o certs.crt
    CURLARGS="--cacert $(PWD)/certs.crt"
	;;
    *)
	die "OS '${OSNAME}' not supported." 10;;
esac

PREFIX="$2"
test -n "$PREFIX" || die "missing second argument PREFIX (installation directory)" 10

#------------------------------------------------------------
# start installation
#------------------------------------------------------------

mkdir -p "$PREFIX" && cd "$PREFIX" || die "Failed to create and cd to $PREFIX" 1
if ! test -f ${TARFILE}; then
    echo "Downloading ${TARFILE} from ${DOWNLOAD_URL}..."
    # fixing curl on travis/anaconda, see http://stackoverflow.com/a/31060428/334357
    export CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt
    curl -L ${DOWNLOAD_URL} -o ${TARFILE} $CURLARGS || \
    die "Failed to download ${TARFILE} from ${DOWNLOAD_URL}" 1
fi

# only install the executables in HOLE_EXE_DIR 
tar xvf ${TARFILE} ${HOLE_EXE_DIR} || die "Failed 'tar xvf ${TARFILE} ${HOLE_EXE_DIR}'" 1

MDA_HOLE_BINDIR="${PWD}/${HOLE_EXE_DIR}"
HOLE_EXE="${MDA_HOLE_BINDIR}/hole"

test -d "${MDA_HOLE_BINDIR}" || { \
    echo "Content of ${PWD}:";
    ls -la;
    die "no HOLE exe dir ${MDA_HOLE_BINDIR} in $PWD" 2; }
test -f "${HOLE_EXE}" || die "hole executable ${HOLE_EXE} not installed" 2
is_native_executable ${HOLE_EXE} ${OFORMAT} || \
    { file ${HOLE_EXE}; \
      die "${HOLE_EXE} will not run on ${OSNAME} (object format ${OFORMAT})" 3; }

echo "HOLE executables were installed into ${MDA_HOLE_BINDIR}"
echo ">>> export PATH=\${PATH}:${MDA_HOLE_BINDIR}"

# repeat this line in .travis.yml
export PATH=${PATH}:${MDA_HOLE_BINDIR}
