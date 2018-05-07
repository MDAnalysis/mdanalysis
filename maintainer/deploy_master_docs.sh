#!/bin/bash
# Deploying docs from travis-ci.
# See https://github.com/MDAnalysis/mdanalysis/issues/386
# Script based on https://github.com/steveklabnik/automatically_update_github_pages_with_travis_example

# Run this script from the top-level of the checked out git
# repository. A github OAuth token must be available in the evironment
# variable GH_TOKEN and is set up through the .travis.yml
# env:global:secure parameter (encrypted with travis-ci's public key)/
#
# Additional environment variables set in .travis.yml
#  GH_REPOSITORY     repo to full from and push to
#  GH_DOC_BRANCH     branch from which the docs are built
#  GIT_CI_USER       name of the user to push docs as
#  GIT_CI_EMAIL      email of the user to push docs as
#  MDA_DOCDIR        path to the docdir from top of repo
#
# NOTE: If any of these environment variables are not set or 
#       empty then the script will exit with and error (-o nounset).

# TESTING TO PUSH TO THE docs repo
# Assume that this is MANUALLY run from top level of repo
# and the docs were built with
#
#  (cd package && python setup.py build_sphinx)
#

# The release sitemap.xml is generated from the development
# sitemap.xml by changing MDA_DEVELOPMENT_DOCS_URL --> MDA_RELEASE_DOCS_URL
#
# See comment in package/doc/sphinx/source/conf.py for site_url.
#
MDA_RELEASE_DOCS_URL="https://www.mdanalysis.org/docs/"
MDA_DEVELOPMENT_DOCS_URL="https://www.mdanalysis.org/mdanalysis/"

GH_REPOSITORY=github.com/MDAnalysis/docs.git
#MDA_DOCDIR=${TRAVIS_BUILD_DIR}/package/doc/html/html
MDA_DOCDIR=package/doc/html/html
GIT_CI_USER="TravisCI"
GIT_CI_EMAIL="TravisCI@mdanalysis.org"

# for informational purposes at the moment
GH_DOC_BRANCH=$(git rev-parse --abbrev-ref HEAD)


set -o errexit -o nounset

function die () {
    local msg="$1" err=${2:-1}
    echo "ERROR: $msg [$err]"
    exit $err
}

rev=$(git rev-parse --short HEAD)
rootdir="$(git rev-parse --show-toplevel)"

maintainer_dir="${rootdir}/maintainer"

# the following tests should be superfluous because of -o nounset
test -n "${GH_TOKEN}" || die "GH_TOKEN is empty: need OAuth GitHub token to continue" 100
test -n "${GH_REPOSITORY}" || die "GH_REPOSITORY must be set in .travis.yml" 100
test -n "${MDA_DOCDIR}" || die "MDA_DOCDIR must be set in .travis.yml" 100

cd ${MDA_DOCDIR} || die "Failed to 'cd ${MDA_DOCDIR}'. Run from the top level of the repository"

# for local testing
test -d .git && rm -rf .git

git init
git config user.name "${GIT_CI_USER}"
git config user.email "${GIT_CI_EMAIL}"

git remote add docs "https://${GH_TOKEN}@${GH_REPOSITORY}"
git fetch --depth 10 docs master
git reset docs/master

touch .
touch .nojekyll

git add -A .
git commit -m "rebuilt html docs from branch ${GH_DOC_BRANCH} with sphinx at MDAnalysis/mdanalysis@${rev}"

# fix sitemap.xml for release docs
${maintainer_dir}/adapt_sitemap.py --search ${MDA_DEVELOPMENT_DOCS_URL} --replace ${MDA_RELEASE_DOCS_URL} -o sitemap.xml sitemap.xml

git add sitemap.xml
git commit -m "adjusted sitemap.xml for release docs at ${MDA_RELEASE_DOCS_URL}" 

git push -q docs HEAD:master


