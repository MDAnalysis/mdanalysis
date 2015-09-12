#!/bin/bash
# Deploying docs from travis-ci.
# See https://github.com/MDAnalysis/mdanalysis/issues/386
# Script based on https://github.com/steveklabnik/automatically_update_github_pages_with_travis_example

# Run this script from the top-level of the checked out git
# repository. A github OAuth token must be available in the evironment
# variable GH_TOKEN and is set up through the .travis.yml
# env:global:secure parameter (encrypted with travis-ci's public key)/

set -o errexit -o nounset

function die () {
    local msg="$1" err={$2:-1}
    echo "ERROR: $msg [$err]"
    exit $err
}

DOCDIR="package/doc/html"
REPOSITORY="github.com/MDAnalysis/mdanalysis.git"

rev=$(git rev-parse --short HEAD)

test -n "${GH_TOKEN}" || die "GH_TOKEN is empty: need OAuth GitHub token to continue" 100
cd $DOCDIR || die "Failed to 'cd $DOCDIR'. Run from the top level of the repository"

git init
git config user.name "Travis CI"
git config user.email "TravisCI@mdanalysis.org"

git remote add upstream "https://${GH_TOKEN}@${REPOSITORY}"
git fetch upstream
git reset upstream/gh-pages

touch .
touch .nojekyll

git add -A .
git commit -m "rebuild html docs with sphinx at ${rev}"
git push upstream HEAD:gh-pages


