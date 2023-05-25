# See https://cirrus-ci.org/guide/programming-tasks/ for more information on
# starlark CirrusCI files
# Inspired by scipy's .cirrus.star script https://github.com/scipy/scipy/blob/main/.cirrus.star
# In the spirit of ensuring this can also be freely used in derviative works by
# the community, we release the contents of this file under the MIT license.


load("cirrus", "env", "fs")


def main(ctx):
    # Default case: don't do anything if not in the core repo
    # or if you're not targetting the develop branch
    if ((env.get("CIRRUS_REPO_FULL_NAME") != "MDAnalysis/mdanalysis")
        or (env.get("CIRRUS_BASE_BRANCH") != "develop")):
        return []

    return fs.read("maintainer/ci/cirrus-deploy.yml")
    #return fs.read("maintainer/ci/cirrus-ci.yml") + fs.read("maintainer/ci/cirrus-deploy.yml")
