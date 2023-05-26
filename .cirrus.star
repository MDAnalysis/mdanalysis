# See https://cirrus-ci.org/guide/programming-tasks/ for more information on
# starlark CirrusCI files
# Inspired by scipy's .cirrus.star script https://github.com/scipy/scipy/blob/main/.cirrus.star
# In the spirit of ensuring this can also be freely used in derviative works by
# the community, we release the contents of this file under the MIT license.


load("cirrus", "env", "fs")


def main(ctx):
    # Default case: don't do anything if not in the core repo
    if env.get("CIRRUS_REPO_FULL_NAME") != "MDAnalysis/mdanalysis":
        return []

    print(env.get("CIRRUS_TAG") == None)
    print(env.get("CIRRUS_RELEASE") == None)
    print(env.get("CIRRUS_PR") != None)

    #if (env.get("CIRRUS_BASE_BRANCH" == "develop") and (env.get("CIRRUS_PR") != ""):
    #    fs.read("maintainer/ci/cirrus-ci.yml")

    #if (env.get("CIRRUS_TAG" != "") or (env.get("CIRRUS_RELEASE") != ""):
    #    fs.read("maintainer/ci/cirrus-deploy.yml")

    #return []
    return fs.read("maintainer/ci/cirrus-deploy.yml")
