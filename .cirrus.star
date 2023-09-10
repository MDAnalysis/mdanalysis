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

    # Some debugging to know what state you are in, by default everything is None but PR
    print(env.get("CIRRUS_TAG") == None)
    print(env.get("CIRRUS_RELEASE") == None)
    print(env.get("CIRRUS_PR") != None)
    print(env.get("CIRRUS_BASE_BRANCH") == "develop")

    # If it's a CRON job named twiceweekly
    if ((env.get("CIRRUS_CRON") == "twiceweekly") and (env.get("CIRRUS_BRANCH") == "develop")):
        return fs.read("maintainer/ci/cirrus-ci.yml")

    # If you've tagged a package or released something, deploy
    if ((env.get("CIRRUS_TAG") != None) or (env.get("CIRRUS_RELEASE") != None)):
        return fs.read("maintainer/ci/cirrus-deploy.yml")

    return []
