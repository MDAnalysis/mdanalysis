# See https://cirrus-ci.org/guide/programming-tasks/ for more information on
# starlark CirrusCI files
# Inspired by scipy's .cirrus.star script
# In the spirit of ensuring this can also be freely used in derviative works by
# the community, we release the contents of this file under the MIT license.


load("cirrus", "env", "fs")


def main(ctx):
    # TODO: Add a summary of what this file does
    if env.get("CIRRUS_REPO_FULL_NAME") != "MDAnalysis/mdanalysis":
        return []

    fs.read("maintainer/ci/cirrus-ci.yml")
