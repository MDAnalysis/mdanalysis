#!/usr/bin/env python3

# MIT License

# Copyright (c) 2023 MDAnalysis Development Team

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import argparse
import os
from urllib import request
import json
from github import Github


parser = argparse.ArgumentParser(
    description="Write PR comment for failed darker linting",
)


parser.add_argument(
    "--main_stat",
    type=str,
    help="Status of main code linting",
)


parser.add_argument(
    "--test_stat",
    type=str,
    help="Status of test code linting",
)


def bool_outcome(outcome: str) -> bool:
    return True if (outcome == 'success') else False


if __name__ == "__main__":
    args = parser.parse_args()

    git = Github(os.environ['GITHUB_TOKEN'])
    repo = git.get_repo("MDAnalysis/mdanalysis")

    run_id = os.environ['GITHUB_RUN_ID']
    job_id = os.environ['GITHUB_RUN_NUMBER']

    def get_pull_requests(repo):
        pulls = [pull for pull in repo.get_pulls()]
        # somehow get PR via the PR number here
        return pr

    def do_comment(pr):

        # check if the comment exists
        comments = [comm for comm in pr.get_comments() if "Linter Bot" in comm.body]
        if len(comments) > 0:
            # update -- note: probably should fail if there's more than 1 comment, but let's ignore this for now
            comments[0].edit(body="Linter Bot: this is an edit message")
        else:
            # add
            pr.create_issue_comment("Linter Bot: this is a new message")

    def get_job_run(repo, pr, run_id):
        lint_wkflow = [wf for wf in repo.get_workflows() if wf.name == 'linters'][0]
        run = [r for r in linters.get_runs(branch=pr.head.ref) if r.id == run_id][0]

        with request.urlopen(run.jobs_url) as url:
            data = json.load(url)

        for job in data['jobs']:
            if job['name'] == 'darker_lint':
                return job['html_url']

        return 'N/A'


    print(f"Linting - code: {bool_outcome(args.main_stat)}, "
          f"tests: {bool_outcome(args.test_stat)}, "
          f"action: https://www.github.com/MDAnalysis/mdanalysis/actions/runs/{run_id}/jobs/{job_id}")
