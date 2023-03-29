#!/usr/bin/env python3

# MIT License

# Copyright (c) 2023 Irfan Alibay

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
    "--json",
    type=str,
    help="Input JSON file with status results",
)


def get_pull_request(repo, pr_num):
    """
    Simple method to get a PyGithub PR object from a PR number

    Parameters
    ----------
    repo : PyGithub.Repository
      Github Repository API object.
    pr_num : int
      Pull request number.

    Returns
    -------
    PyGithub.PullRequest
      Pull Request object corresponding to the input PR number.
    """
    # Can get PR directly from PR number
    return repo.get_pull(pr_num)


def get_action_url(repo, pr, run_id, workflow_name, job_name):
    """
    Mix PyGithub & Github V3 REST API method to extract the url path to a
    github actions workflow job corresponding to a given GITHUB_RUN_ID.

    Parameters
    ----------
    repo : PyGithub.Repository
      Github Repository API object.
    run_id : str
      Github actions RUN ID as defined by the github actions environment
      variable `GITHUB_RUN_ID`.
    workflow_name : str
      Name of the workflow to extract a job url for.
    job_name : str
      Name of the job within the workflow to extract a url for.


    Returns
    -------
    str
      URL to github actions workflow job or 'N/A' if a corresponding job name
      could not be found.
    """
    # Accessing get_workflow directly currently fails when passing a name
    # Instead do the roundabout way by getting list of all workflows
    linters = [wf for wf in repo.get_workflows()
               if wf.name == workflow_name][0]

    # Extract the gh action run
    run = [r for r in linters.get_runs(branch=pr.head.ref)
           if r.id == int(run_id)][0]

    # The exact job url can't be recovered via the Python API
    # Switch over to using the REST API instead
    with request.urlopen(run.jobs_url) as url:
        data = json.load(url)

    for job in data['jobs']:
        if job['name'] == job_name:
            return job['html_url']

    return 'N/A'


def bool_outcome(outcome):
    """
    Converts github action job status outcome to bool.

    Parameters
    ----------
    outcome : str
      Github action job step outcome message.

    Returns
    -------
    bool
      Whether or not the job step was successful.
    """
    return True if (outcome == 'success') else False


def gen_message(pr, main_stat, test_stat, action_url):
    """
    Generate a user facing message on the status of the darker linting
    action.

    Parameters
    ----------
    pr : PyGithub.PullRequest
      Pull Request object representing the target for this message.
    main_stat : str
      Outcome of darker linting of main package code.
    test_stat : str
      Outcome of darker linting of testsuite code.
    action_url : str
      URL pointing to darker linting job log.


    Returns
    -------
    str
      Message to be posted to PR author.
    """

    def _format_outcome(stat):
        if bool_outcome(stat):
            return "✅ Passed"
        else:
            return "⚠️  Possible failure"

    msg = ('### Linter Bot Results:\n\n'
           f'Hi @{pr.user.login}! Thanks for making this PR. '
           'We linted your code and found the following: \n\n')

    # If everything is ok
    if bool_outcome(main_stat) and bool_outcome(test_stat):
        msg += ('There are currently no issues detected! 🎉')
    else:
        msg += ('Some issues were found with the formatting of your code.\n'
                '| Code Location | Outcome |\n'
                '| --- | --- |\n'
                f'| main package | {_format_outcome(main_stat)}|\n'
                f'| testsuite | {_format_outcome(test_stat)}|\n'
                '\nPlease have a look at the `darker-main-code` and '
                '`darker-test-code` steps here for more details: '
                f'{action_url}\n\n'
                '---\n'
                '_**Please note:** The `black` linter is purely '
                'informational, you can safely ignore these outcomes if '
                'there are no flake8 failures!_')
    return msg


def post_comment(pr, message, match_string):
    """
    Post a comment in a Pull Request.

    If a comment with text matching `match_string` is found in the
    Pull Request, the comment will be edited.

    Parameters
    ----------
    pr : PyGithub.PullRequest
      Pull Request object representing the target for this message.
    message : str
      The message to post as a comment.
    match_string : str
      A matching string to recognise if the comment already exists.
    """
    # Extract a list of matching comments from PR
    comments = [comm for comm in pr.get_issue_comments() if match_string in comm.body]

    if len(comments) > 0:
        # Edit the comment in-place
        # By default we assume that the bot can write faster than anyone else
        comments[0].edit(body=message)
    else:
        # If the comment doesn't exist, create one
        pr.create_issue_comment(message)


if __name__ == "__main__":
    args = parser.parse_args()

    git = Github(os.environ['GITHUB_TOKEN'])
    repo = git.get_repo("MDAnalysis/mdanalysis")

    with open(args.json, 'r') as f:
        status = json.load(f)

    run_id = status['RUN_ID']
    print(f"debug run_id: {run_id}")

    # Get Pull Request
    pr_num = int(status['PR_NUM'])
    print(f"debug pr_num: {pr_num}")
    pr = get_pull_request(repo, pr_num)

    # Get the url to the github action job being pointed to
    action_url = get_action_url(repo, pr, run_id,
                                workflow_name='linters',
                                job_name='darker_lint')

    # Get the message you want to post to users
    with open(args.json, 'r') as f:
        results_dict = json.load(f)

    message = gen_message(pr,
                          status['main_stat'],
                          status['test_stat'],
                          action_url)

    # Post your comment
    post_comment(pr, message, match_string='Linter Bot Results:')
