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
    return True if (outcome == 'passed') else False


if __name__ == "__main__":
    args = parser.parse_args()

    run_id = os.environ['GITHUB_RUN_ID']
    job_id = os.environ['GITHUB_RUN_NUMBER']

    def bool_outcome

    print(f"Linting - code: {bool_outcome(args.main_stat)}, "
          f"tests: {bool_outcome(args.test_stat)}, "
          f"action: https://www.github.com/MDAnalysis/mdanalysis/actions/runs/{run_id}/jobs/{job_id}")
