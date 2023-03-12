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
    "--run_id",
    type=str,
    help="Github action run id",
)


parser.add_argument(
    "--job_id",
    type=str,
    help="Github action job id",
)


parser.add_argument(
    "--main_stat",
    type=bool,
    help="Status of main code linting",
)


parser.add_argument(
    "--test_stat",
    type=bool,
    help="Status of test code linting",
)


if __name__ == "__main__":
    args = parser.parse_args()

    print(f"Linting - code: {args.main_stat}, tests: {args.test_stat}, "
          f"action: https://www.github.com/MDAnalysis/mdanalysis/actions/runs/{args.run_id}/jobs/{args.job_id}")
