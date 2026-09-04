## Contributing to MDAnalysis

Thanks for contributing to MDAnalysis!

All members of the MDAnalysis community adhere to our [Code of Conduct](https://github.com/MDAnalysis/mdanalysis/?tab=coc-ov-file). By contributing code and interacting with us on GitHub, the forums, or by any other means you consent to follow the Code of Conduct.

#### Reporting issues

If you've found a defect with MDAnalysis we'd love to know so we can fix it.  Please follow the Issue template so we can quickly diagnose the problem, in particular the piece of code that causes the problem.

If your issue isn't a defect with the code and instead you require help using MDAnalysis, drop by in [discord #users or GitHub discussions](https://www.mdanalysis.org/community/#ask-questions--get-help).

#### Contributing code

We welcome your code contributions! Please check out the [Contributing](https://www.mdanalysis.org/contribute/) section on our web page for how you can contribute to the whole MDAnalysis project and in particular the [User Guide's section on Contributing](https://userguide.mdanalysis.org/stable/contributing_code.html) for the MDAnalysis package in particular. Please also familiarize yourself with the current version of our [AI Policy](https://github.com/MDAnalysis/mdanalysis/blob/develop/AI_POLICY.md)

MDAnalysis devs are most easily reached in the [discord #developers channel](https://www.mdanalysis.org/community/#ask-questions--get-help).

#### Local quality checks

Before opening a pull request, please run the same style checks locally using pre-commit:

```bash
python -m pip install pre-commit
pre-commit install
pre-commit run --all-files
```

Running these checks before pushing helps keep pull requests focused and reduces CI turnaround time.

