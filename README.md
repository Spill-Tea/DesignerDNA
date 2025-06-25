# DesignerDNA
[![build status][buildstatus-image]][buildstatus-url]

[buildstatus-image]: https://github.com/Spill-Tea/DesignerDNA/actions/workflows/python-app.yml/badge.svg?branch=main
[buildstatus-url]: https://github.com/Spill-Tea/DesignerDNA/actions?query=branch%3Amain

DesignerDNA - Design DNA sequences with intent.

<!-- omit in toc -->
## Table of Contents
- [DesignerDNA](#designerdna)
  - [Installation](#installation)
  - [For Developers](#for-developers)
  - [License](#license)


## Installation
Clone the repository and pip install.

```bash
git clone https://github.com/Spill-Tea/DesignerDNA.git
cd DesignerDNA
pip install .
```

Alternatively, you may install directly from github.
```bash
pip install git+https://github.com/Spill-Tea/DesignerDNA@main
```


## For Developers
After cloning the repository, create a new virtual environment and run the following
commands:

```bash
pip install -e ".[dev]"
pre-commit install
pre-commit run --all-files
```

Running unit tests locally is straightforward with tox. Make sure you have all python
versions available required for your project. The `p` flag is not required, but it runs
tox environments in parallel.
```bash
tox -p
```
Be sure to run tox before creating a pull request.

## License
[BSD-3](LICENSE)
