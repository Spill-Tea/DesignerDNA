# PyTemplate
[![build status][buildstatus-image]][buildstatus-url]

[buildstatus-image]: https://github.com/Spill-Tea/PyTemplate/actions/workflows/python-app.yml/badge.svg?branch=main
[buildstatus-url]: https://github.com/Spill-Tea/PyTemplate/actions?query=branch%3Amain

Python Project Template. Be sure to create a template directly
from github.

<!-- omit in toc -->
## Table of Contents
- [PyTemplate](#pytemplate)
  - [Installation](#installation)
  - [For Developers](#for-developers)
  - [License](#license)


## Installation
Clone the repository and pip install.

```bash
git clone https://github.com/Spill-Tea/PyTemplate.git
cd PyTemplate
pip install .
```

Alternatively, you may install directly from github.
```bash
pip install git+https://github.com/Spill-Tea/PyTemplate@main
```


## For Developers
After cloning the repository, create a new virtual environment and run the following commands:

```bash
pip install -e ".[dev]"
pre-commit install
pre-commit run --all-files
```

Running unit tests locally is straightforward with tox. Make sure
you have all python versions available required for your project
The `p` flag is not required, but it runs tox environments in parallel.
```bash
tox -p
```
Be sure to run tox before creating a pull request.

## License
[BSD-3](LICENSE)
