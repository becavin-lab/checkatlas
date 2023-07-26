# How to develop on this project

CheckAtlas welcomes contributions from the community.

**You need PYTHON3!**

This project is in a very early development phase. All helpers are welcome. Please contact us or submit an issue.

Read the [CONTRIBUTING.md](docs/contributing.md) file.

Checkatlas has two repositories:
- [The checkatlas module](https://github.com/becavin-lab/checkatlas)
- [The checkatlas nextflow workflow](https://github.com/becavin-lab/nf-core-checkatlas)

It has a module on MultiQC
- [MultiQC checkatlas branch](https://github.com/becavin-lab/MultiQC)

The checkatlas package is available on PyPI
- [Checkatlas PyPI](https://pypi.org/project/checkatlas/)

The bioconda recipe has been submitted
- [Checkatlas bioconda recipe](https://github.com/drbecavin/bioconda-recipes)


These instructions are for linux base systems. (Linux, MacOS, BSD, etc.)

## Setting up your own fork of this repo.

- On github interface click on `Fork` button.
- Clone your fork of this repo. `git clone git@github.com:YOUR_GIT_USERNAME/checkatlas.git`
- Enter the directory `cd checkatlas`
- Add upstream repo `git remote add upstream https://github.com/becavin-lab/checkatlas`

## Setting up your own environment

Checkatlas relies on poetry packaging tool.
Install poetry and run `make virtualenv` to create a virtual environment.
Then activate it with `source .venv/bin/activate`.

## Install the project in develop mode

Run `make install --with dev` to install the project in develop mode.

## Run the tests to ensure everything is working

Run `make test` to run the tests.

## Create a new branch to work on your contribution

Run `git checkout -b my_contribution`

## Make your changes

Edit the files using your preferred editor. (we recommend VSCode)

## Format the code

Run `make fmt` to format the code.

## Run the linter

Run `make lint` to run the linter.

## Test your changes

Run `make test` to run the tests.

Ensure code coverage report shows `100%` coverage, add tests to your PR.

## Build the docs locally

Run `make docs` to build the docs.

Ensure your new changes are documented.

## Commit your changes

This project uses [conventional git commit messages](https://www.conventionalcommits.org/en/v1.0.0/).

Example: `fix(package): update setup.py arguments 🎉` (emojis are fine too)

## Push your changes to your fork

Run `git push origin my_contribution`

## Submit a pull request

On github interface, click on `Pull Request` button.

Wait CI to run and one of the developers will review your PR.

## Makefile utilities

This project comes with a `Makefile` that contains a number of useful utility.

```bash 
❯ make
Usage: make <target>

Targets:
help:             ## Show the help.
install:          ## Install the project in dev mode.
fmt:              ## Format code using black & isort.
lint:             ## Run pep8, black, mypy linters.
test: lint        ## Run tests and generate coverage report.
watch:            ## Run tests on every change.
clean:            ## Clean unused files.
virtualenv:       ## Create a virtual environment.
release:          ## Create a new tag for release.
docs:             ## Build the documentation.
switch-to-poetry: ## Switch to poetry package manager.
init:             ## Initialize the project based on an application template.
```
