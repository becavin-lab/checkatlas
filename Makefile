.ONESHELL:
ENV_PREFIX=$(shell python -c "if __import__('pathlib').Path('.venv/bin/pip').exists(): print('.venv/bin/')")
USING_POETRY=$(shell grep "tool.poetry" pyproject.toml && echo "yes")
VERSION=$(shell poetry version | awk '{print $$2}')

.PHONY: help
help:             ## Show the help.
	@echo "Usage: make <target>"
	@echo ""
	@echo "Targets:"
	@fgrep "##" Makefile | fgrep -v fgrep


.PHONY: show
show:             ## Show the current environment.
	@echo "Current environment:"
	$(ENV_PREFIX)poetry env info

.PHONY: install
install:          ## Install the project in dev mode.
	@echo "Run checkatlas install - create poetrry virtual env"
	$(ENV_PREFIX)poetry install

.PHONY: fmt
fmt:              ## Format code using black & isort.
	@echo "Run project file formatting"
	$(ENV_PREFIX)poetry run isort .
	$(ENV_PREFIX)poetry run black -l 79 .

.PHONY: lint
lint:             ## Run pep8, black, mypy linters.
	@echo "Run project linting"
	$(ENV_PREFIX)poetry run flake8 checkatlas/
	$(ENV_PREFIX)poetry run flake8 tests/
	$(ENV_PREFIX)poetry run black -l 79 --check checkatlas/
	$(ENV_PREFIX)poetry run black -l 79 --check tests/
	$(ENV_PREFIX)poetry run mypy --ignore-missing-imports checkatlas/
	$(ENV_PREFIX)poetry run mypy --ignore-missing-imports tests/

.PHONY: test
test:             ## Run tests and generate coverage report.
	$(ENV_PREFIX)poetry run pytest -v --cov-config .coveragerc --cov=checkatlas -l --tb=short --maxfail=1 tests/
	$(ENV_PREFIX)poetry run coverage xml
	$(ENV_PREFIX)poetry run coverage html

.PHONY: watch
watch:            ## Run tests on every change.
	ls **/**.py | entr $(ENV_PREFIX)pytest -s -vvv -l --tb=long --maxfail=1 tests/

.PHONY: docs
docs:             ## Build the documentation.
	@echo "Building documentation ..."
	@$(ENV_PREFIX)poetry run mkdocs build
	URL="site/index.html"; open $$URL || xdg-open $$URL || sensible-browser $$URL || x-www-browser $$URL || gnome-open $$URL

.PHONY: ci
ci:          ## Run a continuous integration : Add every change to git and create a new tag for continuous integration.
	@echo "WARNING: You need first to add changes to git with git add"
	@echo "Push change to github and run continous integration scripts"
	@git add --all
	@git commit -m "Continuous integration 🔄 tests-$(VERSION)"
	@echo "creating git tag : tests-$(VERSION)"
	@git tag tests-$(VERSION)-12
	@git push -u origin HEAD --tags
	@echo "Github Actions will detect the new tag and run the continuous integration process."

.PHONY: release
release:          ## Create a new tag for release.
	@echo "WARNING: This operation will create a version tag and push to github"
	@echo "Reading version $(VERSION) from: pyproject.toml"
	@$(ENV_PREFIX)poetry run gitchangelog > HISTORY.md
	@git add HISTORY.md pyproject.toml
	@git commit -m "release: version $(VERSION) 🚀"
	@echo "creating git tag : release-$(VERSION)"
	@git tag release-$(VERSION)
	@git push -u origin HEAD --tags
	@echo "Github Actions will detect the new tag and release the new version."


.PHONY: clean
clean:            ## Clean unused files.
	@find ./ -name '*.pyc' -exec rm -f {} \;
	@find ./ -name '__pycache__' -exec rm -rf {} \;
	@find ./ -name 'Thumbs.db' -exec rm -f {} \;
	@find ./ -name '*~' -exec rm -f {} \;
	@find ./ -name '.nextflow.log*' -exec rm -f {} \;
	@rm -rf .pytest_cache
	@rm -rf .mypy_cache
	@rm -rf build
	@rm -rf dist
	@rm -rf *.egg-info
	@rm -rf htmlcov
	@rm -rf .tox/
	@rm -rf .nextflow/
	@rm -rf docs/_build
