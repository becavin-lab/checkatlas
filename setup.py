"""Python setup.py for checkatlas package"""
import io
import os

from setuptools import find_packages, setup


def read(*paths, **kwargs):
    """Read the contents of a text file safely.
    read("checkatlas", "VERSION")
    '0.1.0'
    """

    content = ""
    with io.open(
        os.path.join(os.path.dirname(__file__), *paths),
        encoding=kwargs.get("encoding", "utf8"),
    ) as open_file:
        content = open_file.read().strip()
    return content


def read_requirements(path):
    return [
        line.strip()
        for line in read(path).split("\n")
        if not line.startswith(('"', "#", "-", "git+"))
    ]


setup(
    name="checkatlas",
    version=read("checkatlas", "VERSION"),
    description="One liner tool to check the quality of your single-cell atlases",
    url="https://github.com/becavin-lab/checkatlas/",
    long_description=read("README.md"),
    long_description_content_type="text/markdown",
    author="becavin-lab",
    packages=find_packages(exclude=["tests", ".github", "*dask-worker-space*"]),
    package_data={"checkatlas": ["checkatlas/convertSeurat.R"]},
    include_package_data=True,
    install_requires=read_requirements("requirements.txt"),
    entry_points={"console_scripts": ["checkatlas = checkatlas.__main__:main"]},
    extras_require={"test": read_requirements("requirements-test.txt")},
)
