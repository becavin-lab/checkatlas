"""Python setup.py for checkatlas package"""
import io
import os
from setuptools import find_packages, setup


def read(*paths, **kwargs):
    """Read the contents of a text file safely.
<<<<<<< HEAD
    read("checkatlas", "VERSION")
=======
    >>> read("checkatlas", "VERSION")
>>>>>>> b6ed53fdf89beca99eb88ccb03bd83ed951329f4
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
<<<<<<< HEAD
    description="project_description",
    url="https://github.com/becavin-lab/checkatlas",
=======
    description="Awesome checkatlas created by becavin-lab",
    url="https://github.com/becavin-lab/checkatlas/",
>>>>>>> b6ed53fdf89beca99eb88ccb03bd83ed951329f4
    long_description=read("README.md"),
    long_description_content_type="text/markdown",
    author="becavin-lab",
    packages=find_packages(exclude=["tests", ".github"]),
    install_requires=read_requirements("requirements.txt"),
    entry_points={
        "console_scripts": ["checkatlas = checkatlas.__main__:main"]
    },
    extras_require={"test": read_requirements("requirements-test.txt")},
)
