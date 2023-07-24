# ![CheckAtlas](docs/images/checkatlas_logo.png) 


![PyPI](https://img.shields.io/pypi/v/checkatlas)
![PyPI - Downloads](https://img.shields.io/pypi/dw/checkatlas)
![PyPI - License](https://img.shields.io/pypi/l/checkatlas)
![Conda](https://img.shields.io/conda/pn/bioconda/checkatlas)

[![codecov](https://codecov.io/gh/becavin-lab/checkatlas/branch/main/graph/badge.svg?token=checkatlas_token_here)](https://codecov.io/gh/becavin-lab/checkatlas)
[![CI](https://github.com/becavin-lab/checkatlas/actions/workflows/tests.yml/badge.svg)](https://github.com/becavin-lab/checkatlas/actions/workflows/tests.yml)
[![Documentation Status](https://readthedocs.org/projects/checkatlas/badge/?version=latest)](https://checkatlas.readthedocs.io/en/latest/?badge=latest)
[![Gitter](https://badges.gitter.im/checkatlas/checkatlas.svg)](https://app.gitter.im/#/room/!KpJcsVTOlGjwJgtLwF:gitter.im)

![Static Badge](https://img.shields.io/badge/Packaging-Poetry-blue)
![Static Badge](https://img.shields.io/badge/Docs-Mkdocs-red)
![Static Badge](https://img.shields.io/badge/Linting-flake8%20black%20mypy-yellow)

CheckAtlas is a one liner tool to check the quality of your single-cell atlases. For every atlas, it produces the
quality control tables and figures which can be then processed by multiqc. CheckAtlas is able to check the quality of Scanpy, Seurat,
and CellRanger files.

More information on the [read the doc page](https://checkatlas.readthedocs.io/en/latest/)


## Summary

Powered by nextflow, checkatlas can be ran in one command line:

```bash
nextflow run nf-core-checkatlas --path search_folder/
```

The checkatlas workflow start with a fast crawl through your working directory. It detects Seurat (.rds), Scanpy (.h5ad) or cellranger (.h5) atlas files.

Then, it goes through all atlas files and produce summary information:

- All basic QC (nRNA, nFeature, ratio_mito)
- General information (nbcells, nbgenes, nblayers)
- All elements in atlas files (obs, obsm, uns, var, varm)
- Reductions (pca, umap, tsne)
- All metrics (clustering, annotation, dimreduction, specificity)

All tables and figs are saved in the **checkatlas_files** folder in your search folder.

A single html report is produced, using MultiQC, in **checkatlas_files/Checkatlas-MultiQC.html**.

![Checkatlas workflow](docs/checkatlas_workflow.png)


## Examples

- Evaluate and compare different scanpy atlases:
[Example 1](https://checkatlas.readthedocs.io/en/latest/examples/CheckAtlas_example_1/Checkatlas_MultiQC.html)

- Evaluate different version of one atlas:
[Example 2](https://checkatlas.readthedocs.io/en/latest/examples/CheckAtlas_example_2/Checkatlas_MultiQC.html)

- Evaluate Scanpy, Seurat and CellRanger objects in your folder:
[Example 3](https://checkatlas.readthedocs.io/en/latest/examples/CheckAtlas_example_3/Checkatlas_MultiQC.html)

- Evaluate an integrated Scanpy atlas with the corresponding raw CellRanger atlases:
[Example 4](https://checkatlas.readthedocs.io/en/latest/examples/CheckAtlas_example_4/Checkatlas_MultiQC.html)

- Evaluate different Cellranger atlases with multiple chemistry version and cellranger version:
[Example 5](https://checkatlas.readthedocs.io/en/latest/examples/CheckAtlas_example_5/Checkatlas_MultiQC.html)


## Installation

CheckAtlas is in two parts. The checkatlas pythn module which can be downloaded with PyPi, and the checkatlas workflow which can be downloaded with nextflow.

```bash
pip install checkatlas
```

```bash
nextflow pull becavin-lab/nf-core-checkatlas
```

You need also to install a version of MultiQC with checkatlas capability (for the moment). This version of MultiQC is available at checkatlas branch of github.com:becavin-lab/MultiQC.

```bash
git clone git@github.com:becavin-lab/MultiQC.git
cd MultiQC/
git checkout checkatlas
pip install .
```

Finally, checkatlas comes with rpy2 to perform the interface between python and R. But, it does not automatically install Seurat. So if you want to screen Seurat atlases you need to perfrom this last installation

```bash
% R
> install.packages('Seurat')
> library(Seurat)
```


## Development

This project is in a very early development phase. All helpers are welcome. Please contact us or submit an issue.

Read the [CONTRIBUTING.md](docs/contributing.md) file.


Project developed thanks to the project template : (https://github.com/rochacbruno/python-project-template/)

