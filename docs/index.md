# CheckAtlas

![PyPI](https://img.shields.io/pypi/v/checkatlas)
![PyPI - Downloads](https://img.shields.io/pypi/dw/checkatlas)
![PyPI - License](https://img.shields.io/pypi/l/checkatlas)
![Conda](https://img.shields.io/conda/pn/bioconda/checkatlas)

[![codecov](https://codecov.io/gh/becavin-lab/checkatlas/branch/main/graph/badge.svg?token=checkatlas_token_here)](https://codecov.io/gh/becavin-lab/checkatlas)
[![CI](https://github.com/becavin-lab/checkatlas/actions/workflows/tests.yml/badge.svg)](https://github.com/becavin-lab/checkatlas/actions/workflows/tests.yml)
[![Documentation Status](https://readthedocs.org/projects/checkatlas/badge/?version=latest)](https://checkatlas.readthedocs.io/en/latest/?badge=latest)
[![Gitter](https://badges.gitter.im/checkatlas/checkatlas.svg)](https://app.gitter.im/#/room/!KpJcsVTOlGjwJgtLwF:gitter.im)


CheckAtlas is a one liner tool to check the quality of your single-cell atlases. For every atlas, it produces the
quality control tables and figures which can be then processed by multiqc. CheckAtlas is able to load Scanpy, Seurat,
and CellRanger files.


## Summary

### Parse Scanpy, Seurat and CellRanger objects

The checkatlas pipeline start with a fast crawl through your working directory. It detects Seurat (.rds), Scanpy (.h5ad) or cellranger (.h5) atlas files.

### Create checkatlas summary files

   Go through all atlas files and produce summary information:

   - All basic QC (nRNA, nFeature, ratio_mito)
   - General information (nbcells, nbgenes, nblayers)
   - All elements in atlas files (obs, obsm, uns, var, varm)
   - Reductions (pca, umap, tsne)
   - All metrics (clustering, annotation, dimreduction, specificity)

### Parse checkatlas files in MultiQC

   Update MultiQC project to add checkatlas parsing. Dev project in: https://github.com/becavin-lab/MultiQC/tree/checkatlas

https://checkatlas.readthedocs.io/en/latest/

