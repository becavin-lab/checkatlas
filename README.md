
# CheckAtlas

[![codecov](https://codecov.io/gh/becavin-lab/checkatlas/branch/main/graph/badge.svg?token=checkatlas_token_here)](https://codecov.io/gh/becavin-lab/checkatlas)
[![CI](https://github.com/becavin-lab/checkatlas/actions/workflows/main.yml/badge.svg)](https://github.com/becavin-lab/checkatlas/actions/workflows/main.yml)
[![Documentation Status](https://readthedocs.org/projects/checkatlas/badge/?version=latest)](https://checkatlas.readthedocs.io/en/latest/?badge=latest)

CheckAtlas is a one liner tool to check the quality of your single-cell atlases. For every atlas, it produces the
quality control tables and figures which can be then processed by multiqc. CheckAtlas is able to load Scanpy, Seurat,
and CellRanger files.

## Summary

1. Parse Scanpy, Seurat and CellRanger objects
    
    CheckAtlas should be able to load : .rds, .h5 and .h5ad corresponding to single-cell experiment. Need to implement :
      - automatic conversion of Seurat object to Scanpy with SeuratDisk
      - Rapid check-up of files to see if a Seurat or Scanpy can be found
      - Automatic search in Scanpy files of key information = raw data, normalized data, integrated data, reductions, layers, assays, metadatas, etc...


2. Create checkatlas summary files
  
    Go through all Scanpy files and extract summary information. We won't to extract :

      - All basic QC (nRNA, nFeature, ratio_mito)
      - General information (nbcells, nbgenes, nblayers)
      - All elements in scanpy objects (obs, obsm, uns, var, varm)
      - Reductions (pca, umap, tsne)
      - All metrics (clustering, annotation, dimreduction, specificity)

3. Parse checkatlas files in MultiQC
  
    Update MultiQC project to add checkatlas parsing. Dev project in: https://github.com/becavin-lab/MultiQC/tree/checkatlas

https://checkatlas.readthedocs.io/en/stable/

## Examples

1. Evaluate and compare different atlases: https://github.com/becavin-lab/checkatlas/blob/3a4f88e94716c09a3b9c86010f570743a5855461/examples/Atlas_comparison.ipynb

https://checkatlas.readthedocs.io/en/stable/CheckAtlas_example_1/CheckAtlas_example_1.html

2. Evaluate different version of your atlas: https://github.com/becavin-lab/checkatlas/blob/3a4f88e94716c09a3b9c86010f570743a5855461/examples/Version_comparison.ipynb

https://checkatlas.readthedocs.io/en/stable/CheckAtlas_example_2/CheckAtlas_example_2.html

3. Explore Scanpy, Seurat and CellRanger objects in your folder: https://github.com/becavin-lab/checkatlas/blob/main/examples/AtlasType_comparison.ipynb

https://checkatlas.readthedocs.io/en/stable/CheckAtlas_example_3/CheckAtlas_example_3.html

## Installation

CheckAtlas can be downloaded from PyPI. However, the project is in an early development phase. We strongly recommend to use the developmental version.

### Install checkatlas development version

```bash
git clone git@github.com:becavin-lab/checkatlas.git
pip install checkatlas/.
```

Install MultiQC with checkatlas file management. This version of MultiQC is available at checkatlas branch of github.com:becavin-lab/MultiQC.

```bash
git clone git@github.com:becavin-lab/MultiQC.git
cd MultiQC/
git checkout checkatlas
pip install .
```

### Install it from PyPI

```bash
pip install checkatlas
```

## Usage

The one liner way to run checkatlas is the following: 

```bash
$ cd your_search_folder/
$ python -m checkatlas .
#or
$ checkatlas .
```

Or run it inside your python workflow.

```py
from checkatlas import checkatlas
checkatlas.run(path, atlas_list, multithread, n_cpus)
```


## Development

Read the [CONTRIBUTING.md](docs/contributing.md) file.

Project developed thanks to the project template : (https://github.com/rochacbruno/python-project-template/)

