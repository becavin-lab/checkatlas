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

### Install Seurat

To be able to manage seurat file, rpy2 should have Seurat installed. The easiest way is to put all checkatlas requirements in a conda environment and add r-seurat.

```bash
conda create -n checkatlas python=3.9
pip install checkatlas
conda install -c bioconda r-seurat
```

Or, open R in checkatlas environment (the one where you ran 'pip install') and install Seurat.

```bash
% R
> install.packages('Seurat')
> library(Seurat)
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