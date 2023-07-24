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
