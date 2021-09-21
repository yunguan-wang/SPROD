# SPROD
SPROD: De-noising of Spatial Expression Profiling Data Based on Latent Graph Learning of in situ Position and Image Data

## Introduction
Spatial expression profiling (SEP) techniques provide gene expression close to or even superior to single cell resolution, while retaining the physical locations of sequencing and sometimes also provide matched pathological images. However, the expression data captured by spatial sequencing techniques suffer from severe inaccuracies, including but not limited to drop-outs as in regular single cell RNA-sequencing (scRNA-seq) data. To reduce the level of noise in the data, we developped the SPROD tool, which incorporated image information and spot positional information and used latent graph modeling to correct the gene expression data of each spot based on its neighboring spots.

## Graphical abstract
<img src="https://github.com/yunguan-wang/SPROD/blob/master/img/model.png" height="500" width="500">

## Requirement
R >= 4.0.2, Rcpp, distances, dplyr.
Python >= 3.7, pandas, scikit-image, numpy, scipy, scikit-learn, tables.

## Installation
We strongly recommend using a python virtual env to manage the installation of this algorithm. To do this simply run:

```
python3 -m venv [path/to/env]
```
Then, download this repo and install it.
```
source [path/to/env]/bin/activate
git clone [repo_path]
cd [path/to/env]
pip install .
```
After this, install R >= 4.0.2 and the requred R packages.

## Test your installation
We have included a simple testing script `test_examples.sh` to test the environment and installation. It is based on the toy example dataset included in this repo. 
