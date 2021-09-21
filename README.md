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

## Test installation
We have included a simple testing script `test_examples.sh` to test the environment and installation. It is based on the toy example dataset included in this repo. 

## Usage

### Data preparation
Sprod workflow requires two mandatory files, a `Counts.txt` for gene expression data,

||gene1|gene2|gene3|gene4|gene5|
|-----|-----|-----|-----|-----|-----|
|spot1|1.06|2.14|1.36|0.94|1.52|
|spot2|0.97|2.42|1.43|1.21|1.17|
|spot3|0.76|2.16|1.07|1.46|1.47|
|spot4|0.82|2.01|1.25|1.18|2.13|
|spot5|1.01|2.07|1.27|1.22|2.16|

as well as a `Spot_metadata.csv` for positional information of each sequenced spot. 

||X|Y|Spot_radius|
|-----|-----|-----|-----|
|spot1|1|2|0.5|
|spot2|2|4|0.5|
|spot3|6|1|0.5|
|spot4|8|1|0.5|
|spot5|1|7|0.5|

We have included a `data_preprocessing.py` script in this repo for processing raw data in Visium or Slide-seq V2 format. For data from other sources, please refer to the script and process your data into the supported format.

### Quick start
Once the input data has been processed into the supported format, the full sprod work flow can be ran by calling the `denoise_job.py` script.
```
python [path/to/denoise_job.py] [path/to/input_folder] [path/to/output_folder]
```
A few important parameters are shown below.

`--type` or `-y` : For small data set, this should be set to `single`, while for larger datasets, this can be set to `patches` to let Sprod ran on subsampled patches in parallell.

`--warm_start` or `-ws` : If image features were extracted in a previous run, feature extraction step can be skipped by toggling `-ws`. 

For details on the parameters used in denoising, please call the script with a `-h` tag for usage details. 
### Sprod step-by-step
### Example applications
