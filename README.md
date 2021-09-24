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
### Important details in the Sprod workflow
#### Feature extraction
##### Dataset with a matching image
Sprod works best when the spatial dataset contains a matching image. A example of such dataset is from the 10X Genomic Visium platform, where a tissue slide is imaged and then divided into thousands of sequencing spots, each barcoded uniquely so that every spot can be mapped back uniquely to the image. For this type of dataset, Sprod will extract image features using the `feature_extraction.extract_img_features` function, which will look at image patches around each spot and extract spot-wise features based on both raw pixel values (intensity features) and on co-occurrence matricies (texture features). The shape of the image patch can be specified using the `feature_mask_shape` keyword.

Note: for block mask shape, the `Row` and `Col` columns must be present in the `Spot_metadata.csv` file.

##### Dataset without a matching image
Sometimes the spatial expression dataset does not have a matching image, such as those from the Slide-Seq platform and a Visium dataset without high-resolution image. In this case, Sprod will apply soft clustering on the spots and use the cluster probabilities as the input features for denoising, which we call pseudo-image features. Spots with similar overall molecular phenotype will have similar features, which is a quality also shared by the features derived from real images. Sprod does this through calling the 
`pseudo_image_gen.make_pseudo_img` function. Currently, the soft clustering is done using either the dirichlet process clustering (https://cran.r-project.org/web/packages/dirichletprocess/index.html) or the HDBSCAN algorithm (https://github.com/scikit-learn-contrib/hdbscan). 

Note: The HDBSCAN package is not installed by default as we found it causes problems in installation sometimes.

#### Handling very big spatial dataset
Sprod works well with datasets of thousands of sequencing spots. However, for large datasets with tens of thousands of spots, special operations must be performed so that Sprod can run smoothly. Srpod employs a splitting-stitching scheme to facilitate large dataset processing. Each Slide-seq dataset is randomly (not dependent on spatial location) divided into n (10 by default) equal-sized subsets, and this process is repeated b (10 by default) times. Then, Sprod denoising is performed on each of the n * b subsets and the denoised results are concatenated. Each spot is exactly denoised b times, and the concatenated denoised data from the n sampling batches are averaged so that the randomness resulting from the sub-sampling are averaged out. This is done using the `slideseq_make_patches.subsample_patches` and the `slide_seq_stiching.stiching_subsampled_patches` functions.

### Example applications
#### Application on Visium
In this example, a public [Visium dataset on ovarian cancer](https://www.10xgenomics.com/resources/datasets/human-ovarian-cancer-whole-transcriptome-analysis-stains-dapi-anti-pan-ck-anti-cd-45-1-standard-1-2-0)  from 10X Genomics is used. As the matching image is a immunofluorescence image with a CD45 channel, it is possible to directly compare the correlation between the CD45 signal with PTPRC expression, which encodes the CD45 protein. 

<img src="https://github.com/yunguan-wang/SPROD/blob/master/img/visium_raw_scatter.j.JPG" height="400" width="400">

As shown in the figure, many spots with high CD45 expression were quite low in PTPRC, suggesting the presence of noise. We applied Sprod to this dataset and then visualized the overlap between CD45 and PTPRC expression, and the overlap was much better after Sprod denoising.

<img src="https://github.com/yunguan-wang/SPROD/blob/master/img/visium_overlap.JPG" height="400" width="800">

In the next example, we evaluated the gene expression dropout level in common expression data formats including bulk RNA-seq, Single-cell RNA-seq, Visium and Slide-Seq, and applied Sprod denoising on the Slide-Seq data. Sprod improved the quality of Slide-Seq data drastically. 

<img src="https://github.com/yunguan-wang/SPROD/blob/master/img/slideseq_dropout_comp.JPG" height="400" width="800">
