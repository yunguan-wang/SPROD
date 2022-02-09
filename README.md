# Sprod
Sprod: De-noising Spatial Transcriptomics Data Based on Position and Image Information

## Introduction
Spatial Transcriptomics (ST) techniques provide gene expression close to or even superior to single cell resolution, while retaining the physical locations of sequencing and sometimes also provide matched pathological images. However, the expression data captured by ST techniques suffer from severe inaccuracies, including but not limited to drop-outs as in regular single cell RNA-sequencing (scRNA-seq) data. To reduce the level of noise in the data, we developped the Sprod tool, which incorporated image information and spot/bead positional information and used latent graph modeling to correct the gene expression data of each spot/bead.

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
cd [path/to/sprod]
pip install .
```
After this, install R >= 4.0.2 and the requred R packages. （If you don't already have those in your system.）

The total installation time is around 10 mintunes. If error occuors please upgrade pip and try again.

## Test installation
We have included a simple testing script `test_examples.py` to test the environment and installation. It is based on the toy example dataset included in this repo. Please note that the example data in this repo is only for the purpose of testing installation. 

Three tests were included in this testing script covering different use cases. 

`single_with_img` is for use cases where a small dataset with matching tif image is present. The data is processed in one run.

`batch_with_img` is for use cases where a big dataset with matching tif image is present. The data is processed in small subsampled patches.

`single_pseudoimage` is for use cases where matching image is not avalable. Cluster probablities of each spot is used as features for sprod.

After installation, run the following code in and the should conclude in a few minutes.

```
python [path/to/sprod/]test_examples.py
```
If everything goes well, you should see outputs similar to the following:
```
Removing previous testing results...Done!
Testing Sprod in single mode
Sprod test single mode succeeded!
...
Testing Sprod in patches mode
Sprod test patches mode succeeded!
...
Testing Sprod in single, pseudoimage mode
Sprod test single pseudoimage mode succeeded!
```

## Usage

### Quick start
Once the input data has been processed into the supported format, the full sprod work flow can be run by calling the `sprod.py` script. Sprod will first try to locate a single `.tif` image file in the input path. If one is found, sprod will assume it is the matching image and extract features from the image. If sprod cannot find an image, it will perform soft clustering on the spots and use the clustering probabilities as the input features for the denoising model. 

```
python [path/to/sprod.py] [path/to/input_folder] [path/to/output_folder]
```

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

Mandantory columns: `X`, `Y`, and `Spot_radius` stand for X coordinates, Y coordinates and the radius of each sequencing spot.

We have included a `data_preprocessing.py` script in this repo for processing raw data in Visium or Slide-seq V2 format. For data from other sources, please refer to the script and process your data into the supported format.

### Feature extraction
Feature extraction (with matching image) or generation (without matching image) is wrapped up in the `sprod.py` script. The process is summerized briefly in below.

#### Dataset with a matching image
Sprod works best when the spatial dataset contains a matching image (such as 10X Visium). For this type of datasets, Sprod will extract image features using the [extract_img_features](https://github.com/yunguan-wang/SPROD/blob/master/sprod/feature_extraction.py#L29) function, which will look at image patches around each spot and extract intensity and texture features. The shape of the image patch can be specified using the `--feature_mask_shape` parameter. For Dataset in which sequencing spots, the region of interest can be the exact sequencing spot ('spot'), or a rectangular box surrounding each spot ('block). 

Note: for block mask shape, the `Row` and `Col` columns must be present in the `Spot_metadata.csv` file.

#### Dataset without a matching image
Sometimes the ST dataset does not have a matching image, such as those from the Slide-Seq platform and a Visium dataset without high-resolution image. In this case, Sprod will apply soft clustering on the spots based on gene expression and will use the cluster identities/probabilities as the input features for denoising, which we call pseudo-image features. Sprod does this through calling the [make_pseudo_img](https://github.com/yunguan-wang/SPROD/blob/master/sprod/pseudo_image_gen.py#L71) function.

### Handling very big spatial dataset
Sprod works well with datasets of thousands of sequencing spots/beads. However, for large datasets with tens of thousands of spots/beads, special operations must be performed so that Sprod can run smoothly. Srpod employs a splitting-stitching scheme to facilitate large dataset processing. Each Slide-seq dataset is randomly (not dependent on spatial location) divided into n (10 by default) equal-sized subsets, and this process is repeated b (10 by default) times. Then, Sprod denoising is performed on each of the n * b subsets and the denoised results are concatenated. Each spot is exactly denoised b times, and the concatenated denoised data from the n sampling batches are averaged so that the randomness resulting from the sub-sampling procedure is averaged out. 


### List of Parameters
```
positional arguments:

  input_path            Input folder containing all necessary files.
  
  output_path           Output path

optional arguments:
  --help, -h            
                        show this help message and exit.
  --input_type, -y      
                        The input type decides the running mode for sprod, select from {'single','patches'}. (default: single)
  --output_prefix, -p   
                        Output prefix used in the output. (default: sprod)
  --sprod_npc, -sn      
                        Number of PCs to use, positive integers. -1 to use all PCs from the features.(default: -1)
  --sprod_umap, -su     
                        Toggle to use UMAP on top of PCA to represent features. (default: False)
  --sprod_R, -r         
                        Spot neighborhood radius ratio, 0-1, radius=R*min(xmax-xmin,ymax-ymin). (default: 0.08)
  --spord_perplexity, -u
                        Perplexity, used in Sprod to find the proper Gaussian kernal for distant representation. (default: 250)
  --sprod_margin, -g    
                        Margin for bisection search, used in Sprod to find the proper Gaussian kernal. smaller => slower => accurate (default: 0.001)
  --sprod_latent_dim, -k
                        Dimension of the latent space used in sprod to represent spots. (default: 10)
  --sprod_graph_reg, -l
                        Regularization term for spot graph contructed in sprod. (default: 1)
  --sprod_weight_reg, -w
                        Regularization term for the denoising weights. (default: 0.625)
  --sprod_diag, -d      
                        Toggle to force graph weights to be diagnoal, useful in reducing over smoothing (default: False)
  --image_feature_type, -i
                        Type of feature extracted. combination from {'spot', 'block'} and {'intensity', 'texture'} with '_' as
                        delimiter. Only relevant if the input dataset contains an matching tif image. (default: spot_intensity)
  --warm_start, -ws     
                        Toggle for warm start, which will skip all preprocessing steps including feature extraction and patch-making. (default: False)
  --num_of_patches, -pn
                        Number of subsampled patches. Only works when --input_type is patches. (default: 10)
  --num_of_batches, -pb
                        How many times subsampling is ran. Only works when --input_type is patches. (default: 10)
  -ci, --img_type
                        Input image type. {'he', 'if'}. The 'if' mode is only tested on Visium-assocaited data. (default: he)
```

A few additional notes for the parameters:

`--type` or `-y` : For small datasets with a few thousands spots, this should be set to `single`, while for larger datasets with tens of thousands of spots, this should be set to `patches` to let Sprod run on subsampled patches in parallell to avoid memory problems.

`--warm_start` or `-ws` : This is helpful if image features were extracted in a previous run and a new run is desired with new parameter sets only for the denoising steps. Sprod will automatically look for the needed files in the input directory. 

## Example applications
### Application on Visium
Expression drop-outs are one important source of the noises that we tackle. We evaluated the gene expression dropout levels of data from bulk RNA-seq, Single-cell RNA-seq, Visium and Slide-Seq, and the Sprod-denoised Slide-Seq data. Sprod improved the quality of Slide-Seq data drastically, in terms of drop-outs. The X-axis is the average RNA expression levels of each gene profiled by each technique/dataset, and the Y axis shows the percentages of counts of exactly 0 for each gene.

<img src="https://github.com/yunguan-wang/SPROD/blob/master/img/slideseq_dropout_comp.JPG" height="250" width="600">

In this following example (shown in our paper as well), a public [Visium dataset on ovarian cancer](https://www.10xgenomics.com/resources/datasets/human-ovarian-cancer-whole-transcriptome-analysis-stains-dapi-anti-pan-ck-anti-cd-45-1-standard-1-2-0) from 10X Genomics is used. As the matching image is an immunofluorescence image with a CD45 channel, it is possible to directly compare the correlation between the CD45 signal with raw/denoised PTPRC expression, which encodes the CD45 protein (the CD45 channel wasn't used during Sprod's denoising). 

<img src="https://github.com/yunguan-wang/SPROD/blob/master/img/visium_raw_scatter.j.JPG" height="300" width="300">

As shown in the figure, many spots with high CD45 staining were quite low in PTPRC expression, suggesting the presence of noise. We applied Sprod to this dataset and then visualized the overlap between CD45 and raw/denoised PTPRC expression, and the overlap was much better after Sprod denoising.

<img src="https://github.com/yunguan-wang/SPROD/blob/master/img/visium_overlap.JPG" height="300" width="600">

For more information regarding the rationale of our study and the performance/usage of Sprod, please refer to our paper (link pending)
