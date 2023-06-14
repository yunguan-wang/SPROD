# Sprod
Sprod: De-noising Spatially Resolved Transcriptomics Data Based on Position and Image Information

## Introduction
Spatial Resolved Transcriptomics (SRT) techniques provide gene expression close to or even superior to single cell resolution, while retaining the physical locations of sequencing and sometimes also provide matched pathological images. However, the expression data captured by SRT techniques suffer from severe inaccuracies, including but not limited to drop-outs as in regular single cell RNA-sequencing (scRNA-seq) data. To reduce the level of noise in the data, we developped the Sprod tool, which incorporated image information and spot/bead positional information and used latent graph modeling to impute the gene expression data of each spot/bead.

Sprod is published in Nature Methods: https://www.nature.com/articles/s41592-022-01560-w.
## Graphical abstract
<img src="https://github.com/yunguan-wang/SPROD/blob/master/img/model.png" height="300" width="300">

## Installation
We strongly recommend using conda to manage the installation of all dependencies. To do this, simply run:

```
conda create --name sprod
conda activate sprod
# conda config --add channels conda-forge ##If you do not have this channel added before#
conda install r-base=4.0 r-optparse r-distances r-dplyr r-gtools r-mvtnorm r-ggplot2 r-umap
conda install python=3.7 pandas scikit-image=0.17 numpy scipy scikit-learn umap-learn
```
Then, download this repo and install it.
```
git clone [repo_path]
cd [path/to/sprod]
pip install .
```

The total installation time is around 10 mintunes. If error occuors, please upgrade pip and try again.

#### Note:
In a rare scenario where: 1) Sprod is executed under a HPC or cloud environment; and 2) you are using a shared `R` installation; and 3) you do not have write access to the shared `R` path; and 4) you have never installed any packages using the shared 'R', you will need to run `R` in interactive mode and install the R package 'dirichletprocess' manually. 

## Test installation
We have included a simple testing script `test_examples.py` to test the environment and installation. It is based on the toy example dataset included in this repo. Please note that the example data in this repo is only for the purpose of testing installation. 

Three tests were included in this testing script covering different use cases. 

`single_with_img` is for use cases where a small dataset with matching tif image is present. The data is processed in one run.

`batch_with_img` is for use cases where a big dataset with matching tif image is present. The data is processed in small subsampled runs.

`single_pseudoimage` is for use cases where matching image is not avalable. The pseudoimage approach is taken,

After installation, run the following code and the processing should conclude in a few minutes.

```
python [path/to/sprod/]test_examples.py
```
If everything goes well, you should see outputs similar to the following:
```
Removing previous testing results...Done!
Testing Sprod in single mode
Sprod test single mode succeeded!
...
Testing Sprod in batch mode
Sprod test batch mode succeeded!
...
Testing Sprod in single, pseudoimage mode
Sprod test single pseudoimage mode succeeded!
```

## Usage

### Quick start
Once the input data have been processed into the supported format, the full sprod workflow can be run by calling the `sprod.py` script. Sprod will first try to locate a single `.tif` image file in the input path. If one is found, sprod will assume it is the matching image and extract features from the image. If sprod cannot find an image, it will perform soft clustering on the spots and use the clustering probabilities as the input features for the denoising model. 

```
python [path/to/sprod.py] [path/to/input_folder] [path/to/output_folder]
```

### Data preparation
Sprod workflow requires two mandatory files, a `Counts.txt` (with "\t" as the delimiter) for gene expression data,

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

Mandantory columns: 

`X`, `Y`, which stand for the X coordinates and Y coordinates of the center of each sequencing spot。

`Spot_radius` is required when the matching image is offered, which stands for the radius of each sequencing spot and is used in feature extraction; The unit of coordinates must be the same between the matching image and the metadata. 

Optional columns: "Z" for Z coordinates if your spatial information has three dimensions. 

We have included a `data_preprocessing.py` script in this repo for processing raw data in Visium or Slide-seq V2 format. For data from other sources, please refer to the script and process your data into the supported format.

### Feature extraction
Feature extraction (with matching image) or generation (without matching image) is wrapped up in the `sprod.py` script. The feature extraction process is carried out automatically on the fly, and the details are also summerized briefly below.

#### Dataset with a matching image
Sprod works best when the spatial dataset contains a matching pathological image (such as 10X Visium). For this type of datasets, Sprod will extract image features using the [extract_img_features](https://github.com/yunguan-wang/SPROD/blob/master/sprod/feature_extraction.py#L29) function, which will look at image regions around each spot and extract intensity and texture features. The shape of the image region can be specified using the `--feature_mask_shape` parameter. The region of interest can be the exact sequencing spot ('spot'), or a rectangular box surrounding each spot ('block). 

Note: for block mask shape, the `Row` and `Col` columns must be present in the `Spot_metadata.csv` file. These two columns correspond to the row or column indices (starting from 0) of the sequencing spots, and are used to determine the size of the bounding box around each spot.  

#### Dataset without a matching image
Sometimes the SRT dataset does not have a matching image, such as those from the Slide-Seq platform and a Visium dataset without high-resolution image. In this case, Sprod will apply clustering on the spots based on gene expression and will use the cluster identities/probabilities as the input features for denoising, which we call pseudo-image features. Sprod does this through calling the [make_pseudo_img](https://github.com/yunguan-wang/SPROD/blob/master/sprod/pseudo_image_gen.py#L71) function.

### Handling very big spatial dataset
Sprod works well with datasets of thousands of sequencing spots/beads. However, for large datasets with tens of thousands of spots/beads, special operations must be performed so that Sprod can run smoothly. Srpod employs a splitting-stitching scheme to facilitate large dataset processing. Each dataset is randomly (not dependent on spatial location) divided into n (10 by default) equal-sized subsets, and this process is repeated b (10 by default) times. Then, Sprod denoising is performed on each of the n * b subsets and the denoised results are concatenated. Each spot is exactly denoised b times, and the concatenated denoised data from the n sampling batches are averaged so that the randomness resulting from the sub-sampling procedure is averaged out. In the `batch` mode, the computational time is linearly proportional to the total number of sequncing spots/cells. 

<img src = "https://github.com/yunguan-wang/SPROD/blob/QBRC_SS/img/TimexSpots.png" height="270" width="300">

### List of Parameters
```
positional arguments:

  input_path            Input folder containing all necessary files.
  
  output_path           Output path

optional arguments:
  --help, -h            
                        show this help message and exit.
  --input_type, -y      
                        The input type decides the running mode for sprod, select from {'single','batch'}. (default: single)
  --output_prefix, -p   
                        Output prefix used in the output. (default: sprod)
  --sprod_npc, -sn      
                        Number of Principal Components (PCs) to use for the image features, positive integers. -1 to use all PCs from the features. The image features are summarized through PCA before being used by sprod. (default: -1)
  --sprod_umap, -su     
                        Toggle to use UMAP on top of PCA to represent features. (default: False)
  --sprod_R, -r         
                        Spot neighborhood radius ratio, 0-1, radius=R*min(xmax-xmin,ymax-ymin). (default: 0.08)
  --spord_perplexity, -u
                        Perplexity, as in tSNE, during the graph building process. (default: 250)
  --sprod_margin, -g    
                        Margin for the bisection method for solving the perpexity equation. smaller => slower => accurate (default: 0.001)
  --sprod_latent_dim, -k
                        Dimension of the latent space used in sprod to represent spots. (default: 10)
  --sprod_graph_reg, -l
                        Regularization term for the similarity graph contruction in sprod. (default: 1)
  --sprod_weight_reg, -w
                        Regularization term for the denoising weights. (default: 0.625)
  --sprod_diag, -d      
                        Toggle to force graph weights to be diagnoal, useful in reducing over smoothing (default: False)
  --image_feature_type, -i
                        Type of feature extracted. combination of {'spot', 'block'} and {'intensity', 'texture'} with '_' as
                        the delimiter. Only relevant if the input dataset contains a matching tif image. (default: spot_intensity)
  --warm_start, -ws     
                        Toggle for warm start, which will skip all preprocessing steps in feature extraction and batch mode preparation. (default: False)
  --num_of_batches, -pb
                        How many times subsampling is run. Only works when --input_type is "batch". (default: 10)
  --img_type, -ci
                        Input image type. {'he', 'if'}. The 'if' mode is only tested on Visium-associated data. (default: he)
  --custom_feature, -cf
                        Option for advanced users. A custom spot by feature csv file can be used together with sprod extracted
                        features. Must have matching spot names and order. The rows names should be spot names and the columns should
                        be feature names. (default: None)
```

A few additional notes for the parameters:

(1) `--type` or `-y` : For small datasets with a few thousands spots, this should be set to `single`, while for larger datasets with tens of thousands of spots, this should be set to `batch` to let Sprod run on subsamples in parallell to avoid memory problems.

(2) `--warm_start` or `-ws` : This is helpful if the image features were extracted in a previous run and a new run is desired with new parameter sets only for the denoising steps. Sprod will automatically look for the needed files in the input directory. In the warm start mode, the input folder must contains either the real image features or the pseudo image features. The real image features should be named as  `[spot/block]_level_[intensity/texture]_features.csv`, and the pseudo image features should be called `pseudo_image_features.csv`. If both sets of features exist, Sprod will only use the real image features.

(3) `--custom_features` or `-cf`: This optional parameter is for advanced users who wish to use their own feature matrix in addition to sprod's extracted features. This is also where users can import results from cell type deconvolution software such as SpatialDWLS (columns are the cell types and their proportions/probabilities/abundances). Please make sure the input feature matrix is of the correct format (a csv file with spot ids as row names and feature names as column names, and the row names should match other inputs on a one-to-one basis).

An example custom feature matrix can be found at `test_example/input/custom_features_deconvoluted_cell_types.csv` 

### Automatic parameter selection
We have included a script for automatic grid search of the tuning parameters, which is `parameter_selection.py`. User-specified paramater combinations are tested and the qualities of the constructed latent graphs are evaluated, by two quantitative metrics (see our paper for details). The the diagnostic plots for each parameter combination will be saved in the output folder. In addition, this script will output a table showing the parameter combinations, rankable by these two metrics and the average of their ranks. An example of this table is shown below (rank refers to the average of rank(graph distance) and rank(image distance)).

||Average Graph Distance|Average Image Distance|Rank|
|-----|-----|-----|-----|
|R-0.1_K-5_L-15|759.78|2.51|7.0|
|R-0.1_K-5_L-10|759.63|2.70|7.5|
|R-0.1_K-3_L-5|754.38|2.82|8.0|
|R-0.1_K-10_L-15|754.94|2.88|9.0|

An example of how to use this script, including both the inputs and outputs of this functionality, can be found in `parameter_selection_example/outputs`

We prefer smaller values of these two metrics and smaller ranks. But we suggest the users to check the actual diagnostic plots of at least the top few parameter combinations, to select the optimal parameter set. 

### Contact Us
If you have any suggestions/ideas for Sprod or are having issues trying to use it, please don't hesitate to reach out to us.

Yunguan Wang, yunguan[dot]wang@utsouthestern[dot]edu

Bing Song, bing[dot]song@utsouthestern[dot]edu

Tao Wang, tao[dot]wang@utsouthestern[dot]edu
