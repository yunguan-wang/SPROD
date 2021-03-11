'''
This script will generate a matrix of pseudo image features from expression matrix.
A conda env is available in the project folder. 

# input 1 : counts matrix, raw counts, sample by gene name
# input 2 : spot metadata, MUST have 'row','col', 'X', 'Y', and 'Spor_radius' columns.
# input 3 : output path, optional, default is ./intermediate
# input 4 : clustering algorithm, optional. {'hdbscan','dp'}. Defaul 'hdbscan'

# output 1: pseudo_image_features.csv
                pseudo image feature matrix csv, sample by gene
# output 2: dimention_reduced_data_for_clustering.csv
                dimension reduced data, intermediate data needed for dp.r
# output 3: pseudo image tif, for visulization only
# output 4: Clustering plot, for debug purpose only.

# example usage:
conda activate /project/shared/xiao_wang/projects/MOCCA/conda_env
python /project/shared/xiao_wang/projects/MOCCA/code/pseudo_image_gen.py \
    /project/shared/xiao_wang/projects/MOCCA/data/Visium/LN/Counts.txt \
    /project/shared/xiao_wang/projects/MOCCA/data/Visium/LN/Spot_metadata.csv
   
'''
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import sys
from skimage import io, morphology
from sklearn.decomposition import PCA
from sklearn.preprocessing import minmax_scale, scale
from umap import UMAP
import hdbscan
import seaborn as sns

def cal_norm_dispersion(cts):
    '''
    Adapted from Scanpy _highly_variable_genes_single_batch.
    https://github.com/theislab/scanpy/blob/f7279f6342f1e4a340bae2a8d345c1c43b2097bb/scanpy/preprocessing/_highly_variable_genes.py
    '''
    mean, var = cts.mean(),  cts.std()
    mean[mean == 0] = 1e-12  # set entries equal to zero to small value
    dispersion = var / mean
    dispersion[dispersion == 0] = np.nan
    dispersion = np.log(dispersion)
    mean = np.log1p(mean)

    # all of the following quantities are "per-gene" here
    df = pd.DataFrame()
    df['means'] = mean
    df['dispersions'] = dispersion
    df['mean_bin'] = pd.cut(df['means'], bins=20)
    disp_grouped = df.groupby('mean_bin')['dispersions']
    disp_mean_bin = disp_grouped.mean()
    disp_std_bin = disp_grouped.std(ddof=1)
    # retrieve those genes that have nan std, these are the ones where
    # only a single gene fell in the bin and implicitly set them to have
    # a normalized disperion of 1
    one_gene_per_bin = disp_std_bin.isnull()
    disp_std_bin[one_gene_per_bin.values] = disp_mean_bin[
        one_gene_per_bin.values
    ].values
    disp_mean_bin[one_gene_per_bin.values] = 0
    # actually do the normalization
    df['dispersions_norm'] = (
        df['dispersions'].values  # use values here as index differs
        - disp_mean_bin[df['mean_bin'].values].values
    ) / disp_std_bin[df['mean_bin'].values].values
    return df.means.values, df.dispersions_norm.values

# Parsing arguments.
cts_fn = sys.argv[1]
spot_metadata_fn = sys.argv[2]
try:
    output_path = sys.argv[3]
except IndexError:
    output_path = './intermediate'

try:
    algo = sys.argv[4]
except IndexError:
    algo = 'hdbscan'
if not os.path.exists(output_path):
    os.makedirs(output_path)

# Dimension reduction, can be skipped if it is already done.
if not os.path.exists(output_path + '/dimention_reduced_data_for_clustering.csv'):
    # Reads data, and align two matrices.
    cts = pd.read_csv(cts_fn,sep='\t',index_col=0)
    spot_meta = pd.read_csv(spot_metadata_fn, index_col=0)
    cm_samples = [x for x in cts.index if x in spot_meta.index]
    if len(cm_samples) != cts.shape[0]:
        print('Counts matrix and spot metadata have unequal number of samples!')
    cts = cts.loc[cm_samples]
    spot_meta = spot_meta.loc[cm_samples]

    # Filterout genes expressed in less than 10 spots.
    # Todo: Evaluate if this cutoff makes sense in practice.
    valid_genes = cts.columns[((cts > 0).sum() >= 10)]
    # CPM normalization
    cts = 1e4 * cts.apply(lambda x: x/x.sum(), axis=1)
    cts = cts[valid_genes]

    # Identify highly variable genes.
    # Todo: Evaluate if this cutoff makes sense. Or make it a parameter.
    means, disp = cal_norm_dispersion(cts)
    hvg = cts.columns[(means>0.01) & (disp>0.01)]
    cts_hvg = np.log2(cts[hvg]+1)
    cts_hvg = cts_hvg.reindex(spot_meta.index)

    # Dimention reduction using UMAP on top 25 PCs. 5 Umap components are used.
    # Todo: Evaluate if these works well in practice. A smaller number of PCs will tend to capture more distinct clusters.
    clustering_data = UMAP(
        n_neighbors = 15,n_components=5, min_dist=0.1,spread=1,random_state=0
        ).fit_transform(PCA(25).fit_transform(scale(cts_hvg)))
    pd.DataFrame(
        clustering_data, index = cts.index
        ).to_csv(output_path + '/dimention_reduced_data_for_clustering.csv')
else: # If the dimention reduced data exists, use it.
    clustering_data = pd.read_csv(
        output_path + '/dimention_reduced_data_for_clustering.csv', index_col=0)
    spot_meta = pd.read_csv(spot_metadata_fn, index_col=0)
    
# Clustering the expression data.
if algo == 'hdbscan':
    clusterer = hdbscan.HDBSCAN(
        min_cluster_size = 25, min_samples = 30, prediction_data=True
        ).fit(clustering_data)
    soft_clusters = hdbscan.all_points_membership_vectors(clusterer)
elif algo == 'DP':
    # Todo: call DP.r or implement a pure python version.
    pass
    # os.system('Rscript dp.r umap_comps.csv')
    # soft_clusters = pd.read_csv('dp_clusters.csv', index_col=0)

# Generating pseudo image for visulizing purpose.
top3 = np.argsort(soft_clusters.sum(axis=0))[::-1][:3] # Select only the top three clusters.
soft_clusters_top3 = soft_clusters[:,top3]
pseudo_img = np.zeros(
    [spot_meta.Row.max() + 1,spot_meta.Col.max() + 1, 3], dtype='uint8')
pseudo_img_clusters = np.zeros(pseudo_img.shape, dtype='uint8')
spot_r = spot_meta.Spot_radius[0]
spot_mask = morphology.disk(spot_r)
for i in range(spot_meta.shape[0]):
    row, col, X, Y = spot_meta.iloc[i,1:5].astype(int)
    pseudo_img_clusters[row,col,:] = 255 * soft_clusters_top3[i,:]
# Scale the pseudo image so that each channel range from 0-255
pseudo_img_clusters = (
    255 * minmax_scale(pseudo_img_clusters.reshape(-1,3)).reshape(pseudo_img.shape)
    ).astype('uint8')
io.imsave(output_path + '/Pseudo_image_cluster_probabilities.tif',pseudo_img_clusters)

# Output soft cluster probabilities as the pseudo image features.
pd.DataFrame(
    soft_clusters, index = clustering_data.index).to_csv(output_path + '/pseudo_image_features.csv')

##########
# Debug purpose codes, output clusters in UMAP 2D space.
_ = plt.figure(figsize=(6,6))
sns.scatterplot(
    clustering_data.iloc[:,0], clustering_data.iloc[:,1], s=5, hue = clusterer.labels_)
plt.savefig(output_path + '/UMAP_clustering.pdf')


