'''
This script will generate a matrix of pseudo image features from expression matrix.
A conda env is available in the project folder. 

# input 1 : counts matrix, raw counts, sample by gene name
# input 2 : spot metadata, MUST have 'row','col' columns.
# input 3 : output path, optional, default is ./intermediate
# input 4 : clustering algorithm, optional. {'hdbscan','dp'}. Defaul 'hdbscan'
# input 5 : input dataset type, optional. {'slideseq','visium'}. Defaul 'slideseq'

# output 1: pseudo_image_features.csv
                pseudo image feature matrix csv, sample by gene
# output 2: dimention_reduced_data_for_clustering.csv
                dimension reduced data, intermediate data needed for dp.r
# output 3: pseudo image tif, for visulization only
# output 4: Clustering plot, for debug purpose only.

# example usage:
module load R
conda activate /project/shared/xiao_wang/projects/MOCCA/conda_env
python /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/script/pseudo_image_gen.py \
    /project/shared/xiao_wang/projects/MOCCA/data/Visium/LN/Counts.txt \
    /project/shared/xiao_wang/projects/MOCCA/data/Visium/LN/Spot_metadata.csv \
    dp visium
   
'''
import pandas as pd
import numpy as np
import os
import sys
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale
# from skimage import io
# import seaborn as sns

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

def make_pseudo_img(
    cts_fn, spot_metadata_fn, output_path, algo, 
    dp_script_path = '/home2/s190548/work_xiao/projects/MOCCA/code/Spatial_denoise/script/dirichlet_process_clustering.R'):
    if not os.path.exists(output_path + '/dimention_reduced_data_for_clustering.csv'):
        from umap import UMAP
        # Reads data, and align two matrices.
        print('Loading counts data and performing dimension reduction using UMAP.')
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
        cts = cts.dropna() # some visium data has blank spots

        # Identify highly variable genes.
        # Todo: Evaluate if this cutoff makes sense. Or make it a parameter.
        if cts.shape[1] > 3000: # targeted sequencing libarary
            means, disp = cal_norm_dispersion(cts)
            hvg = cts.columns[(means>0.01) & (disp>0.01)]
        else:
            hvg = cts.columns
        cts_hvg = cts[hvg]
        cts_hvg.loc[:,:] = scale(np.log2(cts_hvg+1))
        spot_meta = spot_meta.loc[cts_hvg.index]


        # Dimention reduction using UMAP on top 25 PCs. 15 Umap components are used.
        # Todo: Evaluate if these works well in practice. A smaller number of PCs will tend to capture more distinct clusters.
        clustering_data = UMAP(
            n_neighbors = 15,n_components=5, min_dist=0.1,spread=1,random_state=0
            ).fit_transform(PCA(25).fit_transform(cts_hvg))
        clustering_data = pd.DataFrame(clustering_data, index = cts.index)
        clustering_data.to_csv(output_path + '/dimention_reduced_data_for_clustering.csv')
    else: # If the dimention reduced data exists, use it.
        print('Found existing dimension reduced data.')
        clustering_data = pd.read_csv(
            output_path + '/dimention_reduced_data_for_clustering.csv', index_col=0)
        spot_meta = pd.read_csv(spot_metadata_fn, index_col=0).loc[clustering_data.index]
        
    # Clustering the expression data.
    if algo == 'hdbscan':
        import hdbscan
        clusterer = hdbscan.HDBSCAN(
            min_cluster_size = 25, min_samples = 30, prediction_data=True
            ).fit(clustering_data)
        soft_clusters = hdbscan.all_points_membership_vectors(clusterer)
    elif algo == 'dp':
        if os.path.exists(output_path + '/dp_cluster_prob.csv'):
            print('Found existing soft clusters, skipping DP clustering.')
        else:
            # hardcoded dp script location.
            os.system(
                'Rscript {} {} {}'.format(
                    dp_script_path,
                    output_path + '/dimention_reduced_data_for_clustering.csv',
                    output_path
                )
            )
        soft_clusters = pd.read_csv(output_path + '/dp_cluster_prob.csv', index_col=0).values

    # if input_type == 'visium': 
    #     # Generating pseudo image for visulizing purpose.
    #     top3 = np.argsort(soft_clusters.sum(axis=0))[::-1][:3] # Select only the top three clusters.
    #     soft_clusters_top3 = soft_clusters[:,top3]
    #     pseudo_img = np.zeros(
    #         [spot_meta.Row.max() + 1,spot_meta.Col.max() + 1, 3], dtype='uint8')
    #     pseudo_img_clusters = np.zeros(pseudo_img.shape, dtype='uint8')
    #     spot_r = spot_meta.Spot_radius[0]
    #     for i in range(spot_meta.shape[0]):
    #         row, col = spot_meta.iloc[i][['Row','Col']].astype(int)
    #         pseudo_img_clusters[row,col,:] = 255 * soft_clusters_top3[i,:]
    #     # Scale the pseudo image so that each channel range from 0-255
    #     pseudo_img_clusters = (
    #         255 * minmax_scale(pseudo_img_clusters.reshape(-1,3)).reshape(pseudo_img.shape)
    #         ).astype('uint8')
    #     io.imsave(output_path + '/Pseudo_image_cluster_probabilities.tif',pseudo_img_clusters)
    # elif input_type == 'slideseq':
    #     _ = plt.figure(figsize=(10,10))
    #     sns.scatterplot(
    #         x = spot_meta['X'], y =  spot_meta['Y'], s=2,
    #         hue = np.argmax(soft_clusters, axis=1).astype(str))
    #     plt.savefig(output_path + '/Pseudo_image_cluster_probabilities.pdf')
    #     plt.close()

    # Output soft cluster probabilities as the pseudo image features.
    pd.DataFrame(
        soft_clusters, index = clustering_data.index).to_csv(output_path + '/pseudo_image_features.csv')

    ##########
    # Debug purpose codes, output clusters in UMAP 2D space.
    # _ = plt.figure(figsize=(6,6))
    # cluster_labels = np.argmax(soft_clusters, axis=1)
    # sns.scatterplot(
    #     x = clustering_data.iloc[:,0], y = clustering_data.iloc[:,1], s=5,
    #     hue = cluster_labels)
    # plt.savefig(output_path + '/UMAP_clustering.pdf')
    # plt.close()

if __name__ == '__main__':
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
    print("Processing {}".format(cts_fn[:-10]))
    make_pseudo_img(cts_fn, spot_metadata_fn, output_path, algo)