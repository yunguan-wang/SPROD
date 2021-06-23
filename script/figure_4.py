#%%
import sys
sys.path.insert(1,'/home2/s190548/codes/quality_of_life_improvement')
from scrna_utils import *
from collections import defaultdict
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndi
from scipy.interpolate import UnivariateSpline
from scipy.spatial.distance import cdist
#%%
def merge_common(lists):
    neigh = defaultdict(set)
    visited = set()
    for each in lists:
        for item in each:
            neigh[item].update(each)
    def comp(node, neigh = neigh, visited = visited, vis = visited.add):
        nodes = set([node])
        next_node = nodes.pop
        while nodes:
            node = next_node()
            vis(node)
            nodes |= neigh[node] - visited
            yield node
    for node in neigh:
        if node not in visited:
            yield sorted(comp(node))

def find_cluster_spline(
    cluster, spot_meta, x_min = None, x_max = None, num_bins = 75):
    '''
    Fit a spline in a target cluster population
    '''
    target_cells = spot_meta[spot_meta.cluster==cluster].copy()
    if x_min is not None:
        target_cells = target_cells[(target_cells.X>=x_min)]
    if x_max is not None:
        target_cells = target_cells[(target_cells.X<=x_max)]
    x_ori = target_cells['X'].values
    y_ori = target_cells['Y'].values
    x_bins = pd.cut(target_cells.X, num_bins)
    target_cells = target_cells.groupby(x_bins).median().dropna()
    x_train = target_cells['X'].values
    y_train = target_cells['Y'].values
    spl = UnivariateSpline(x_train, y_train,k=3)
    xs = np.arange(min(x_train),max(x_train),1)
    ys = spl(xs)
    # plt.plot(xs, ys, 'g', lw=3)
    # plt.scatter(x_ori,y_ori, c = 'r', s=10)
    return xs, ys, x_ori, y_ori
# %%
np.random.seed(0)
slideseq_counts = pd.read_hdf(
    '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/Counts.hdf',
    index_col=0)
spot_meta = pd.read_csv(
    '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/Spot_metadata.csv',
    index_col=0)
output_path = '/project/shared/xiao_wang/projects/MOCCA/figures/'
# %%
# Processing raw data
slideseq_counts = slideseq_counts.set_index('Unnamed: 0')
adata = sc.AnnData(slideseq_counts)
sc.pp.calculate_qc_metrics(adata, inplace=True)
sc.pp.filter_cells(adata,min_counts = 100)
sc.pp.filter_cells(adata,max_counts= 10000)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_per_cell(adata, 1000)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_disp=0.05, min_mean=0.05)
sc.pl.highly_variable_genes(adata)
adata.raw = adata
sc.pp.scale(adata)
sc.tl.pca(adata, n_comps=100)
sc.pp.neighbors(adata, n_pcs=25)
sc.tl.leiden(adata, resolution=1)
sc.tl.umap(adata)
adata.write_h5ad('/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/inprogress_cia.h5ad')

#%%
adata = sc.read_h5ad(
    '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/inprogress_cia.h5ad'
)
spot_meta = pd.read_csv(
    '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/Spot_metadata.csv',
    index_col=0)
spot_meta = spot_meta.loc[adata.obs_names]
spot_meta['cluster'] = 'c_' + adata.obs['leiden'].astype(str)
spot_meta['X'] = spot_meta['X'].astype(int)
spot_meta['Y'] = spot_meta['Y'].astype(int)
output_path = '/project/shared/xiao_wang/projects/MOCCA/figures/CIA/'
cts = adata.raw.to_adata().to_df()
MTgenes = [x for x in adata.var_names if x[:3] == 'mt-']
rps_genes = [x for x in adata.var_names if x[:3] in ['rps','rpl']]
valid_genes = [x for x in adata.var_names if x not in (MTgenes+rps_genes)]
adata.obs['deg'] = adata.obs.leiden.copy()
adata.obs['deg'] = 'c_' + adata.obs['deg'].astype(str)
sc.tl.rank_genes_groups(adata,'deg', method='wilcoxon', n_genes=-1)

# %%
ca1_genes = list(set(pd.read_csv(
    '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/CA1_genes.tsv',
    index_col=0, sep='\t').iloc[:,0]))
ca1_genes = [x for x in ca1_genes if x in adata.var_names]
ca1_genes_slideseq = sc.get.rank_genes_groups_df(adata, 'c_10')
ca1_genes_slideseq = ca1_genes_slideseq[
    (ca1_genes_slideseq.logfoldchanges>=1) & 
    (ca1_genes_slideseq.pvals_adj<=0.05)].names.values
ca1_genes = [x for x in ca1_genes if x in ca1_genes_slideseq]
ca1_genes = list(set(ca1_genes) | set(ca1_genes_slideseq[:100]))
#%%
_ = plt.figure(figsize=(8,8))
sns.scatterplot(
    x = spot_meta.loc[adata.obs_names,'X'],
    y = spot_meta.loc[adata.obs_names,'Y'],
    hue = adata.raw.to_adata().to_df()[ca1_genes].mean(axis=1), s = 1, linewidth=0)

# %%
_ = plt.figure(figsize=(8,8))
sns.scatterplot(
    x = spot_meta.loc[adata.obs_names,'X'],
    y = spot_meta.loc[adata.obs_names,'Y'],
    hue = adata.obs.leiden=='10', s = 1, linewidth=0)
#%%
#%%
# %%
# %%
xs,ys,x_ori, y_ori = find_cluster_spline('c_10', spot_meta, 1500,4500, num_bins=50)
plt.plot(xs, ys, 'g', lw=3)
plt.scatter(x_ori,y_ori, c = 'r', s=10)
# %%
spline = np.array([xs,ys]).T
# %%
spot_dist = cdist(spot_meta.iloc[:,:2],spline)
spot_spline_pos = np.argmin(spot_dist, axis=1)
spot_dist = np.min(spot_dist, axis=1)
# %%
spot_meta['dist_to_spline'] = spot_dist
spot_meta['spline_y'] = ys[spot_spline_pos]
spot_meta['dist_to_spline'] = spot_meta['dist_to_spline'] * (((spot_meta['spline_y'] - spot_meta['Y'])>0)*-2 + 1)
# %%
plot_meta = spot_meta[
    (spot_meta.dist_to_spline<=400)&(spot_meta.dist_to_spline>=-100)].copy()
plot_meta = plot_meta[(plot_meta.X>=1500)&(plot_meta.X<=4500)]
plot_meta['CA1_layer'] = 'soma'
plot_meta.loc[plot_meta.dist_to_spline<=-32.5, 'CA1_layer'] = 'basal neuropil'
plot_meta.loc[plot_meta.dist_to_spline>=32.5, 'CA1_layer'] = 'proximal neuropil'
# %%
_ = plt.figure(figsize=(8,8))
sns.scatterplot(
    data = plot_meta,
    x = 'X',
    y = 'Y',
    hue = 'CA1_layer', s = 5, linewidth=0)
plt.plot(xs, ys, 'r', lw=3)
# %%
genes = ['Hpca']
c_vec = adata.raw.to_adata().to_df()[genes].mean(axis=1)[plot_meta.index]
plot_meta['ca1_gene_exp'] = c_vec
# %%
sns.regplot(
    data = plot_meta, x='dist_to_spline', y='ca1_gene_exp', order=10,
    scatter_kws={"s": 2, 'color':'grey'})

# %%
denoised_cts = pd.read_hdf(
    '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/denoised_cts.h5df',
    index_col=0)
denoised_adata = denoised_cts.loc[adata.obs_names]
denoised_adata.drop(['a','T'], inplace=True, axis=1)
cols = [
    x[1:] if (x[0] == 'X') & (x[1].isnumeric()) else x for x in denoised_adata.columns]
cm_cols = [x for x in cols if x in adata.var_names]
denoised_adata.columns = cols
denoised_adata = sc.AnnData(denoised_adata[cm_cols])
sc.pp.normalize_per_cell(denoised_adata, 1000)
sc.pp.log1p(denoised_adata)
sc.pp.highly_variable_genes(denoised_adata, min_disp=0.05, min_mean=0.05)
sc.pl.highly_variable_genes(denoised_adata)
denoised_adata.raw = denoised_adata
sc.pp.scale(denoised_adata)
sc.tl.pca(denoised_adata, n_comps=100)
sc.pp.neighbors(denoised_adata, n_pcs=25)
sc.tl.leiden(denoised_adata, resolution=0.5)
#%%
adata_raw = adata.raw.to_adata().to_df()
denoised_raw = denoised_adata.to_df()
# %%
gene = ca1_genes
plot_meta['ca1_gene_exp_raw'] = adata_raw.loc[plot_meta.index, gene].values
plot_meta['ca1_gene_exp_denoised'] = denoised_raw.loc[plot_meta.index, gene].values
# %%
sns.regplot(
    data = plot_meta, x='dist_to_spline', y='ca1_gene_exp_raw', order=5,
    scatter_kws={"s": 2, 'color':'grey'}, label='Raw')
sns.regplot(
    data = plot_meta, x='dist_to_spline', y='ca1_gene_exp_denoised', order=10,
    scatter_kws={"s": 2, 'color':'pink'}, label = 'Denoised')
plt.legend(markerscale=5)
plt.ylim((0,1))
# %%
# %%
_ = plt.figure(figsize=(8,8))
sns.scatterplot(
    x = spot_meta.loc[denoised_adata.obs_names,'X'],
    y = spot_meta.loc[denoised_adata.obs_names,'Y'],
    hue = denoised_adata.obs['leiden'], s = 1, linewidth=0)
# %%
