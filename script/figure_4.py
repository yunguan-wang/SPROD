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
import gseapy as gs
import scipy.cluster.hierarchy as sch
import os
from sklearn.cluster import SpectralCoclustering
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

# %%
np.random.seed(0)
slideseq_counts = pd.read_hdf(
    '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/Counts.hdf',
    index_col=0)
spot_meta = pd.read_csv(
    '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/Spot_metadata.csv',
    index_col=0)
denoised_cts = pd.read_hdf(
    '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/denoised_cts.h5df',
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
sc.tl.leiden(adata, resolution=0.5)
sc.tl.umap(adata)
adata.write_h5ad('/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/inprogress_cia.h5ad')
# %%
denoised_adata = denoised_cts.loc[adata.obs_names]
denoised_adata = sc.AnnData(denoised_adata)
sc.pp.normalize_per_cell(denoised_adata, 1000)
#%%
adata = sc.read_h5ad(
    '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/inprogress_cia.h5ad'
)
spot_meta = pd.read_csv(
    '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/Spot_metadata.csv',
    index_col=0)
spot_meta = spot_meta.loc[adata.obs_names]
spot_meta['cluster'] = 'c_' + adata.obs['leiden'].astype(str)
output_path = '/project/shared/xiao_wang/projects/MOCCA/figures/CIA/'
cts = adata.raw.to_adata().to_df()
MTgenes = [x for x in adata.var_names if x[:3] == 'mt-']
rps_genes = [x for x in adata.var_names if x[:3] in ['rps','rpl']]
valid_genes = [x for x in adata.var_names if x not in (MTgenes+rps_genes)]
adata.obs['deg'] = adata.obs.leiden.copy()
adata.obs['deg'] = 'c_' + adata.obs['deg'].astype(str)
sc.tl.rank_genes_groups(adata,'deg', method='wilcoxon', n_genes=-1)
#%%
_ = plt.figure(figsize=(8,8))
sns.scatterplot(
    x = spot_meta.loc[adata.obs_names,'X'],
    y = spot_meta.loc[adata.obs_names,'Y'],
    hue = adata.raw.to_adata().to_df()['Arhgap12'], s = 1, linewidth=0)
# %%
_ = plt.figure(figsize=(8,8))
sns.scatterplot(
    x = spot_meta.loc[adata.obs_names,'X'],
    y = spot_meta.loc[adata.obs_names,'Y'],
    hue = denoised_adata.to_df()['Arhgap12'], s = 1, linewidth=0)
# %%
