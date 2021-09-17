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

def calculate_neighbor_radius(
    spot_meta, sample_size=1000, target_n_neighbors = 10, margin = 10,
    r_min = 1, r_max = 1000):
    r_samples = spot_meta.loc[
        np.random.choice(spot_meta.index,sample_size, False)].copy().loc[:,['X','Y']]
    r_sample_dist = cdist(r_samples, spot_meta.loc[:,['X','Y']])
    n_steps = 1
    while n_steps <= 100:
        r_next = (r_max + r_min)/2
        n_neighbors = np.median(np.sum((r_sample_dist <= r_next),axis=1))
        if n_neighbors == target_n_neighbors:
            break
        elif (r_max - r_min)/2 < margin:
            break 
        n_steps += 1
        nn_r_min = np.median(np.sum((r_sample_dist <= r_min),axis=1))
        nn_r_max = np.median(np.sum((r_sample_dist <= r_max),axis=1))
        if (n_neighbors>target_n_neighbors) == (nn_r_max > target_n_neighbors):
            r_max = r_next
        elif (n_neighbors>target_n_neighbors) == (nn_r_min> target_n_neighbors):
            r_min = r_next
    return r_next

def find_neighbor_clusters(spot_meta, R_max, R_min):
    neighborhood = spot_meta.copy()
    euc_dist = cdist(neighborhood.iloc[:,:2], neighborhood.iloc[:,:2])
    # remove self neighbor
    np.fill_diagonal(euc_dist, np.inf)
    R_criteria = np.logical_and((euc_dist <= R_max),(euc_dist >= R_min))
    rows, cols = np.nonzero(R_criteria)
    neighbors = pd.Series(cols, index=rows, dtype=object)
    neighbors = [neighbors.loc[[x]].tolist() for x in neighbors.index.unique()]
    neighbors = pd.Series(neighbors, index = list(set(rows)))

    neighbors[:] = [
        [neighborhood.index[x]] if not isinstance(x, list) else neighborhood.index[x] for x in neighbors.values]
    # Update spotmedata data with neighborhood information
    neighborhood['cluster'] = 'c_' + adata.obs.loc[neighborhood.index,'leiden'].astype(str)
    neighborhood.loc[neighborhood.index[neighbors.index],'neighbor_idx'] = neighbors.values
    neighborhood = neighborhood.dropna(subset = ['neighbor_idx'])
    neighborhood['neighbor_cluster'] = neighborhood.neighbor_idx.apply(
        lambda x: list(set(neighborhood.loc[x,'cluster']))).values
    neighborhood['num_neighbor_clusters'] = neighborhood['neighbor_cluster'].apply(len)
    return neighborhood

def find_interface_cluster(cluster, spot_meta, linkage_R = 20):
    # Identify the major region in the cluster that constitute the interface.
    cluster_df = spot_meta[spot_meta.cluster==cluster].index.values
    lm = sch.linkage(spot_meta.loc[cluster_df, ['X','Y']])
    clusters = sch.fcluster(lm, linkage_R, 'distance')
    unique_clusters, cluster_counts = np.unique(clusters, return_counts=True)
    cluster_fraction = cluster_counts/cluster_counts.sum()
    large_clusters = np.where(cluster_fraction>=0.05)
    # print(cluster_counts[large_clusters], unique_clusters[large_clusters])
    large_clusters = unique_clusters[large_clusters]
    cluter_area_idx = cluster_df[np.isin(clusters, large_clusters)]
    return cluter_area_idx.tolist()

def find_interface(interface_area_idx, interface_cluster, spot_meta):
    area1 = interface_area_idx[0]
    area2 = interface_area_idx[1]
    cluster1 = interface_cluster[0]
    cluster2 = interface_cluster[1]
    interface_idx = []
    # find area1 interface in area2
    interface_df = spot_meta.loc[area1]
    interface_df = interface_df[interface_df.neighbor_cluster.apply(
        lambda x: cluster2 in x)]
    interface_df = interface_df[
        interface_df.neighbor_idx.apply(lambda x: any([i for i in x if i in area2]))]
    interface_idx += interface_df.index.tolist()
    # find area 2 interface in area1
    interface_df = spot_meta.loc[area2]
    interface_df = interface_df[interface_df.neighbor_cluster.apply(
        lambda x: cluster1 in x)]
    interface_df = interface_df[
        interface_df.neighbor_idx.apply(lambda x: any([i for i in x if i in area1]))]
    interface_idx += interface_df.index.tolist()
    return interface_idx

def distance_pairing(spot_meta, celltypeA, celltypeB):
    cell_dist = cdist(spot_meta.loc[celltypeB,['X','Y']], spot_meta.loc[celltypeA,['X','Y']])
    if cell_dist.shape[0] >= cell_dist.shape[1]:
        rowwise_mins = np.argmin(cell_dist,axis=1)
        cell_pairs = pd.DataFrame(
            np.array(celltypeA, dtype=object)[rowwise_mins], index = celltypeB, columns = ['CelltypeA'])
        cell_pairs.index.name = 'CelltypeB'
    else:
        colwise_mins = np.argmin(cell_dist,axis=0)
        cell_pairs = pd.DataFrame(
            np.array(celltypeB, dtype=object)[colwise_mins], index = celltypeA, columns = ['CelltypeB'])
        cell_pairs.index.name = 'CelltypeA'

    cell_pairs['dist'] = np.diag(cdist(
        spot_meta.loc[cell_pairs.index,['X','Y']],
        spot_meta.loc[cell_pairs.values.flatten(),['X','Y']]))

    return cell_pairs

def plot_interface(interface_idx, spot_meta, output_fn = None):
    non_interface_data = spot_meta.loc[
    [x for x in spot_meta.index if x not in interface_idx]]
    interface_area = spot_meta.loc[interface_idx]
    _ = plt.figure(figsize=(8,8))
    hue_order = sorted(spot_meta.cluster.unique())
    sns.scatterplot(
        data = interface_area, x = 'X', y = 'Y', s=2, hue='cluster', linewidth=0,
        hue_order=hue_order)
    sns.scatterplot(
        data = non_interface_data, x = 'X', y = 'Y', s=1, hue='cluster', alpha=0.3,
        linewidth=0, legend=None,hue_order=hue_order)
    if output_fn is not None:
        plt.savefig(output_fn)
        plt.close()

def find_pairs(
    spot_meta, celltype1, celltype2, min_n = 3, target_n = 20, max_linkage_dist = 40):
    R_min = calculate_neighbor_radius(spot_meta,target_n_neighbors=min_n)
    R_max = calculate_neighbor_radius(spot_meta,target_n_neighbors=target_n)
    neighborhood = find_neighbor_clusters(spot_meta, R_max, R_min)
    c1_area = find_interface_cluster(celltype1, neighborhood, max_linkage_dist)
    c2_area = find_interface_cluster(celltype2, neighborhood, max_linkage_dist)
    c1_c2_interface = find_interface([c1_area, c2_area], [celltype2,celltype1], neighborhood)
    c1_interface = [x for x in c1_c2_interface if x in c1_area]
    non_interface_c1 = [x for x in c1_area if x not in c1_interface]
    c2_interface = [x for x in c1_c2_interface if x in c2_area]
    non_interface_c2 = [x for x in c2_area if x not in c2_interface]
    interface_pairs = distance_pairing(spot_meta, c1_interface, c2_interface)
    noninterface_pairs = distance_pairing(spot_meta, non_interface_c1, non_interface_c2)
    return interface_pairs, noninterface_pairs,c1_c2_interface

def cal_correlations(
    interface_pairs, noninterface_pairs, cts, min_noninterface_dist=500):
    interface_corr = 1 - cdist(
        cts.loc[interface_pairs.index].transpose(),
        cts.loc[interface_pairs.iloc[:,0].values].transpose(),metric='correlation'
    )
    non_interface_cells = noninterface_pairs[noninterface_pairs['dist']>=min_noninterface_dist].index
    non_interface_cell_t = noninterface_pairs.loc[non_interface_cells].iloc[:,0].values
    noninterface_corr = 1 - cdist(
        cts.loc[non_interface_cells].transpose(),
        cts.loc[non_interface_cell_t].transpose(),metric='correlation'
    )
    interface_corr = np.nan_to_num(interface_corr)
    noninterface_corr = np.nan_to_num(noninterface_corr)
    interface_corr = pd.DataFrame(
        interface_corr, index = cts.columns, columns=cts.columns)
    noninterface_corr = pd.DataFrame(
        noninterface_corr, index = cts.columns, columns=cts.columns)
    return interface_corr, noninterface_corr

def find_gene_modules(interface_corr, noninterface_corr, corr_cutoff = 0.8):
    np.random.seed(0)
    gene_dist = pd.DataFrame(interface_corr)
    crit = pd.DataFrame(np.logical_and(
        interface_corr > noninterface_corr, interface_corr>=corr_cutoff))
    # gene_dist = pd.DataFrame(interface_corr - noninterface_corr)
    # crit = pd.DataFrame(gene_dist>=corr_cutoff)
    valid_rows = crit[crit.sum(axis=1)>0].index
    valid_cols = crit.columns[crit.sum(axis=0)>0]
    gene_dist = gene_dist.loc[valid_rows, valid_cols]
    # Using bi-clustering to get rid of gene not forming modules.
    if gene_dist.shape[0] > 1:
        num_clusters = np.min([int(np.sqrt(gene_dist.shape[0])),10])
        cluster_model = SpectralCoclustering(n_clusters=num_clusters, n_jobs=16)
        cluster_model.fit(gene_dist)
        fit_data = gene_dist.iloc[np.argsort(cluster_model.row_labels_)]
        fit_data = fit_data.iloc[:, np.argsort(cluster_model.column_labels_)]
        rows_to_remove = []
        cols_to_remove = []
        next_round_num_clusters = num_clusters
        for gc in range(num_clusters):
            rows = gene_dist.index[cluster_model.rows_[gc]]
            cols = gene_dist.columns[cluster_model.columns_[gc]]
            gene_module_corr = gene_dist.loc[rows, cols].values.flatten().mean()
            if gene_module_corr < (corr_cutoff/5):
                rows_to_remove += rows.tolist()
                cols_to_remove += cols.tolist()
                next_round_num_clusters -= 1
        fit_data.drop(rows_to_remove, inplace = True)
        fit_data.drop(cols_to_remove, axis=1, inplace = True)
        # 2nd round of biclustering
        while True:
            if next_round_num_clusters == 1:
                print('No clusters found')
                break
            cluster_model = SpectralCoclustering(
                n_clusters=next_round_num_clusters, n_jobs=16)
            cluster_model.fit(fit_data)
            row_clusters = cluster_model.row_labels_
            col_clusters = cluster_model.column_labels_
            if sorted(np.unique(row_clusters)) != sorted(np.unique(col_clusters)):
                next_round_num_clusters -= 1
            else:
                break
        # return gene module
        gene_modules = []
        for c in range(next_round_num_clusters):
            row_genes = fit_data.index[row_clusters==c].tolist()
            col_genes = fit_data.columns[col_clusters==c].tolist()
            gene_modules.append(row_genes + col_genes)
        # return gene module heatmap data.
        fit_data = fit_data.iloc[np.argsort(cluster_model.row_labels_)]
        fit_data = fit_data.iloc[:, np.argsort(cluster_model.column_labels_)]
    else:
        gene_modules = [gene_dist.index.tolist() + gene_dist.columns.tolist()]
        fit_data = gene_dist
    return gene_modules, fit_data

def module_enrichment(
    gene_modules,
    libs = [
        'Allen_Brain_Atlas_up','Allen_Brain_Atlas_10x_scRNA_2021',
        'WikiPathways_2019_Mouse']):
    gm_pathways = pd.DataFrame()
    for lib in libs:
        for i, gl in enumerate(gene_modules):     
            res = gs.enrichr(
                list(gl),
                gene_sets=lib, no_plot=True)
            res = res.res2d[res.res2d['Adjusted P-value']<=0.05]
            res['gene_module'] = i
            gm_pathways = gm_pathways.append(res)
    return gm_pathways
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
sc.tl.leiden(adata, resolution=0.5)
sc.tl.umap(adata)
adata.write_h5ad('/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/inprogress_cia.h5ad')
# %%
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
_ = plt.figure(figsize=(6,6))
cluster_labels = adata.obs.leiden.astype(str)
sns.scatterplot(
    x = spot_meta.loc[adata.obs_names,'X'], y = spot_meta.loc[adata.obs_names,'Y'], s=2,
    hue = cluster_labels, linewidth=0)
plt.savefig(
    '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/pseudo_img_leiden.pdf')
plt.close()
#%%
# Differential analysis
for cluster in adata.obs.deg.unique():
    _deg = sc.get.rank_genes_groups_df(adata, cluster)
    _deg = _deg[_deg.names.isin(valid_genes)]
    _deg.to_excel(
        '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/{}_DEG.xlsx'.format(cluster),
        index=None)
#%%
interface = ['c_0','c_3']
# get genes highly expressed in the interface clusters.
adata_hvg = adata.raw.to_adata()
adata_hvg = adata_hvg[spot_meta[spot_meta.cluster.isin(interface)].index]
sc.pp.highly_variable_genes(adata_hvg,min_disp=0.5)
hvg = adata_hvg.var_names[adata_hvg.var.highly_variable]
hvg = [x for x in hvg if x in valid_genes]

interface_pairs,noninterface_pairs, interface_cells = find_pairs(
    spot_meta, interface[0], interface[1])
interface_corr, noninterface_corr = cal_correlations(interface_pairs, noninterface_pairs, cts[hvg])
gene_modules,gm_data = find_gene_modules(interface_corr, noninterface_corr, 0.5)
gm_pathways = module_enrichment(
    gene_modules,
    libs = [
        'Allen_Brain_Atlas_up','GO_Molecular_Function_2018','WikiPathways_2019_Mouse',
        'GO_Cellular_Component_2018','GO_Biological_Process_2018'])
#%%
# outputting results
plot_interface(
    interface_cells, spot_meta,
    output_fn=output_path + '{}_{}_interface.pdf'.format(*interface))
# Interface vs non-interface gene correlations
_ = plt.figure(figsize=(9,9))
_ = plt.hist(
    interface_corr.values.flatten(), bins=100, label='Interface', alpha=0.75)
_ = plt.hist(
    noninterface_corr.values.flatten(), bins=100, label='Non-Interface', alpha=0.75)
plt.legend()
plt.savefig(output_path + '{}_{}_interface_non_interface_gene_correlations.pdf'.format(*interface))
# Gene module heatmap
_ = plt.figure(figsize=(18,18))
ax = sns.heatmap(gm_data, cmap=plt.cm.Blues, xticklabels=True, yticklabels = True)
plt.savefig(output_path + '{}_{}_gene modules.pdf'.format(*interface))
plt.close()
# Save modules
pd.DataFrame(gene_modules).transpose().to_excel(
    output_path +  '{}_{}_interface_gene_modules.xlsx'.format(*interface))
# Pathway heatmap
gm_pathways = gm_pathways.sort_values('Adjusted P-value')
pathway_hmap = gm_pathways.groupby('gene_module').head(10)
pathway_hmap = pathway_hmap.pivot(index = 'Term', columns = 'gene_module', values = 'Adjusted P-value')
pathway_hmap = pathway_hmap.fillna(1)
pathway_hmap = -np.log10(pathway_hmap)
_ = plt.figure(figsize=(12,24))
sns.clustermap(
    pathway_hmap, yticklabels=True)
plt.savefig(
    output_path +  '{}_{}_interface_gene_modules_pathways.pdf'.format(*interface),
    bbox_inches = 'tight')
plt.close()
#%%
# # cols_as_rowmin = pd.Series.value_counts(row_mins)
# multi_cols = cols_as_rowmin[cols_as_rowmin > 1].index.values
# col_interface_cdist = interface_cdist[:,multi_cols]
# # %%
# meta_interface = neighborhood.loc[c3_c0_interface,'X':'cluster']
# #%%
# pca = PCA(n_components=1).fit(meta_interface.iloc[:,:2])
# pc_coord = pca.transform(meta_interface.iloc[:,:2])
# meta_interface['PC_1'] = pc_coord
# meta_interface.loc[c0_interface, 'bin'] = pd.cut(
#     meta_interface.loc[c0_interface, 'PC_1'], 20,
#     labels = ['bin_' + str(x) for x in np.arange(1,21)]).values
# meta_interface.loc[c3_interface, 'bin'] = pd.cut(
#     meta_interface.loc[c3_interface, 'PC_1'], 20,
#     labels = ['bin_' + str(x) for x in np.arange(1,21)]).values

# raw_cts = adata.raw.to_adata().to_df()
# # Calculate interface score for the interface
# c0_bined = raw_cts.loc[c0_interface].groupby(
#     meta_interface.loc[c0_interface, 'bin']).mean()
# c3_bined = raw_cts.loc[c3_interface].groupby(
#     meta_interface.loc[c3_interface, 'bin']).mean()
# c0_bined = c0_bined.sort_index()
# c3_bined = c3_bined.sort_index()
# c0_c3_intercept = c0_bined.clip(upper = c3_bined, axis=None)
# gene_interface_score = 2 * c0_c3_intercept.sum() / (c0_bined.sum() + c3_bined.sum())

# meta_non_interface = neighborhood.loc[non_interface_c3 + non_interface_c0,['X','Y']]
# pc_coord = pca.transform(meta_non_interface)
# meta_non_interface['PC_1'] = pc_coord
# meta_non_interface.loc[non_interface_c0, 'bin'] = pd.cut(
#     meta_non_interface.loc[non_interface_c0, 'PC_1'], 20,
#     labels = ['bin_' + str(x) for x in np.arange(1,21)]).values
# meta_non_interface.loc[non_interface_c3, 'bin'] = pd.cut(
#     meta_non_interface.loc[non_interface_c3, 'PC_1'], 20,
#     labels = ['bin_' + str(x) for x in np.arange(1,21)]).values
# # Calculate non interface score
# ni_c0_bined = raw_cts.loc[non_interface_c0].groupby(
#     meta_non_interface.loc[non_interface_c0, 'bin']).mean()
# ni_c3_bined = raw_cts.loc[non_interface_c3].groupby(
#     meta_non_interface.loc[non_interface_c3, 'bin']).mean()
# ni_c0_bined = ni_c0_bined.sort_index()
# ni_c3_bined = ni_c3_bined.sort_index()
# c0_c3_intercept = ni_c0_bined.clip(upper = ni_c3_bined, axis=None)
# gene_noninterface_score = 2 * c0_c3_intercept.sum() / (ni_c0_bined.sum() + ni_c3_bined.sum())
# cia_score = (gene_interface_score - gene_noninterface_score).sort_values(ascending=False)

# #%%
# from sklearn.linear_model import LinearRegression
# interface_idx = c3_c0_interface
# bins = meta_interface['bin']
# non_interface_data = neighborhood.loc[
#     [x for x in neighborhood.index if x not in interface_idx]]
# interface_area = neighborhood.loc[interface_idx]
# _ = plt.figure(figsize=(8,8))
# hue_order = sorted(neighborhood.cluster.unique())
# sns.scatterplot(
#     data = interface_area, x = 'X', y = 'Y', s=2, hue=bins, linewidth=0,
#     palette='tab20', hue_order=['bin_{}'.format(x) for x in range(1,21)])
# sns.scatterplot(
#     data = non_interface_data, x = 'X', y = 'Y', s=1, hue='cluster', alpha=0.3,
#     linewidth=0, legend=None,hue_order=hue_order)
# # plot interface proximate line
# x = meta_interface.iloc[:,0].values.reshape(-1,1)
# y_true = meta_interface.iloc[:,1].values.reshape(-1,1)
# lm = LinearRegression().fit(x, y_true)
# y_pred = lm.predict(x)
# sns.lineplot(x=x.flatten(), y=y_pred.flatten())
# plt.legend(bbox_to_anchor = (1,0.5), loc='center left')
# plt.savefig(output_path + 'Interface binned.pdf', bbox_inches = 'tight')
# plt.close()


# # %%
# # %%
# top16_genes = cia_score.index[:16].tolist()
# bot16_genes = cia_score.sort_values().index[:16].tolist()
# # %%
# sns.set_theme(font_scale=2)
# fig, ax = plt.subplots(4,4, figsize=(32,32))
# ax = ax.ravel()
# for i in range(16):
#     sns.lineplot(
#         y = c0_bined[top16_genes[i]], x = np.arange(20), ax = ax[i], color = 'b')
#     sns.lineplot(
#         y = c3_bined[top16_genes[i]], x = np.arange(20), ax = ax[i], color = 'r')
# plt.savefig(output_path + 'top 16 interface genes binned expression.pdf')
# # %%
# sns.set_theme(font_scale=2)
# fig, ax = plt.subplots(4,4, figsize=(32,32))
# ax = ax.ravel()
# for i in range(16):
#     sns.lineplot(
#         y = c0_bined[bot16_genes[i]], x = np.arange(20), ax = ax[i], color = 'b')
#     sns.lineplot(
#         y = c3_bined[bot16_genes[i]], x = np.arange(20), ax = ax[i], color = 'r')
# plt.savefig(output_path + 'bot 16 interface genes binned expression.pdf')
# %%
