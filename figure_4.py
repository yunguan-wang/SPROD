#%%
import sys
sys.path.insert(1,'/home2/s190548/codes/quality_of_life_improvement')
from scrna_utils import *
from collections import defaultdict
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import scipy.ndimage as ndi
from scipy.interpolate import UnivariateSpline
from scipy.spatial.distance import cdist
from sklearn.preprocessing import minmax_scale
from sklearn.cluster import AgglomerativeClustering
from matplotlib_venn import venn2
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
sc.tl.leiden(adata, resolution=0.9)
sc.tl.umap(adata)
adata.write_h5ad('/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/inprogress_cia.h5ad')

denoised_cts = pd.read_hdf(
    '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/' + 
        'subsample_patches/denoised_l30/denoised_stiched.hdf',
    index_col=0)
# denoised_cts = pd.read_hdf(
#     '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/subsample_patches/denoised_counts.hdf',
#     index_col=0)
denoised_adata = denoised_cts.loc[adata.obs_names]
cm_cols = [x for x in denoised_adata.columns if x in adata.var_names]
denoised_adata = denoised_adata[cm_cols].dropna()
# denoised_adata[denoised_adata<=0.02] = 0
denoised_adata = sc.AnnData(denoised_adata)
denoised_adata.var_names_make_unique()
sc.pp.normalize_per_cell(denoised_adata, 1000)
sc.pp.log1p(denoised_adata)
denoised_adata.raw = denoised_adata
sc.pp.highly_variable_genes(denoised_adata, min_disp=0.5, min_mean=0.05)
# sc.pl.highly_variable_genes(denoised_adata)
sc.pp.scale(denoised_adata)
sc.tl.pca(denoised_adata, n_comps=50)
sc.pp.neighbors(denoised_adata, n_pcs=25)
sc.tl.leiden(denoised_adata, resolution=0.5)

saver_cts = pd.read_hdf(
    '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/saver_cts.h5df',
    index_col=0)
# saver_cts = pd.read_hdf(
#     '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/subsample_patches/denoised_counts.hdf',
#     index_col=0)
saver_adata = saver_cts.loc[adata.obs_names]
cm_cols = [x for x in saver_adata.columns if x in adata.var_names]
saver_adata = saver_adata[cm_cols].dropna()
saver_adata = sc.AnnData(saver_adata)
saver_adata.var_names_make_unique()
sc.pp.normalize_per_cell(saver_adata, 1000)
sc.pp.log1p(saver_adata)
saver_adata.raw = saver_adata

denoised_adata.write_h5ad(
    '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/denoised_cts.h5ad'
)
saver_adata.write_h5ad(
    '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/saver_cts.h5ad'
)
#%%
output_path = '/project/shared/xiao_wang/projects/MOCCA/figures/'
np.random.seed(0)
sns.set_theme(font_scale=3, style='white', context='paper')
sc.settings.figdir = output_path
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
denoised_adata = sc.read_h5ad(
    '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/denoised_cts.h5ad'
)
saver_adata = sc.read_h5ad(
    '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/saver_cts.h5ad'
)

#%%
# annotate cells near CA1 region
# plt.plot(xs, ys, 'g', lw=3)
# plt.scatter(x_ori,y_ori, c = 'r', s=10)
spline_raw = spot_meta.copy()
spline_raw['cluster'] = 'c_' + adata.obs['leiden'].astype(str)
xs,ys,x_ori, y_ori = find_cluster_spline('c_9', spline_raw, 1500,4500, num_bins=50)
spline = np.array([xs,ys]).T
spot_dist = cdist(spline_raw.iloc[:,:2],spline)
spot_spline_pos = np.argmin(spot_dist, axis=1)
spot_dist = np.min(spot_dist, axis=1)
spline_raw['Soma-Proximal Axis (um)'] = spot_dist
spline_raw['spline_y'] = ys[spot_spline_pos]
spline_raw['Soma-Proximal Axis (um)'] = spline_raw['Soma-Proximal Axis (um)'] * (((spline_raw['spline_y'] - spline_raw['Y'])>0)*-2 + 1)
spline_raw = spline_raw[
    (spline_raw['Soma-Proximal Axis (um)']<=665)&(spline_raw['Soma-Proximal Axis (um)']>=-265)].copy()
spline_raw = spline_raw[
    (spline_raw.X>=1500) & (spline_raw.X<=4500) & (spline_raw.Y<=4200)]
spline_raw['CA1_layer'] = 'soma'
spline_raw.loc[spline_raw['Soma-Proximal Axis (um)']<=-65, 'CA1_layer'] = 'basal neuropil'
spline_raw.loc[spline_raw['Soma-Proximal Axis (um)']>=65, 'CA1_layer'] = 'proximal neuropil'
spline_raw = spline_raw[~spline_raw.cluster.isin(['c_5','c_8','c_12','c_13','c_14'])]
spline_raw['Soma-Proximal Axis (um)'] = spline_raw['Soma-Proximal Axis (um)']/2
#%%
# annotate cells near CA1 region
# plt.plot(xs, ys, 'g', lw=3)
# plt.scatter(x_ori,y_ori, c = 'r', s=10)
spline_denoised = spot_meta.copy()
spline_denoised['cluster'] = 'c_' + denoised_adata.obs['leiden'].astype(str)
xs,ys,x_ori, y_ori = find_cluster_spline('c_5', spline_denoised, 1500,4500, num_bins=50)
spline = np.array([xs,ys]).T
spot_dist = cdist(spline_denoised.iloc[:,:2],spline)
spot_spline_pos = np.argmin(spot_dist, axis=1)
spot_dist = np.min(spot_dist, axis=1)
spline_denoised['Soma-Proximal Axis (um)'] = spot_dist
spline_denoised['spline_y'] = ys[spot_spline_pos]
spline_denoised['Soma-Proximal Axis (um)'] = spline_denoised['Soma-Proximal Axis (um)'] * (((spline_denoised['spline_y'] - spline_denoised['Y'])>0)*-2 + 1)
spline_denoised = spline_denoised[
    (spline_denoised['Soma-Proximal Axis (um)']<=665)&(spline_raw['Soma-Proximal Axis (um)']>=-265)].copy()
spline_denoised = spline_denoised[
    (spline_denoised.X>=1500) & (spline_denoised.X<=4500) & (spline_denoised.Y<=4200)]
spline_denoised['CA1_layer'] = 'soma'
spline_denoised.loc[spline_denoised['Soma-Proximal Axis (um)']<=-65, 'CA1_layer'] = 'basal neuropil'
spline_denoised.loc[spline_denoised['Soma-Proximal Axis (um)']>=65, 'CA1_layer'] = 'proximal neuropil'
spline_denoised = spline_denoised[~spline_denoised.cluster.isin(['c_8','c_9','c_11','c_12','c_13'])]
spline_denoised['Soma-Proximal Axis (um)'] = spline_denoised['Soma-Proximal Axis (um)']/2
# %%
for i, spline in enumerate([spline_raw, spline_denoised]):
    _ = plt.figure(figsize=(8,8))
    sns.scatterplot(
        data = spline,
        x = 'X',
        y = 'Y',
        hue = 'CA1_layer', s = 10, linewidth=0,
        palette=sns.color_palette('spring',3))
    plt.plot(xs, ys, 'r', lw=5)
    plt.legend(bbox_to_anchor = (1,0.5), loc='center left', markerscale=3)
    plt.savefig(
        output_path + 'CA1 basal to proximal axis {}.pdf'.format(['Raw','Sprod'][i]),
        bbox_inches = 'tight')

# background = spot_meta[
#     (spot_meta.X>=1500) &
#     (spot_meta.X<=4500) &
#     (spot_meta.Y<=4600) &
#     (spot_meta.Y>=2000)]
# background.drop(spline_raw.index, axis=0, inplace=True)
# sns.scatterplot(
#     data = background,
#     x = 'X',
#     y = 'Y',
#     color = 'grey', s = 3, linewidth=0)
#%%
# patch_meta_folder = '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/patches'
# meta_fns = [x for x in os.listdir(patch_meta_folder) if 'Spot_metadata.csv' in x]
# denoised_adata.obs['batch'] = ''
# for meta_fn in meta_fns:
#     patch_name = meta_fn.replace('_Spot_metadata.csv','')
#     meta_fn = os.path.join(patch_meta_folder, meta_fn)
#     patch_meta = pd.read_csv(meta_fn, index_col=0)
#     patch_bc = patch_meta[patch_meta.patch_core == patch_name].index
#     patch_bc = [x for x in patch_bc if x in denoised_adata.obs_names]
#     denoised_adata.obs.loc[patch_bc, 'batch'] = patch_name
# sc.external.pp.harmony_integrate(denoised_adata, 'batch')
# sc.pp.neighbors(denoised_adata, n_pcs=25, use_rep='X_pca_harmony')
# sc.tl.leiden(denoised_adata, resolution=0.25)
# %%
ca1_genes = list(set(pd.read_csv(
    '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/CA1_genes.tsv',
    index_col=0, sep='\t').iloc[:,0]))
ca1_genes = [x for x in ca1_genes if x in adata.var_names]
ca1_genes_slideseq = sc.get.rank_genes_groups_df(adata, 'c_9')
ca1_genes_slideseq = ca1_genes_slideseq[
    (ca1_genes_slideseq.logfoldchanges>=1) & 
    (ca1_genes_slideseq.pvals_adj<=0.05)].names.values
ca1_genes = [x for x in ca1_genes if x in ca1_genes_slideseq]
ca1_genes = list(set(ca1_genes) | set(ca1_genes_slideseq[:100]))
#%%
adata_raw = adata.raw.to_adata().to_df()
denoised_raw = denoised_adata.raw.to_adata().to_df()
saver_raw = saver_adata.raw.to_adata().to_df()
# %%
for gene in ['Camk2a','Hpca','Map2']:
    _, axes = plt.subplots(3, 1, figsize = (6,18), sharex=True)
    axes = axes.ravel()
    for i, df in enumerate([adata_raw, saver_raw, denoised_raw]):
        spline = [spline_raw, spline_raw, spline_denoised][i].copy()
        genevec = df.loc[spline.index, gene]
        title = gene
        if isinstance(gene, list):
            genevec = df[gene].mean(axis=1)
            title = 'Mean of CA1 genes'
        spline['log1p CPK'] = genevec.values
        _ = sns.regplot(
            data = spline, x='Soma-Proximal Axis (um)', y='log1p CPK', order=10,
            scatter_kws={"s": 2, 'color':'grey'}, ci=95, ax = axes[i], color='magenta')
        axes[i].set_title(title + '_' + ['Raw', 'Saver', 'Sprod'][i])
        if i != 2:
            axes[i].set_xlabel('')
            axes[i].set_xticks([-100,0,100,200,300])
            axes[i].set_xticklabels([-100,0,100,200,300], fontsize=12)
    # plt.savefig(output_path + 'CA1 gene profile {}.pdf'.format(title))

# %%
spline_raw, spline_denoised = spline_raw.align(
    spline_denoised, join='inner', axis=0)
# Define colors
for gene in ['Camk2a','Hpca','Map2']:
    _, axes = plt.subplots(3, 1, figsize = (6,18), sharex=True)
    axes = axes.ravel()
    for i, df in enumerate([adata_raw, saver_raw, denoised_raw]):
        spline = [spline_raw, spline_raw, spline_denoised][i].copy()
        genevec = df.loc[spline.index, gene]
        title = gene
        if isinstance(gene, list):
            genevec = df[gene].mean(axis=1)
            title = 'Mean of CA1 genes'
        spline['log1p CPK'] = genevec.values
        sns.scatterplot(
            data = spline,
            x = 'X',
            y = 'Y',
            hue = 'log1p CPK', linewidth=0, s=10, ax = axes[i])
        axes[i].set_title(title + '_' + ['Raw', 'Saver', 'Sprod'][i])
    # plt.savefig(output_path + 'CA1 gene expression {}.pdf'.format(title))

#%%
# ==============================================================================
# ==============================================================================
# ==============================================================================
# Dendrite enrichment
for i, deg_adata in enumerate([adata, denoised_adata]):
    spline = [spline_raw, spline_denoised][i].copy()
    deg_adata.obs['dendrite_enrichment'] = ''
    deg_adata.obs.loc[
        spline[spline.CA1_layer == 'proximal neuropil'].index,
        'dendrite_enrichment'
        ] = 'proximal neuropil'
    deg_adata.obs.loc[
        spline[spline.CA1_layer == 'soma'].index,
        'dendrite_enrichment'
        ] = 'soma'
    sc.tl.rank_genes_groups(
        deg_adata, groupby = 'dendrite_enrichment',
        groups=['proximal neuropil'], reference='soma',
        n_genes=-1, method='wilcoxon')
    dendrite_df = sc.get.rank_genes_groups_df(
        deg_adata, 'proximal neuropil')
    output_prefix = ['Raw', 'Sprod'][i]
    dendrite_df[
        (abs(dendrite_df.logfoldchanges) >= 1) & 
        (dendrite_df.pvals_adj <= 0.05)
        ].to_csv('~/temp/{}_dendrite_enriched.csv'.format(output_prefix))
#%%
# ==============================================================================
# ==============================================================================
# ==============================================================================
# SpatialDE analysis
# %%
roi_raw = spline_raw.index
roi_de = spline_denoised.index
sde_raw = adata.raw.to_adata().to_df().loc[roi_raw]
sde_denoised = denoised_adata.raw.to_adata().to_df().loc[roi_de]
#%%
res_raw = pd.read_csv(
    "/project/shared/xiao_wang/projects/MOCCA/data/" + 
    "Sprod_ready_data/slideseq/Puck_200115_08/SpatialDE_raw.csv", index_col=0)
res_raw = res_raw[res_raw.g.apply(lambda x: 'mt-' not in x)]
res_de = pd.read_csv(
    "/project/shared/xiao_wang/projects/MOCCA/data/" + 
    "Sprod_ready_data/slideseq/Puck_200115_08/SpatialDE_de.csv",index_col=0)
#%%
gene_clusters = []
for i, sde_res in enumerate([res_raw, res_de]):
    np.random.seed(0)
    sig_genes = sde_res[sde_res.qval<0.01].g.values
    sig_gene_df = [sde_raw, sde_denoised][i][sig_genes].T
    cluster_model = AgglomerativeClustering(16, affinity='correlation',linkage='complete')
    _clusters = cluster_model.fit_predict(sig_gene_df)
    gene_clusters.append(_clusters)
    fig, ax = plt.subplots(4,4, sharex=True, sharey=True, figsize = (15,15))
    ax = ax.ravel()
    for j,gc in enumerate(np.unique(_clusters)):
        gc_df = sig_gene_df[_clusters==gc].mean()
        sns.scatterplot(
            data = spot_meta.loc[sig_gene_df.columns],
            x = 'X',
            y = 'Y',
            hue = gc_df, linewidth=0, palette='Reds', s=5,
            legend=None, ax = ax[j])
        ax[j].set_title('Gene cluster {}'.format(j))
    plt.savefig(
        output_path + 'SpatialDE_gene_clusters_pattern_{}.pdf'.format(
            ['Raw','Sprod'][i]
        ))
    plt.close()
# %%
de_raw_genes = gene_clusters[0]
de_raw_genes = [True if x in [7] else False for x in de_raw_genes ]
de_sprod_genes = gene_clusters[1]
de_sprod_genes = [True if x in [9] else False for x in de_sprod_genes ]
de_sprod_genes = (res_de[res_de.qval<0.01].g[de_sprod_genes]).values
de_raw_genes = (res_raw[res_raw.qval<0.01].g[de_raw_genes]).values
venn2([set(de_raw_genes), set(de_sprod_genes)],['SpatialDE_raw','SpatialDE_sprod'])
plt.savefig(output_path + 'Venn_raw_sprod_spatialde.pdf', bbox_inches = 'tight')
plt.close()
pd.DataFrame(
    [de_raw_genes,de_sprod_genes], index=['spatialde_raw', 'spatialde_sprod']
    ).T.to_csv("/project/shared/xiao_wang/projects/MOCCA/data/" + 
    "Sprod_ready_data/slideseq/Puck_200115_08/spatialde_dendrite_enriched.csv")
# %%
de_enriched = pd.read_csv("/project/shared/xiao_wang/projects/MOCCA/data/" + 
    "Sprod_ready_data/slideseq/Puck_200115_08/Dendrite_enriched_genes.csv")
g1 = de_enriched.Spatialde_sprod.dropna().values
g2 = de_enriched.Dendrite_DE_sprod.dropna().values
taget_genes = list(set(g1) |set(g2))
sns.set_theme(font_scale=1.5, style='white', context='paper')

for _adata, _spline in zip([adata, denoised_adata],[spline_raw, spline_denoised]):
    _spline = _spline.sort_values('Soma-Proximal Axis (um)')
    plot_adata = _adata[_spline.index]
    plot_adata.obs['Region'] = _spline['CA1_layer']
    plot_adata.obs['Region'] = plot_adata.obs['Region'].astype('category')
    plot_adata.obs['Region'].cat.reorder_categories(
        ['basal neuropil', 'soma', 'proximal neuropil'], inplace=True)
    prefix = 'Raw' if i==0 else 'Sprod'
    sc.pl.tracksplot(
        plot_adata, g1, 'Region', dendrogram=False,yticklabels=[],
        figsize=(8,16), show=False, 
        save = 'Dendritic_genes_' + prefix + '_DE.pdf')
    sc.pl.tracksplot(
        plot_adata, g2, 'Region', dendrogram=False,yticklabels=[],
        figsize=(8,16),show=False, 
        save = 'Dendritic_genes_' + prefix + '_spatialDE.pdf')
    sc.pl.tracksplot(
        plot_adata, taget_genes, 'Region', dendrogram=False,yticklabels=[],
        figsize=(8,24),show=False, 
        save = 'Dendritic_genes_' + prefix + '_union.pdf')
    i += 1
#%%
from scipy.stats import hypergeom
from collections import defaultdict
gl = {}
for col in de_enriched:
    genes = de_enriched[col].drop_duplicates().dropna().values
    gl[col] = genes
dict_keys = list(gl.keys())

for k1 in dict_keys[:4]:
    g = gl[k1]
    N = len(g)
    for k2 in dict_keys[4:]:
        g_r = gl[k2]
        n = len(g_r)
        x = len(set(g) & set(g_r))
        print(
            'Pvalue of {} vs {} : {:.2e}'.format(
                k1,k2, hypergeom(20000, n, N).sf(x-1)))
# %%
import gseapy as gs
import scipy.cluster.hierarchy as sch
import matplotlib.lines as mlines
from matplotlib.collections import PatchCollection

def dotted_heatmap(
    s, c, xlab, ylab, swap_axis = False, normalize_c = True, cmap='coolwarm',
    figsize=None):
    '''
    Dotted heatmap showing two layers of information as dot sizes and dot colors. \
        heatmap has a shape of xlab by ylab. s, c are dataframes of xlab by ylab. \
        clustering is always done on the x axis.
    '''
    # sort pathways based on c matrix, which is NES
    lm = sch.linkage(c, 'complete')
    orders = sch.dendrogram(lm, no_plot=True)['leaves']
    xlab = c.index[orders]
    c = c.loc[xlab]
    s = s.loc[xlab]
    # Normalize colors
    if normalize_c:
        c_max1 = c.max().max()
        c_max2 = abs(c.min().min())
        c[c > min(c_max1, c_max2)] = min(c_max1, c_max2)
        c[c < -min(c_max1, c_max2)] = -min(c_max1, c_max2)
    # normalize sizes
    R = s/s.max().max()/2
    if swap_axis == True:
        xlab, ylab = ylab, xlab
        x, y = np.meshgrid(np.arange(len(xlab)), np.arange(len(ylab)), indexing='xy')
        circles = [plt.Circle((i,j), radius=r) for r, i, j in zip(R.values.flat, x.flat, y.flat)]
        fig, ax = plt.subplots(figsize= figsize or (len(xlab)/2,max(5,len(ylab)/3)))
    else:
        x, y = np.meshgrid(np.arange(len(xlab)), np.arange(len(ylab)), indexing='ij')
        circles = [plt.Circle((i,j), radius=r) for r, i, j in zip(R.values.flat, x.flat, y.flat)]
        fig, ax = plt.subplots(figsize=figsize or (len(xlab)/3,max(5,len(ylab)/3)))
    ax.axis('equal')
    col = PatchCollection(circles, array=c.values.flatten(), cmap=cmap)
    ax.add_collection(col)
    ax.set_xticks(np.arange(len(xlab)))
    if swap_axis:
        ax.set_xticklabels(xlab, fontsize='large')
    else:
        ax.set_xticklabels(xlab, fontsize='large', rotation=90)
    ax.set_yticks(np.arange(len(ylab)))
    ax.set_yticklabels(ylab, fontsize='large')
    ax.grid(which='minor')
    ax.set_xlim(-0.5, len(xlab))
    ax.set_ylim(-0.5, len(ylab))
    fig.colorbar(
        col, use_gridspec = False, label = '-log10 (adj-pval)', panchor=(1.2,0.4))
    # make legends
    handles = []
    for j in [0.25, 0.5, 1]:
        handles.append(
            mlines.Line2D(
                [], [], color='red', marker='.', linestyle='None',
                markersize=5 * j/0.25, label=''))
    labels = [int(x) for x in [s.min().min(), s.median().median(), s.max().max()]]
    ax.legend(
        handles=handles , labels = labels,  title = 'Number of genes',
        bbox_to_anchor = (1.4,0.9,0.05,0.1), loc = 'upper left')
    return ax

# %%
raw = de_enriched.iloc[:,[0,2]].values.flatten()
raw = [x for x in raw if isinstance(x,str)]
sprod = de_enriched.iloc[:,[1,3]].values.flatten()
sprod = [x for x in sprod if isinstance(x,str)]
# %%
libs = ['GO_Biological_Process_2021', 'GO_Cellular_Component_2021']
res = [pd.DataFrame(), pd.DataFrame()]
for i, gl in enumerate([raw, sprod]):
    for lib in libs:
        gl = [x for x in gl if ('Rps' not in x) and ('Rpl' not in x)]
        _res = gs.enrichr(gl, lib, 'mouse', no_plot=True).res2d
        _res = _res.sort_values('Combined Score', ascending=False)
        _res = _res[_res['Adjusted P-value']<=0.05].set_index('Term')
        _res['Overlap'] = [int(x.split('/')[0]) for x in _res['Overlap']]
        res[i] = res[i].append(_res)

#%%
dot_hmap_data = res[1].sort_values(
    'Combined Score', ascending=False).iloc[:20][['Adjusted P-value','Overlap']]
dot_hmap_data = dot_hmap_data.merge(
    res[0][['Adjusted P-value','Overlap']], 
    left_index = True, right_index = True, how='left').fillna(1)
dot_hmap_data.columns = ['Sprod','s_overlap','Raw','r_overlap']
dot_hmap_data.iloc[:,[0,2]] = -np.log10(dot_hmap_data.iloc[:,[0,2]])
dot_hmap_data.index = [x.split(' (')[0] for x in dot_hmap_data.index]
#
# %%
_ = plt.figure(figsize=(8,16))
ax = dotted_heatmap(
    dot_hmap_data.iloc[:,[1,3]], dot_hmap_data.iloc[:,[0,2]],
    dot_hmap_data.index, 
    ['Srpod','Raw'], normalize_c=False, 
    swap_axis=True,cmap='Reds',
    )
ax.set_xticks([-0.5,0.5,1])
ax.set_xticklabels(['Sprod', '', 'Raw'], rotation=45)
plt.savefig(
    "/project/shared/xiao_wang/projects/MOCCA/figures/Sprod_raw_GoCC_terms.svg",
    bbox_inches='tight')
#%%
libs = ['PanglaoDB_Augmented_2021']
res = [pd.DataFrame(), pd.DataFrame()]
for i, gl in enumerate([raw, sprod]):
    for lib in libs:
        gl = [x for x in gl if ('Rps' not in x) and ('Rpl' not in x)]
        _res = gs.enrichr(gl, lib, 'mouse', no_plot=True).res2d
        _res = _res.sort_values('Combined Score', ascending=False)
        _res = _res[_res['Adjusted P-value']<=0.05].set_index('Term')
        _res['Overlap'] = [int(x.split('/')[0]) for x in _res['Overlap']]
        res[i] = res[i].append(_res)
dot_hmap_data = res[1].sort_values(
    'Combined Score', ascending=False).iloc[:10][['Adjusted P-value','Overlap']]
dot_hmap_data = dot_hmap_data.merge(
    res[0][['Adjusted P-value','Overlap']], 
    left_index = True, right_index = True, how='left').fillna(1)
dot_hmap_data.columns = ['Sprod','s_overlap','Raw','r_overlap']
dot_hmap_data.iloc[:,[0,2]] = -np.log10(dot_hmap_data.iloc[:,[0,2]])
dot_hmap_data.index = [x.split(' (')[0] for x in dot_hmap_data.index]
dot_hmap_data = dot_hmap_data.rename({'M?ller Cells':'Muller Cells'})
ax = dotted_heatmap(
    dot_hmap_data.iloc[:,[1,3]], dot_hmap_data.iloc[:,[0,2]],
    dot_hmap_data.index, 
    ['Srpod','Raw'], normalize_c=False, 
    swap_axis=True,cmap='Reds',
    figsize=(3,8),
    )
ax.set_xticks([-0.5,0.5,1])
ax.set_xticklabels(['Sprod', '', 'Raw'], rotation=45)
plt.savefig(
    "/project/shared/xiao_wang/projects/MOCCA/figures/Sprod_raw_PangloDB.svg",
    bbox_inches='tight')
# %%
libs = ['DisGeNET']
res = [pd.DataFrame(), pd.DataFrame()]
for i, gl in enumerate([raw, sprod]):
    for lib in libs:
        gl = [x for x in gl if ('Rps' not in x) and ('Rpl' not in x)]
        _res = gs.enrichr(gl, lib, 'mouse', no_plot=True).res2d
        _res = _res.sort_values('Combined Score', ascending=False)
        _res = _res[_res['Adjusted P-value']<=0.05].set_index('Term')
        _res['Overlap'] = [int(x.split('/')[0]) for x in _res['Overlap']]
        res[i] = res[i].append(_res)
dot_hmap_data = res[1].sort_values(
    'Adjusted P-value').iloc[:10][['Adjusted P-value','Overlap']]
dot_hmap_data = dot_hmap_data.merge(
    res[0][['Adjusted P-value','Overlap']], 
    left_index = True, right_index = True, how='left').fillna(1)
dot_hmap_data.columns = ['Sprod','s_overlap','Raw','r_overlap']
dot_hmap_data.iloc[:,[0,2]] = -np.log10(dot_hmap_data.iloc[:,[0,2]])
dot_hmap_data.index = [x.split(' (')[0] for x in dot_hmap_data.index]
dot_hmap_data = dot_hmap_data.rename({'M?ller Cells':'Muller Cells'})
ax = dotted_heatmap(
    dot_hmap_data.iloc[:,[1,3]], dot_hmap_data.iloc[:,[0,2]],
    dot_hmap_data.index, 
    ['Srpod','Raw'], normalize_c=False, 
    swap_axis=True,cmap='Reds',
    figsize=(3,8),
    )
ax.set_xticks([-0.5,0.5,1])
ax.set_xticklabels(['Sprod', '', 'Raw'], rotation=45)
plt.savefig(
    "/project/shared/xiao_wang/projects/MOCCA/figures/Sprod_raw_DisGeNET.svg",
    bbox_inches='tight')
# %%
import holoviews as hv
hv.extension('matplotlib')
# %%
edges = pd.read_csv("/project/shared/xiao_wang/projects/MOCCA/data/" + 
    "Sprod_ready_data/slideseq/Puck_200115_08/rbp_edges.csv")
sankey = {
    'Sprod_Not RBP Binding': 0,
    'Raw_Not RBP Binding': 0,
    'Raw_RBP Binding': 0,
    'Sprod_RBP Binding': 0}
rbp_genes = []
raw_genes = []
sprod_genes = []
for rn, row in edges.set_index('Target').iterrows():
    if row['Node_attr'] == 'RBP':
        rbp = rn
        gene = row['Source']
        if row['Type'] == 'Sprod':
            if gene not in sprod_genes:
                sankey['Sprod_RBP Binding'] += 1
                sprod_genes.append(gene)
        else:
            if gene not in raw_genes:
                sankey['Raw_RBP Binding'] += 1
                raw_genes.append(gene)
        key_v = 'RBP Binding_' + rbp
        if key_v not in sankey.keys():
            sankey['RBP Binding_' + rbp] = 1
        else:
            sankey['RBP Binding_' + rbp] +=1
        rbp_genes.append(gene)
    else:
        gene = rn
        if gene in rbp_genes:
            continue
        if row['Source'] == 'SPROD':
            sankey['Sprod_Not RBP Binding'] += 1
        else:
            sankey['Raw_Not RBP Binding'] += 1
sankey['Not RBP Binding_ '] = 0.001
sankey = pd.DataFrame.from_dict(sankey, orient = 'index')
sankey.columns = ['value']
sankey['target'] = [x.split('_')[1] for x in sankey.index]
sankey = sankey[['target','value']]
sankey.index = [x.split('_')[0] for x in sankey.index]
g = hv.Sankey(sankey.reset_index(), label='RNA-RBP interactions')
g.opts(label_position='left', edge_color='target', node_color='index', cmap='tab10')
hv.save(g, "/project/shared/xiao_wang/projects/MOCCA/figures/RBP-interactions.svg")

# %%
edges = pd.read_csv("/project/shared/xiao_wang/projects/MOCCA/data/" + 
    "Sprod_ready_data/slideseq/Puck_200115_08/rbp_edges.csv")
df = pd.DataFrame(index = ['Raw','Sprod'])
df.loc['Raw','RBP Binding'] = edges[
    (edges.Node_attr=='RBP') & (edges.Type=='Raw')]['Source'].nunique()
df.loc['Raw','Not RBP Binding'] = (
    edges[(edges.Node_attr=='Gene') & (edges.Source=='RAW')].shape[0] - 
    df.loc['Raw','RBP Binding'])
df.loc['Sprod','RBP Binding'] = edges[
    (edges.Node_attr=='RBP') & (edges.Type=='Sprod')]['Source'].nunique()
df.loc['Sprod','Not RBP Binding'] = (
    edges[(edges.Node_attr=='Gene') & (edges.Source=='SPROD')].shape[0] - 
    df.loc['Sprod','RBP Binding'])
rbp_edges = edges[edges.Node_attr=='RBP']
for rbp in rbp_edges.Target.unique():
    df_rbp = rbp_edges[rbp_edges.Target==rbp]
    for j in ['Raw','Sprod']:
        df.loc[j,rbp] = (df_rbp['Type']==j).sum()

_, ax = plt.subplots(1,2, figsize=(16,8))
df.iloc[:,:2].plot.bar(stacked=True,fontsize = 16, ax=ax[0])
df.iloc[:,2:].T.plot.bar(stacked=True,fontsize = 16, ax=ax[1])
plt.savefig(
    "/project/shared/xiao_wang/projects/MOCCA/figures/RBP-interactions stacked bar.svg",
    bbox_inches = 'tight')
# %%
