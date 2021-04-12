import pandas as pd
import seaborn as sns
import os
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import numpy as np
from sklearn.preprocessing import scale
#%%
def compare_gene_dropouts(raw_cts, denoised_cts, gene):
    _, ax = plt.subplots(1,2,figsize=(8,4))
    np.log2(1+raw_cts[gene]).hist(bins=32, ax=ax[0])
    np.log2(1+denoised_cts[gene]).hist(bins=32, ax=ax[1])
    ax[0].set_title('Raw counts')
    ax[1].set_title('Denoised')

def compare_correlation(raw_cts, denoised_cts, gene1, gene2):
    _, ax = plt.subplots(1,2,figsize=(8,4))
    x1 = np.log2(1+raw_cts[gene1])
    x2 = np.log2(1+denoised_cts[gene1])
    y1 = np.log2(1+raw_cts[gene2])
    y2 = np.log2(1+denoised_cts[gene2])
    sns.scatterplot(x = x1, y = y1, s = 1, ax = ax[0])
    sns.scatterplot(x = x2, y = y2, s = 1, ax = ax[1])
    ax[0].set_title('Raw counts')
    ax[1].set_title('Denoised')

def compare_dropouts(*cts_dfs, labels = None):
    for i, cts in enumerate(cts_dfs):
        if labels is not None:
            label = labels[i]
        else:
            label = ''
        zero_ratio = (cts<=0.01).sum()/cts.shape[0]
        cpm = np.log2(1+cts.apply(lambda x: 1e4*x/x.sum(),axis=1).fillna(0))
        gene_means = cpm.mean()
        # gene_means = np.log2(1+cts.mean())
        plt.scatter(gene_means, zero_ratio,  s=5, label=label, alpha=0.5)
    plt.legend(bbox_to_anchor = (1, 0.5), loc = 'center left', markerscale=2)
    plt.xlabel('Mean log2(1+CPM)')
    plt.ylabel('Zero ratio')

#%%
sc_rna_path = '/project/shared/xiao_wang/projects/MOCCA/data/10x_scrna_seq/'
sc_rna_cts = pd.read_hdf(os.path.join(sc_rna_path, '10k_pbmc.h5df'))
slideseq_counts = pd.read_hdf(
    '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_190921_19/Counts.hdf')
visium_counts = pd.read_csv(
    '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/visium/LN/Counts.txt', 
    sep='\t', index_col=0)
bulk_cts = pd.read_csv(
    '/project/shared/xiao_wang/projects/MOCCA/data/TCGA_BRCA_raw_counts.txt', 
    sep= '\t', index_col=0, header=None)
bulk_cts = bulk_cts.loc[:,bulk_cts.iloc[1] == 'raw_counts']
bulk_cts.columns = bulk_cts.iloc[0]
bulk_cts = bulk_cts.iloc[2:,:].astype(int).transpose()
#%%
# Compare dropout in datasets.
compare_dropouts(
    sc_rna_cts, slideseq_counts, visium_counts, bulk_cts,
    labels = ['SingleCell', 'SlideSeq', 'Visium', 'TCGA_BRCA_bulk'])
plt.savefig(
    '/project/shared/xiao_wang/projects/MOCCA/figures/Dropouts comparison.pdf',
    bbox_inches = 'tight')
#%%

#%%
# norm_cts = raw_counts.apply(lambda x: x/x.sum(), axis=1)
# denoised_counts = pd.read_hdf('denoised_cts.h5df')
# #%%n_genes
# n_samples, n_genes = raw_counts.shape
# zero_ratio_raw = (raw_counts<=0.01).sum()/n_samples
# gene_means_raw = raw_counts.mean()
# zero_ratio_dn = (denoised_counts<=0.01).sum()/n_samples
# gene_means_dn = denoised_counts.mean()
# dropout_df = pd.DataFrame(
#     [zero_ratio_raw.values,gene_means_raw.values, ['Raw']*n_genes], 
#     index = ['Zero ratio', 'Mean', 'Dataset']).transpose()
# dropout_df = dropout_df.append(
#     pd.DataFrame(
#         [zero_ratio_dn.values,gene_means_dn.values, ['Denoised']*n_genes], 
#         index = ['Zero ratio', 'Mean', 'Dataset']).transpose()
# )

# sns.scatterplot(
#     data = dropout_df, x = 'Mean', y = 'Zero ratio', s=5, hue = 'Dataset')
# #%%
# total_counts = pd.DataFrame(index = raw_counts.index, columns=['Raw', 'Denoised'])
# total_counts['Raw'] = raw_counts.sum(axis=1)
# total_counts['Denoised'] = denoised_counts.sum(axis=1)
# #%%
# random_genes = np.random.choice(raw_counts.columns, 10, replace=False)
# exp_comparison = raw_counts[random_genes].stack()
# exp_comparison = pd.DataFrame(exp_comparison, columns = ['Raw'])
# exp_comparison['Denoised'] = denoised_counts[random_genes].stack()
# exp_comparison = exp_comparison.reset_index()
# #%%
# # random_idx = np.random.choice(exp_comparison.index, 100000, replace=False)
# # exp_comparison = exp_comparison.iloc[random_idx]
# exp_comparison['diff'] = exp_comparison.Raw - exp_comparison.Denoised
# sns.boxplot(
#     data = exp_comparison, x = 'level_1', y = 'diff')
# sns.scatterplot(
#     data = exp_comparison, x = 'Raw', y = 'Denoised', s=5)
#%%
# pc_comps = PCA(n_components=10).fit(scale(norm_cts)).components_

# for i in range(10):
#     top_genes = np.argsort(pc_comps[i,:])[::-1][:10]
#     print(norm_cts.columns[top_genes])

# compare_dropouts(raw_counts, denoised_counts, 'Ttr')

# compare_correlation(raw_counts, denoised_counts, 'Ttr', 'Treh')

# raw_counts[['Ttr','Treh']].corr()

# from sklearn.cluster import k_means
# from sklearn.metrics import adjusted_rand_score
# sim_E = pd.read_csv(
#     '/project/shared/xiao_wang/projects/MOCCA/code/Simulation/data/E_noised_simu.txt',
#      sep=' ').transpose()
# c_E = pd.read_csv(
#     '/project/shared/xiao_wang/projects/MOCCA/code/Simulation/data/C_simu.txt', 
#     sep=' ')
# c_E['Y'] = [x[1] for x in c_E.index]
# c_E.index = [x[0] for x in c_E.index]

# clusters = k_means(PCA(0.8).fit_transform(scale(sim_E)),3)[1]

# adjusted_rand_score(c_E.celltype, clusters)

# sns.scatterplot(x = c_E.iloc[:,0], y = c_E.iloc[:,1], hue=c_E.celltype, s=5)
# sns.scatterplot(x = c_E.iloc[:,0], y = c_E.iloc[:,1], hue=clusters.astype(str), s=5)

