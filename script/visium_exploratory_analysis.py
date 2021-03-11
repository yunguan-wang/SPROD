#%%
from numpy.core.fromnumeric import mean
import pandas as pd
from skimage import io, morphology, img_as_float
import numpy as np
import matplotlib.pyplot as plt
import os
from sklearn.preprocessing import minmax_scale

def marker_overlay(gene, spot_meta, cts, img, normalizing = False, img_type='HE'):
    gene_exp_mask = np.zeros(img.shape[:2], 'uint8')
    r = spot_meta['Spot_radius'][0]
    bbox_mask = morphology.disk(r)
    # c1_median = spot_meta['C1_median'].median()
    for id, row in spot_meta.iterrows():
        x, y = row[['X','Y']].astype(int)
        if gene in ['C1_median','C2_median']:
            gene_exp_level = row[gene].astype(int)
        else:
            # if normalizing:
            #     c1 = row.C1_median
            #     gene_exp_level = cts.loc[id,gene]
            #     gene_exp_level = np.max((0,gene_exp_level + cd1c_ratio * (c1-c1_median)))
            # else:
            gene_exp_level = cts.loc[id,gene].astype(int)
        gene_exp_mask[x-r:x+r+1,y-r:y+r+1] = gene_exp_level
        gene_exp_mask[x-r:x+r+1,y-r:y+r+1] = gene_exp_mask[x-r:x+r+1,y-r:y+r+1]*bbox_mask
    if normalizing:
        output_fn = '{} overlay normalized.pdf'.format(gene)
        gene_exp_mask = 255*minmax_scale(gene_exp_mask.flatten()).reshape(gene_exp_mask.shape)
        gene_exp_mask = gene_exp_mask.astype('uint8')
    else:
        output_fn = '{} overlay.pdf'.format(gene)
    if img_type == 'IF':
        exp_mask = np.zeros(
            (gene_exp_mask.shape[0],gene_exp_mask.shape[1],3), dtype='uint8')
        exp_mask[:,:,1] = gene_exp_mask
        gene_exp_mask = exp_mask.copy()
    _, ax = plt.subplots(figsize=(64,36))
    io.imshow(img, ax=ax)
    io.imshow(gene_exp_mask, ax=ax,alpha=0.5)
    plt.savefig(output_fn)
    plt.close()
    return gene_exp_mask
#%%
# working_dir = '/project/shared/xiao_wang/projects/MOCCA/data/Visium/breast_cancer/'
# os.chdir(working_dir)

# img_tif = [x for x in os.listdir() if 'tif' in x][0]
# cts = pd.read_csv('Counts.txt',sep='\t', index_col=0)
# cts = cts/cts.sum()*1e4
# cts = np.log2(cts+1)
# spot_meta = pd.read_csv('Spot_metadata.csv', index_col=0)
# img = io.imread(img_tif)
# # alightment
# spot_meta = spot_meta.loc[cts.index]

# tc = cts.sum().sum()
# t_cd1c = cts['CD1C'].sum()
# cd1c_ratio = t_cd1c/tc*1e4

# spot_meta.C1_median - spot_meta.C1_median.median()

# gm = marker_overlay('VIM',spot_meta, cts)
# gm = marker_overlay('C1_median',spot_meta, cts)


# gene_exp_mask = np.zeros(img.shape[:2], 'uint8')
# r = spot_meta['Spot_radius'][0]
# bbox_mask = morphology.disk(r)
# for id, row in spot_meta.iterrows():
#     x, y = row[['X','Y']].astype(int)
#     gene_exp_level = row.C1_median
#     gene_exp_mask[x-r:x+r+1,y-r:y+r+1] = gene_exp_level
#     gene_exp_mask[x-r:x+r+1,y-r:y+r+1] = gene_exp_mask[x-r:x+r+1,y-r:y+r+1]*bbox_mask
# gene_exp_mask = (
#     255*(gene_exp_mask-gene_exp_mask.min())/np.ptp(gene_exp_mask)).astype('uint8')
# _ = plt.figure(figsize=(64,36))
# io.imshow(img)
# io.imshow(gene_exp_mask, alpha=0.5, cmap='gray')
# plt.savefig('C1 overlay.pdf')
# plt.close()

# io.imshow(gene_exp_mask)
#%%
oc = '/project/shared/xiao_wang/projects/MOCCA/data/Visium/ovarian_cancer_immune/'
oc_wt = '/project/shared/xiao_wang/projects/MOCCA/data/Visium/ovarian_cancer_WT/'
cts_oc = pd.read_csv(os.path.join(oc, 'Counts.txt'),sep='\t',index_col=0)
cts_oc_wt = pd.read_csv(os.path.join(oc_wt, 'Counts.txt'),sep='\t',index_col=0)
spot_meta = pd.read_csv(os.path.join(oc,'Spot_metadata.csv'), index_col=0)
os.chdir(oc)
img_tif = [x for x in os.listdir(oc) if 'tif' in x][0]

# %%
oc = '/project/shared/xiao_wang/projects/MOCCA/data/Visium/LN/'
os.chdir(oc)
cts = pd.read_csv(os.path.join(oc, 'Counts.txt'),sep='\t',index_col=0)
spot_meta = pd.read_csv(os.path.join(oc,'Spot_metadata.csv'), index_col=0)
img_tif = [x for x in os.listdir(oc) if 'tif' in x][0]
img = img_as_float(io.imread(img_tif))
img = (255 * img).astype('uint8')
cts = cts/cts.sum()*1e4
cts = np.log2(cts+1).fillna(0)
# cts_oc_wt = cts_oc_wt/cts_oc_wt.sum()*1e4
# cts_oc_wt = np.log2(cts_oc_wt+1)
gm = marker_overlay('IGHD',spot_meta, cts, img, normalizing=True, img_type='IF')
gm = marker_overlay('CD3D',spot_meta, cts, img, normalizing=True, img_type='IF')
gm = marker_overlay('MS4A1',spot_meta, cts, img, normalizing=True, img_type='IF')
gm = marker_overlay('CD1C',spot_meta, cts, img, normalizing=True, img_type='IF')
#%%
from matplotlib.backends.backend_pdf import PdfPages
os.chdir(oc)
n_pages = 4
cutoffs=(0.1,0.2,0.3,0.4)
with PdfPages('CD45 correlations.pdf') as pdf:
    for _, cutoff in enumerate(cutoffs):
        fig, axes = plt.subplots(2,2, figsize=(8,8), sharex=True, sharey=True)
        for i, cts in enumerate([cts_oc, cts_oc_wt]):
            keep = cts.mean(axis=1)>cutoff
            kept_cts = cts[keep].copy()
            kept_spot_meta = spot_meta[keep].copy()
            
            kept_cts = (kept_cts.T/kept_cts.T.sum() * 1e6).T.fillna(0)
            kept_cts = np.log2(kept_cts+1)

            cors = np.array([
                np.corrcoef(
                    kept_cts["PTPRC"], kept_cts[x]
                    )[0,1] for x in kept_cts.columns])
            genes = kept_cts.columns[cors>0.1]
            # Correlation with PTPRC correlated genes
            cor = np.corrcoef(
                kept_spot_meta['C1_median'],
                kept_cts[genes].mean(axis=1))[0,1]
            axes[i,0].scatter(
                kept_spot_meta.C1_median,kept_cts[genes].mean(axis=1), s=1)
            axes[i,0].set_title(
                'N_counts cutoff={}, mean of PTPRC cor genes\ncor={:.2f}'.format(
                    cutoff, cor))
            if i == 0:
                axes[i,0].set_ylabel('Targeted Log2(CPM+1)')
            else:
                axes[i,0].set_ylabel('Whole transcriptome Log2(CPM+1)')
            # Correlation with PTPRC itself
            cor = np.corrcoef(
                kept_spot_meta['C1_median'],
                kept_cts["PTPRC"])[0,1]
            axes[i,1].scatter(kept_spot_meta.C1_median,kept_cts.PTPRC, s=1)
            axes[i,1].set_title(
                'N_counts cutoff={}, PTPRC\ncor={:.2f}'.format(
                    cutoff, cor))
        fig.text(0.5, 0.05, 'CD45 Median Intensity', ha='center')
        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()
#%%
immune_genes = [
    'CD3D','BCL11B','TRAT1','CD3G','ITK','CD2','THEMIS','UBASH3A','TRBC1',
    'CD6','CD28','CD247','TRAC','SIRPG','SH2D1A']
immune_genes = [x for x in immune_genes if x in cts_oc.columns][1:7]
#%%
geneset = immune_genes
cpm_cts = cts_oc_wt/cts_oc_wt.sum()*1e6
cpm_cts = np.log2(cpm_cts+1)
_, axes = plt.subplots(2,3, figsize = (9,6), sharex=True, sharey=True)
axes = axes.ravel()
for ax, gene in zip(axes, geneset):
    cor = np.corrcoef(cpm_cts['CD3D'],cpm_cts[gene])[0,1]
    ax.scatter(cpm_cts['CD3D'], cpm_cts[gene], s=1)
    ax.set_xlabel('CD3D')
    ax.set_title('Correlation:{:.2f}'.format(cor))
    ax.set_ylabel(gene)
#%%
geneset = immune_genes
cpm_cts = cts_oc/cts_oc.sum()*1e6
cpm_cts = np.log2(cpm_cts+1)
_, axes = plt.subplots(2,3, figsize = (9,6), sharex=True, sharey=True)
axes = axes.ravel()
for ax, gene in zip(axes, geneset):
    cor = np.corrcoef(cpm_cts['CD3D'],cpm_cts[gene])[0,1]
    ax.scatter(cpm_cts['CD3D'], cpm_cts[gene], s=1)
    ax.set_xlabel('CD3D')
    ax.set_title('Correlation:{:.2f}'.format(cor))
    ax.set_ylabel(gene)
# %%
# %%

# %%
