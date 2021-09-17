#%%
import pandas as pd
import scanpy as sc
import os
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import minmax_scale
import numpy as np
#%%
os.chdir('/home2/s190548/work_xiao/projects/MOCCA/data/Sprod_ready_data/visium/breast_cancer/')
cts = pd.read_csv('Counts.txt',sep='\t', index_col=0)
spot_meta = pd.read_csv('Spot_metadata.csv', index_col=0)
immune_genes = [
    'CCL13','IL1RN','CXCL9','CCL11','EBI3','IL24','CXCR4','CXCL13','IL2RG',
    'CCL8','CCL5','CCL4','CCL3','TNFRSF17','CCR7','CCL19','CCL18','IL32','NGFR',
    'IL33','CCL22','CCL21','BMP8A','IL16','INHBA','NGF','CXCL10','CXCL11',
    'CXCL12','IL1B','IL2RB','CD27','LTB','IL7R','CCL13','CXCL9','CCL11','CCL22',
    'CCL21','IL24','CXCR4','CXCL13','IL2RG','CXCL10','CXCL11','CXCL12','CCL8',
    'CCL5','CCL4','IL2RB','CCL3','CCR7','CCL19','CCL18','ITK','CCL13','CXCL9',
    'CCL11','CCL22','CCL21','ADCY4','CXCR4','CXCL13','RASGRP2','CXCL10',
    'CXCL11','CXCL12','CCL8','CCL5','CCL4','PLCG2','RAC2','CCL3','CCR7','CCL19',
    'CCL18','JAK3','CD79A','ZAP70','LCK','CD19','TAP1','CD3E','IL2RG','IL7R',
    'CD3D','JAK3','COMP','SELP','IL1B','PECAM1','HBB','HBA2','CD36','HBA1',
    'ACKR1','THBS2','SELE','ITK','MMP2','RHOH','CXCR4','THY1','MMP9','CDH5',
    'CLDN5','CXCL12','PLCG2','RAC2','PECAM1','ESAM','MYL9','ITGA1','CD3E',
    'CD3D','CD2','CD19','IL1B','CD38','CD37','CD36','IL7R','CD34','MS4A1','CD22',
    'ZAP70','CCL13','CXCL12','GADD45B','CCL21','LCK','IL1B','CCL4','PLCG2',
    'LTB','CCL19','DLL4','ZAP70','LCK','IL2RB','FOS','CD247','CD3E','IL2RG',
    'CD3D','JAK3','DUSP4','MAP4K1','PDGFRB','NGFR','FLT1','GADD45B','HSPA6',
    'VEGFC','FOS','RASGRP2','NGF','AREG','CACNA1H','NR4A1','IL1B','KDR','RAC2',
    'TEK','ZAP70','LCK','IL1B','IL2RB','FOS','CD247','CD3E','IL2RG','CD3D',
    'JAK3','CD79B','CD79A','IFITM1','CD19','RAC2','PLCG2','FOS','CD22','PDGFRB',
    'NGFR','PLA2G2D','FLT1','PLA2G2A','VEGFC','RASAL3','RASGRP2','NGF','ZAP70',
    'KDR','PLCG2','RAC2','TEK']

stromal_genes = [
    'DCN','SFRP4','THBS2','CXCL14','COL10A1','ACTG2','SH2D1A','SULF1','BGN',
    'CXCL12','COL1A2','VCAM1','SCUBE2','LUM','CDH5','RAMP3','EDIL3','COL15A1',
    'PLXDC1','COL6A3','COL3A1','COMP','ITM2A','ASPN','LRRC15','COL8A2','CD93',
    'GREM1','LMOD1','KCNJ8','KDR','MFAP5','ISLR','OLFML2B','CILP','COL14A1',
    'EMCN','RARRES2','PDGFRB','COL5A3','TNN','ENPEP','HDC','RSAD2','MMP3',
    'CXCL9','FBLN2','IL1B']

graph = pd.read_csv('/home2/s190548/temp/BC_trial_Detected_graph.txt',sep='\t')
features = pd.read_csv('spot_level_intensity_features.csv',index_col=0)
# %%
adata = sc.AnnData(cts)
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)
sc.pp.scale(adata)
# %%
sc.tl.score_genes(adata,immune_genes,ctrl_size=200,score_name='immune')
sc.tl.score_genes(adata,stromal_genes,ctrl_size=200,score_name='stromal')
scores_pre = adata.obs[['immune','stromal']]
#%%
sc.pp.regress_out(adata,['immune','stromal'])
sc.pp.scale(adata)
#%%
sc.tl.score_genes(adata,immune_genes,ctrl_size=200,score_name='immune')
sc.tl.score_genes(adata,stromal_genes,ctrl_size=200,score_name='stromal')
scores_post = adata.obs[['immune','stromal']]
# %%
_ = plt.figure(figsize=(8,8))
sns.scatterplot(x=spot_meta.X, y=spot_meta.Y, hue=scores_pre, linewidth=0)
# %%
_ = plt.figure(figsize=(8,8))
sns.scatterplot(x=spot_meta.X, y=spot_meta.Y, hue=scores_post, linewidth=0)
# %%
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.leiden(adata)
sc.tl.umap(adata)
#%%
sc.pl.umap(adata, color='leiden')
#%%
_ = plt.figure(figsize=(8,8))
sns.scatterplot(
    x=spot_meta.X, y=spot_meta.Y, hue=adata.obs.leiden,linewidth=0)
# %%
edges = graph.stack()
edges = edges[edges>0.5]
edges = edges.reset_index()
x0 = spot_meta.loc[edges.iloc[:,0],'X'].values
x1 = spot_meta.loc[edges.iloc[:,1],'X'].values
y0 = spot_meta.loc[edges.iloc[:,0],'Y'].values
y1 = spot_meta.loc[edges.iloc[:,1],'Y'].values
dx = x1-x0
dy = y1-y0
#%%
spot_meta['celltype'] = adata.obs.leiden
_, ax = plt.subplots(figsize=(5,5))
r1000 = np.random.choice(len(edges), 1000, replace=False)
solid_spots = edges.iloc[r1000,:2].values.flatten()
transparent_spots = [x for x in spot_meta.index if x not in solid_spots]
solid_meta = spot_meta.loc[solid_spots].sort_values('celltype')
transparent_meta = spot_meta.loc[transparent_spots].sort_values('celltype')
sns.scatterplot(
    x = solid_meta['X'], y = solid_meta['Y'], 
    hue = solid_meta['celltype'],s=12, linewidth=0, ax=ax
    )
sns.scatterplot(
    x = transparent_meta['X'], y = transparent_meta['Y'], 
    hue = transparent_meta['celltype'],s=6, linewidth=0,
    legend=None, alpha=0.5, ax=ax)
for i in r1000:
    ax.arrow(x0[i],y0[i],dx[i],dy[i], head_width=0,color='#595959', alpha=0.5,
    width=0.0005)
plt.legend(bbox_to_anchor=(1,0.5), loc = 'center left', markerscale=2)
ax.set_facecolor('#d1d1d1')
# %%
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from umap import UMAP
# %%
f_umap = UMAP(min_dist=0.25).fit_transform(
    PCA(0.95, whiten=True).fit_transform(features)
)

# %%
_ = plt.figure(figsize=(8,8))
sns.scatterplot(
    x=f_umap[:,0], y=f_umap[:,1], hue=adata.obs.leiden,linewidth=0)
# %%
_spot_meta = spot_meta.copy()
_spot_meta['celltype'] = adata.obs.leiden
_spot_meta.X = f_umap[:,0]
_spot_meta.Y = f_umap[:,1]

edges = graph.stack()
edges = edges[edges>0.5]
edges = edges.reset_index()
x0 = _spot_meta.loc[edges.iloc[:,0],'X'].values
x1 = _spot_meta.loc[edges.iloc[:,1],'X'].values
y0 = _spot_meta.loc[edges.iloc[:,0],'Y'].values
y1 = _spot_meta.loc[edges.iloc[:,1],'Y'].values
dx = x1-x0
dy = y1-y0

_, ax = plt.subplots(figsize=(5,5))
r1000 = np.random.choice(len(edges), 1000, replace=False)
solid_spots = edges.iloc[r1000,:2].values.flatten()
transparent_spots = [x for x in _spot_meta.index if x not in solid_spots]
solid_meta = _spot_meta.loc[solid_spots].sort_values('celltype')
transparent_meta = _spot_meta.loc[transparent_spots].sort_values('celltype')
sns.scatterplot(
    x = solid_meta['X'], y = solid_meta['Y'], 
    hue = solid_meta['celltype'],s=12, linewidth=0, ax=ax
    )
sns.scatterplot(
    x = transparent_meta['X'], y = transparent_meta['Y'], 
    hue = transparent_meta['celltype'],s=6, linewidth=0,
    legend=None, alpha=0.5, ax=ax)
for i in r1000:
    ax.arrow(x0[i],y0[i],dx[i],dy[i], head_width=0,color='#595959', alpha=0.5,
    width=0.0005)
plt.legend(bbox_to_anchor=(1,0.5), loc = 'center left', markerscale=2)
ax.set_facecolor('#d1d1d1')
# %%

edges['sc'] = _spot_meta.loc[edges.iloc[:,0],'celltype'].values
edges['tc'] = _spot_meta.loc[edges.iloc[:,1],'celltype'].values
# %%

# %%
sns.boxplot(x=(edges.sc == edges.tc), y = edges.iloc[:,2])
# %%
ys = []
xs = np.arange(0.05,1.05,0.05)
for t in xs:
    edges = graph.stack()
    edges = edges[edges>=t]
    edges = edges.reset_index()
    edges['sc'] = _spot_meta.loc[edges.iloc[:,0],'celltype'].values
    edges['tc'] = _spot_meta.loc[edges.iloc[:,1],'celltype'].values
    te = (edges.sc == edges.tc).sum()
    ys.append(te/len(edges))
# %%
sns.scatterplot(x = xs, y = ys)
plt.ylabel('Fraction of same cluster edges.')
plt.ylabel('Cutoff')
# %%
