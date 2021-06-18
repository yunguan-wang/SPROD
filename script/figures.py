#%%
from numpy.lib.function_base import percentile
import pandas as pd
import seaborn as sns
import os
import matplotlib.pyplot as plt
import numpy as np
from skimage import morphology, io
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize as cm_norm
import matplotlib.cm as cm
from sklearn.preprocessing import minmax_scale
#%%
sns.set_theme(context='paper',style='white', font_scale=2)
output_folder = '/project/shared/xiao_wang/projects/MOCCA/figures/'
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

def compare_dropouts(*cts_dfs, labels = None, output_fn = None):
    for i, cts in enumerate(cts_dfs):
        if labels is not None:
            label = labels[i]
        else:
            label = ''
        zero_ratio = (cts<=1e-8).sum()/cts.shape[0]
        cpm = np.log2(1+cts.apply(lambda x: 1e4*x/x.sum(),axis=1).fillna(0))
        gene_means = cpm.mean()
        # gene_means = np.log2(1+cts.mean())
        plt.scatter(gene_means, zero_ratio,  s=5, label=label, alpha=0.5)
    plt.legend(bbox_to_anchor = (1, 0.5), loc = 'center left', markerscale=2)
    plt.xlabel('Mean log2(1+CPM)')
    plt.ylabel('Zero ratio')
    if output_fn is not None:
        plt.savefig(output_fn, bbox_inches = 'tight')
        plt.close()

def spot_pseudo_img_gen(spot_meta, render_values, spot_r = 2):
    '''
    Make pseudo image with spots.
    '''
    img_scale = 2*spot_r + 1
    render_values = (
        render_values - render_values.min())/(np.percentile(render_values, 99)-render_values.min())
    render_values = render_values.clip(0 , 1)
    render_values = (255*render_values).astype('uint8')
    pseudo_img = np.zeros(
        (img_scale*(spot_meta.Row.max()+1),img_scale*(spot_meta.Col.max()+1)),
        dtype='uint8')
    r = int(1.5*spot_r)
    spot_mask = morphology.disk(r)
    for i in range(spot_meta.shape[0]):
        idx = spot_meta.index[i]
        row, col = spot_meta.iloc[i][['Row','Col']].astype(int)
        pseudo_img[
            img_scale*row-r:img_scale*row+r+1,
            img_scale*col-r:img_scale*col+r+1
            ][spot_mask==1] = np.max([render_values[idx],10])
    # Scale the pseudo image so that each channel range from 0-255
    return pseudo_img

def merge_pseudo_img(
    img_prefix, r_img = None, g_img = None, b_img = None,
    channel_names=[None,None,None], merge_mode = [0,1,2]):
    crit = [x is None for x in [r_img, g_img, b_img]]
    if set(crit) == {None}:
        return
    else:
        x, y = [r_img, g_img, b_img][np.argmin(crit)].shape
        pseudo_img = np.zeros([x, y, 3], dtype='uint8')
        for i, img in enumerate([r_img, g_img, b_img]):
            _tmp_img = np.zeros([x, y, 3], dtype='uint8')
            _tmp_img[:,:,i] = img
            _tmp_img[img==0] = 255
            io.imsave(img_prefix + channel_names[i] + '.tif',_tmp_img)
            if i not in merge_mode:
                continue
            pseudo_img[:,:,i] = img
        pseudo_img[(pseudo_img==0).sum(axis=2)==3] = 255
        io.imsave(
            img_prefix + 'merged.tif',pseudo_img)
        return pseudo_img

def marker_overlay(
    gene, spot_meta, cts, img, img_type='HE', output=None, no_raw_img = False):
    x_min, x_max = spot_meta.X.min()//100 * 100, (1+spot_meta.X.max()//100) * 100
    y_min, y_max = spot_meta.Y.min()//100 * 100, (1+spot_meta.Y.max()//100) * 100
    gene_exp_mask = np.zeros(img.shape[:2])
    max_expression = cts[gene].max()
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
            gene_exp_level = np.max([cts.loc[id,gene], 0.05*max_expression])
        gene_exp_mask[x-r:x+r+1,y-r:y+r+1] = gene_exp_level
        gene_exp_mask[x-r:x+r+1,y-r:y+r+1] = gene_exp_mask[x-r:x+r+1,y-r:y+r+1]*bbox_mask
    gene_exp_mask = 255*minmax_scale(gene_exp_mask.flatten()).reshape(gene_exp_mask.shape)
    gene_exp_mask = gene_exp_mask.astype('uint8')
    tmp_array = np.zeros(img.shape, 'uint8')
    tmp_array[:,:,1] = gene_exp_mask
    gene_exp_mask = tmp_array
    if img_type == 'IF':
        exp_mask = np.zeros(
            (gene_exp_mask.shape[0],gene_exp_mask.shape[1],3), dtype='uint8')
        exp_mask[:,:,1] = gene_exp_mask
        gene_exp_mask = exp_mask.copy()
    _, ax = plt.subplots(figsize=(36,36))
    # merged_img = img[x_min:x_max, y_min:y_max].copy()
    # spots_mask = gene_exp_mask[x_min:x_max, y_min:y_max,1]!=0
    # merged_img[:,:,1][spots_mask] = np.add(
    #     0.7 * merged_img[:,:,1][spots_mask],
    #     0.3 * gene_exp_mask[x_min:x_max, y_min:y_max,1][spots_mask]
    # ).astype('uint8')
    if no_raw_img:
        merged_img = gene_exp_mask
    else:
        merged_img = (
            0.7*img[x_min:x_max, y_min:y_max] + 0.3*gene_exp_mask[x_min:x_max, y_min:y_max]
            ).astype('uint8')
    io.imshow(merged_img, ax=ax)
    plt.axis('off')
    if output is not None:
        plt.savefig(output)
        plt.close()
    return gene_exp_mask


def rgb_to_cmyk(rgba):
    r,g,b = rgba[:3]
    if (r, g, b) == (0, 0, 0):
        # black
        return 0, 0, 0, 1
    # rgb [0,255] -> cmy [0,1]
    c = 1 - r
    m = 1 - g
    y = 1 - b

    # extract out k [0, 1]
    k = min(c, m, y)
    c = (c - k) / (1 - k)
    m = (m - k) / (1 - k)
    y = (y - k) / (1 - k)
    # rescale to the range [0,CMYK_SCALE]
    return (c, m, y, k)

def cmyk_to_rgb(cmyk, a = 1):
    c, m, y, k = cmyk
    r = (1-c) * (1-k)
    g = (1-m) * (1-k)
    b = (1-y) *(1-k)
    return (r, g, b, a)

#%%
# ==============================================================================
# ==============================================================================
# ==============================================================================
# figure 1
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

########

# Slideseq variability check.
#%%
tc = slideseq_counts.sum(axis=1)
cpm = slideseq_counts.apply(lambda x: 1e3*x/x.sum(),axis=1)
n_genes = (slideseq_counts>0).sum(axis=1)
gene_vs_counts_df = pd.DataFrame(n_genes, columns = ['Number of expressed genes'])
tc_cuts = pd.qcut(tc, 4, labels=['Q1','Q2','Q3','Q4']).astype(str)
gene_vs_counts_df['Total read counts'] = tc
gene_vs_counts_df['Total read bins'] = tc_cuts.values
gene_vs_counts_df['Mean of non-zero CPK'] = cpm.apply(
    lambda x: x[x>0].mean(), axis=1
)
#%%
cpm_probs = cpm.stack().reset_index()
cpm_probs = cpm_probs.iloc[:,2]
cpm_probs = cpm_probs[cpm_probs>0]

cpm_probs_single = cpm.apply(
        lambda x: pd.Series(
            [x[x>=i].shape[0]/x[x>0].shape[0] for i in np.arange(0.5,10.5,0.5)],
            index = np.arange(0.5,10.5,0.5)
            ),
        axis=1
        )

cpm_probs_single = cpm_probs_single.stack().reset_index()
cpm_probs_single = cpm_probs_single.set_index('level_0')
cpm_probs_single.loc[:,'Total_counts_bin'] = tc_cuts[cpm_probs_single.index].values
cpm_probs_single.columns = [
    'CPK', 'Ratio of read count above / Total non zero reads', 'Total_counts_bin']
#%%
_ = plt.figure(figsize=(4,4))
sns.relplot(data=gene_vs_counts_df,
    x='Number of expressed genes', y = 'Mean of non-zero CPK',
    col = 'Total read bins', col_order=['Q1','Q2','Q3','Q4'], col_wrap=2,
    s=10, linewidth=0, facet_kws = {'sharex' : False, 'sharey' : False})
plt.savefig(
    '/project/shared/xiao_wang/projects/MOCCA/figures/Slideseq_Puck_190921_19_CPK_inflation.pdf',
    bbox_inches = 'tight')

#%%
_ = plt.figure(figsize=(4,4))
sns.lineplot(
    data=cpm_probs_single,
    x='CPK', y = 'Ratio of read count above / Total non zero reads',
    hue = 'Total_counts_bin', hue_order=['Q1','Q2','Q3','Q4'])
plt.legend(bbox_to_anchor = (1,0.5), loc = 'center left')
plt.savefig(
    '/project/shared/xiao_wang/projects/MOCCA/figures/Slideseq_Puck_190921_19_variability.pdf',
    bbox_inches = 'tight')

########
# Visium Ovarian cancer dropouts figure
#%%
# Visium Ovarian cancer dataset
oc_im_cts = pd.read_csv(
    '/project/shared/xiao_wang/projects//MOCCA/data/Sprod_ready_data/visium/ovarian_cancer_immune/Counts.txt',
    sep='\t', index_col=0)
oc_wt_cts = pd.read_csv(
    '/project/shared/xiao_wang/projects//MOCCA/data/Sprod_ready_data/visium/ovarian_cancer/Counts.txt',
    sep='\t', index_col=0)
img_features = pd.read_csv(
    '/project/shared/xiao_wang/projects//MOCCA/data/Sprod_ready_data/visium/ovarian_cancer_immune/Spot_level_haralick_features.csv',
    index_col=0)
spot_meta = pd.read_csv(
    '/project/shared/xiao_wang/projects//MOCCA/data/Sprod_ready_data/visium/ovarian_cancer_immune/Spot_metadata.csv',
    index_col=0)

high_qual_cells = oc_im_cts.index[oc_im_cts.sum(axis=1)>=oc_im_cts.sum(axis=1).median()]

oc_im_cts = np.log2(1+oc_im_cts.apply(lambda x: 1e4*x/x.sum(),axis=1).fillna(0))
oc_wt_cts = np.log2(1+oc_wt_cts.apply(lambda x: 1e4*x/x.sum(),axis=1).fillna(0))

#%%
scatterplot_data = pd.concat([oc_im_cts['PTPRC'],img_features['f0_median']], axis=1)
scatterplot_data = scatterplot_data.append(
    pd.concat([oc_wt_cts['PTPRC'],img_features['f0_median']], axis=1)
)
scatterplot_data.columns = ['PTPRC log2(CPM+1)', 'CD45']
scatterplot_data['Dataset'] = ['Targeted']*len(oc_wt_cts) + ['Whole transcriptome']*len(oc_wt_cts)

sns.scatterplot(
    data = scatterplot_data,
    x='PTPRC log2(CPM+1)', y='CD45', hue = 'Dataset', linewidth=0, s=10)
plt.legend(bbox_to_anchor=(1,0.5), loc='center left')
plt.savefig(
    '/project/shared/xiao_wang/projects/MOCCA/figures/Visium_OC_CD45_scatterplot.pdf',
    bbox_inches='tight', alpha=0.5)
plt.close()

sns.scatterplot(
    data = scatterplot_data.loc[high_qual_cells],
    x='PTPRC log2(CPM+1)', y='CD45', hue = 'Dataset', linewidth=0, s=10)
plt.legend(bbox_to_anchor=(1,0.5), loc='center left')
plt.savefig(
    '/project/shared/xiao_wang/projects/MOCCA/figures/Visium_OC_CD45_scatterplot_high_quality_spots.pdf',
    bbox_inches='tight', alpha=0.5)
plt.close()
#%%
# Merged plot
for cells, qual in zip(
    [spot_meta.index.tolist(), high_qual_cells],
    ['', 'HighQualitySpots']
    ):
    y_vector = spot_meta.loc[cells,'Row']
    x_vector = spot_meta.loc[cells,'Col']
    y_vector = y_vector.max() - y_vector + 1
    x_vector = x_vector.max() - x_vector + 1

    for color_set, fig_name, c in zip(
        [
            img_features.loc[cells, 'f0_median'],
            oc_im_cts.loc[cells,'PTPRC'],
            oc_wt_cts.loc[cells,'PTPRC']
            ],
        [
            'CD45', 'PTPRC Targetd','PTPRC WholeTranscriptome'
            ],
        [
            'red','green','cyan'
            ]
        ):
        _, ax = plt.subplots(figsize=(9,9))
        sns.scatterplot(
            x = x_vector, y = y_vector, hue=color_set, legend=None, s=50, marker='H',
            palette = sns.light_palette(c, as_cmap=True), linewidth=0)
        ax.set_facecolor('#bfbfbf')
        plt.savefig(
            '/project/shared/xiao_wang/projects/MOCCA/figures/Visium_OC_pseudo_img_{}_{}.pdf'.format(qual, fig_name),
            bbox_inches = 'tight')
        plt.close()

    cd45_norm = cm_norm(
        vmin=np.percentile(img_features.loc[cells, 'f0_median'],5),
        vmax=np.percentile(img_features.loc[cells, 'f0_median'],95),
        clip=True)
    cd45_mapper = cm.ScalarMappable(
        norm=cd45_norm, cmap=sns.light_palette("red", as_cmap=True))
    ptprc_norm = cm_norm(
        vmin=oc_im_cts.loc[cells,'PTPRC'].min(),
        vmax=oc_im_cts.loc[cells,'PTPRC'].max(),clip=True)
    ptprc_mapper = cm.ScalarMappable(
        norm=ptprc_norm, cmap=sns.light_palette("green", as_cmap=True))
    merged_colors = []
    for v in cells:
        v_cd45 = img_features.loc[v, 'f0_median']
        v_ptprc = oc_im_cts.loc[v,'PTPRC']
        c1 = cd45_mapper.to_rgba(v_cd45)
        c2 = ptprc_mapper.to_rgba(v_ptprc)
        if v_ptprc > 0:
            c1_cmyk = rgb_to_cmyk(c1)
            c2_cmyk = rgb_to_cmyk(c2)
            c3 = tuple((c1_cmyk[i]*c2_cmyk[i]) for i in range(4))
            c3 = cmyk_to_rgb(c3)
        else:
            c3 = c1
        merged_colors.append(c3)

    _, ax = plt.subplots(figsize=(9,9))
    ax.scatter(
        x_vector, y_vector, c = merged_colors, linewidth=0, s=50, marker='H')
    ax.set_facecolor('#bfbfbf')
    plt.savefig(
        '/project/shared/xiao_wang/projects/MOCCA/figures/Visium_OC_pseudo_img_{}_merged.pdf'.format(qual),
        bbox_inches = 'tight')  
    plt.close()

########
# Visium LN expression noise figure
#%%
raw_cts = pd.read_csv(
    '/project/shared/xiao_wang/projects//MOCCA/data/Sprod_ready_data/visium/LN/Counts.txt',
    sep='\t', index_col=0)
img_features = pd.read_csv(
    '/project/shared/xiao_wang/projects//MOCCA/data/Sprod_ready_data/visium/LN/Spot_level_haralick_features.csv',
    index_col=0)
spot_meta = pd.read_csv(
    '/project/shared/xiao_wang/projects//MOCCA/data/Sprod_ready_data/visium/LN/Spot_metadata.csv',
    index_col=0)
img_fn = '/project/shared/xiao_wang/projects/MOCCA/data/Visium/LN/V1_Human_Lymph_Node_image.tif'
img = io.imread(img_fn)
ln_cpm = np.log2(1+raw_cts.apply(lambda x: 1e4*x/x.sum(),axis=1).fillna(0))
#%%
igd = marker_overlay(
    'IGHD', spot_meta, ln_cpm, img,
    output = output_folder + 'Visium_LN_IGHD.pdf')
ms4a1 = marker_overlay(
    'MS4A1', spot_meta, ln_cpm, img,
    output = output_folder + 'Visium_LN_MS4A1.pdf')
cd1c = marker_overlay(
    'CD1C', spot_meta, ln_cpm, img,
    output = output_folder + 'Visium_LN_CD1C.pdf')
cd3 = marker_overlay(
    'CD3D', spot_meta, ln_cpm, img,
    output = output_folder + 'Visium_LN_CD3D.pdf')

#%%
# ==============================================================================
# ==============================================================================
# ==============================================================================
# figure 2
# Grid search on simulation data plot.
input_path = '/home2/s190548/work_xiao/projects/MOCCA/data/simulation'
E = pd.read_csv(os.path.join(input_path, 'Es_5k.txt'), sep='\t')
for noise_level in range(1,6):
    E_noised = pd.read_csv(
        os.path.join(input_path, 'Es_0.{}n_5k.txt'.format(noise_level)), sep='\t')
    total_sae = abs(E_noised-E).sum().sum()
    grid_search_res = pd.read_csv(
        os.path.join(input_path,'../grid_search_0.{}n_clean_F.csv'.format(noise_level))
        )
    y = sorted(1- grid_search_res.SAE_vs_True_E/total_sae)
    if len(y)>120:
        y = sorted(np.random.choice(y, 120, replace=False))
    x = np.arange(len(y))
    plt.plot(x, y, label = 'Noise level = 0.{}'.format(noise_level))
plt.legend(loc='center left', bbox_to_anchor = (1,0.5))
plt.xlabel('Grid search runs')
plt.ylabel('Percent of noise removed')
plt.savefig(
    output_folder + 'Grid search on different noise levels.pdf', bbox_inches = 'tight')
#%%
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

input_path = '/home2/s190548/work_xiao/projects/MOCCA/data/grid_search_0.5n_clean_F'
E = pd.read_csv(
    '/home2/s190548/work_xiao/projects/MOCCA/data/simulation/Es_5k.txt', sep='\t')
E_noised = pd.read_csv(
    '/home2/s190548/work_xiao/projects/MOCCA/data/simulation/Es_0.5n_5k.txt', sep='\t')
total_sae = abs(E_noised-E).sum().sum()
r_vec = [x.split('_')[1] for x in os.listdir(input_path)]
k_vec = [x.split('_')[3] for x in os.listdir(input_path)]
lambda_vec = [x.split('_')[5] for x in os.listdir(input_path)]
percent_denoised = []
for folder in os.listdir(input_path):
    E_denoised = pd.read_csv(
        os.path.join(input_path, folder, 'Denoised_matrix.txt'), sep='\t')
    sae = abs(E_denoised-E).sum().sum()
    percent_denoised.append(1- sae/total_sae)
plot_data = pd.DataFrame([r_vec, k_vec, lambda_vec,percent_denoised]).transpose()
plot_data.columns = ['R','K','Lambda','Denoised']
plot_data = plot_data.astype(float)

# R vs Lambda
for z in plot_data.K.unique():
    _plot_data = plot_data[plot_data.K==z]
    X = sorted(_plot_data.R.astype(float))
    Y = sorted(_plot_data.Lambda.values.astype(float))
    X, Y = np.meshgrid(X, Y)

    Z = np.zeros(shape=(X.shape))
    for i in range(len(X)):
        for j in range(len(X)):
            Z[i,j] = _plot_data.loc[
                (_plot_data.R==X[i,j])&(_plot_data.Lambda==Y[i,j]), 'Denoised']

    # Plot the surface.
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(8,10))
    surf = ax.plot_surface(
        X, Y, Z, alpha=0.5, cmap = 'GnBu', linewidth=0, antialiased=False)
    xmax, ymax, zmax = _plot_data.sort_values('Denoised', ascending=False).iloc[0,[0,2,3]]
    ax.plot([0, xmax], [0.9, ymax], [zmax + 0.005, zmax], lw=3, color="r")
    _ = ax.text(
        0, 0.9, zmax + 0.005, 
        'R = {}, Lambda = {}, % removed noise = {:.1f}'.format(xmax, ymax, 100*zmax),
        fontsize=10
        )
    ax.set_xlabel('\nR', linespacing = 3)
    ax.set_ylabel('\nLambda', linespacing = 3)
    ax.set_zlabel('\nPercentage of removed noise', linespacing = 3)
    ax.set_title('K = {}'.format(int(z)))
    plt.savefig(
        output_folder + 'Grid_search_3d_K={}.pdf'.format(z), bbox_inches = 'tight')

for z in plot_data.Lambda.unique():
    _plot_data = plot_data[plot_data.Lambda==z]
    X = sorted(_plot_data.R.astype(float))
    Y = sorted(_plot_data.K.values.astype(float))
    X, Y = np.meshgrid(X, Y)

    Z = np.zeros(shape=(X.shape))
    for i in range(len(X)):
        for j in range(len(X)):
            Z[i,j] = _plot_data.loc[
                (_plot_data.R==X[i,j])&(_plot_data.K==Y[i,j]), 'Denoised']

    # Plot the surface.
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(8,10))
    surf = ax.plot_surface(
        X, Y, Z, alpha=0.5, cmap = 'GnBu', linewidth=0, antialiased=False)
    xmax, ymax, zmax = _plot_data.sort_values('Denoised', ascending=False).iloc[0,[0,1,3]]
    ax.plot([0, xmax], [11, ymax], [zmax + 0.005, zmax], lw=3, color="r")
    _ = ax.text(
        0, 11, zmax + 0.005, 
        'R = {}, K = {}, % removed noise = {:.1f}'.format(xmax, ymax, 100*zmax),
        fontsize=10
        )
    ax.set_xlabel('\nR', linespacing = 3)
    ax.set_ylabel('\nK', linespacing = 3)
    ax.set_zlabel('\nPercentage of removed noise', linespacing = 3)
    ax.set_title('Lambda = {}'.format(z))
    plt.savefig(
        output_folder + 'Grid_search_3d_Lambda={}.pdf'.format(z), bbox_inches = 'tight')
#%%
########
# Visualize graph
graph = pd.read_csv(
    '/project/shared/xiao_wang/projects/MOCCA/data/grid_search_0.5n_clean_F/3.0_0.08_500.0_3.0_1.0_0.625_1.0_0.001/Detected_graph.txt',
    sep='\t')
spot_meta = pd.read_csv(
    '/project/shared/xiao_wang/projects/MOCCA/data/simulation/C_simu_5k.csv',
    index_col=0)
#%%
edges = graph.stack()
edges = edges[edges>0]
edges = edges.reset_index()
x0 = spot_meta.loc[edges.iloc[:,0],'X'].values
x1 = spot_meta.loc[edges.iloc[:,1],'X'].values
y0 = spot_meta.loc[edges.iloc[:,0],'Y'].values
y1 = spot_meta.loc[edges.iloc[:,1],'Y'].values
dx = x1-x0
dy = y1-y0
#%%
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
plt.savefig(output_folder + 'Simulation E_0.5n graph.pdf', bbox_inches = 'tight')
#%%
# ==============================================================================
# ==============================================================================
# ==============================================================================
# Figure 3
slideseq_raw = '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/Counts.hdf'
slideseq_denoised = '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/denoised_cts.h5df'
shuffled_denoised = '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/coord_shuffled_denoised_cts.hdf'
saver_cts = '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/saver_cts.h5df'
raw_cts = pd.read_hdf(slideseq_raw)
denoised_cts = pd.read_hdf(slideseq_denoised)
saver_cts = pd.read_hdf(saver_cts)
sh_cts = pd.read_hdf(shuffled_denoised)
raw_cts = raw_cts.set_index(raw_cts.columns[0])
cm_barcodes = list(
    set(raw_cts.index) & set(denoised_cts.index) & set(saver_cts.index & set(sh_cts.index))
    )
raw_cts, denoised_cts, saver_cts, sh_cts = (
    raw_cts.loc[cm_barcodes], denoised_cts.loc[cm_barcodes],
    saver_cts.loc[cm_barcodes], sh_cts.loc[cm_barcodes]
)
# %%
raw_cpm = raw_cts.apply(lambda x: 1e3*x/x.sum(),axis=1).fillna(0)
denoised_cpm = denoised_cts.apply(lambda x: 1e3*x/x.sum(),axis=1).fillna(0)
saver_cpm = saver_cts.apply(lambda x: 1e3*x/x.sum(),axis=1).fillna(0)
sh_cpm = sh_cts.apply(lambda x: 1e3*x/x.sum(),axis=1).fillna(0)
tc = raw_cts.sum(axis=1)
tc_cuts = pd.qcut(tc, 4, labels=['Q1','Q2','Q3','Q4']).astype(str)
#%%
# ========
# dropoff comparison plot
df_zeros_r = pd.DataFrame()
for v,l in zip(
    [raw_cpm, saver_cpm,sh_cpm,denoised_cpm],
    ['Raw','Saver','Scrambled','Sprod']
    ):
    _zr = (v<=1e-8).sum()/v.shape[0]
    _gb = v.mean()
    _gb = pd.qcut(
        _gb, 8,
        labels = ['Q1','Q2','Q3','Q4','Q5','Q6','Q7','Q8'])
    _tmp = pd.DataFrame(_zr, columns = ['ZR'])
    _tmp['cpm_bin'] = _gb.values
    _tmp['Dataset'] = l
    df_zeros_r = df_zeros_r.append(_tmp)

sns.boxplot(
    data = df_zeros_r, x = 'cpm_bin', y = 'ZR', hue='Dataset', palette = 'rainbow')
plt.legend(bbox_to_anchor = (1,0.5), loc = 'center left')
plt.xlabel('Gene expression bins')
plt.ylabel('Zero rations')
plt.savefig(output_folder + 'Dropout_comparison_box.pdf', bbox_inches = 'tight')
# %%
cpm_bins = [1e-4, 1e-3, 1e-2, 1e-1, 1, 2, 3, 4, 5, 10, 100, 1000]
cpm_bins = [np.log2(1+x) for x in cpm_bins]
cpm_bin_labels = ['{:.1e}'.format(x) for x in cpm_bins]

raw_cpm_log = np.log2(1+raw_cpm)
raw_cpm_single = raw_cpm_log.apply(
        lambda x: pd.Series(
            [(x>=i).sum()/(x>0).sum() for i in cpm_bins],
            index = cpm_bin_labels
            ),
        axis=1
        )
raw_cpm_single = raw_cpm_single.stack().reset_index()
raw_cpm_single = raw_cpm_single.set_index('level_0')
raw_cpm_single.loc[:,'Total_counts_bin'] = tc_cuts[raw_cpm_single.index].values
raw_cpm_single.columns = [
    'CPK', 'sum(x>i)/sum(x>0)', 'Total_counts_bin']

denoised_cpm_log = np.log2(1+denoised_cpm)
denoised_cpm_single = denoised_cpm_log.apply(
        lambda x: pd.Series(
            [(x>=i).sum()/(x>0).sum() for i in cpm_bins],
            index = cpm_bin_labels
            ),
        axis=1
        )
denoised_cpm_single = denoised_cpm_single.stack().reset_index()
denoised_cpm_single = denoised_cpm_single.set_index('level_0')
denoised_cpm_single.loc[:,'Total_counts_bin'] = tc_cuts[denoised_cpm_single.index].values
denoised_cpm_single.columns = [
    'CPK', 'sum(x>i)/sum(x>0)', 'Total_counts_bin']

saver_cpm_log = np.log2(1+saver_cpm)
saver_cpm_single = saver_cpm_log.apply(
        lambda x: pd.Series(
            [(x>=i).sum()/(x>0).sum() for i in cpm_bins],
            index = cpm_bin_labels
            ),
        axis=1
        )
saver_cpm_single = saver_cpm_single.stack().reset_index()
saver_cpm_single = saver_cpm_single.set_index('level_0')
saver_cpm_single.loc[:,'Total_counts_bin'] = tc_cuts[saver_cpm_single.index].values
saver_cpm_single.columns = [
    'CPK', 'sum(x>i)/sum(x>0)', 'Total_counts_bin']
# %%
_ = plt.figure(figsize=(4,4))
sns.lineplot(
    data=raw_cpm_single,
    x='CPK', y = 'sum(x>i)/sum(x>0)',
    hue = 'Total_counts_bin', hue_order=['Q1','Q2','Q3','Q4'])
plt.legend(bbox_to_anchor = (1,0.5), loc = 'center left')
plt.xlabel('log2(1+CPK)')
plt.xticks(rotation = 90)
plt.ylim((0,1.05))
plt.savefig(
    '/project/shared/xiao_wang/projects/MOCCA/figures/figure3a_variability_raw.pdf',
    bbox_inches = 'tight')

_ = plt.figure(figsize=(4,4))
sns.lineplot(
    data=denoised_cpm_single,
    x='CPK', y = 'sum(x>i)/sum(x>0)',
    hue = 'Total_counts_bin', hue_order=['Q1','Q2','Q3','Q4'])
plt.legend(bbox_to_anchor = (1,0.5), loc = 'center left')
plt.xlabel('log2(1+CPK)')
plt.xticks(rotation = 90)
plt.ylim((0,1.05))
plt.savefig(
    '/project/shared/xiao_wang/projects/MOCCA/figures/figure3a_variability_denoised.pdf',
    bbox_inches = 'tight')

_ = plt.figure(figsize=(4,4))
sns.lineplot(
    data=saver_cpm_single,
    x='CPK', y = 'sum(x>i)/sum(x>0)',
    hue = 'Total_counts_bin', hue_order=['Q1','Q2','Q3','Q4'])
plt.legend(bbox_to_anchor = (1,0.5), loc = 'center left')
plt.xlabel('log2(1+CPK)')
plt.xticks(rotation = 90)
plt.ylim((0,1.05))
plt.savefig(
    '/project/shared/xiao_wang/projects/MOCCA/figures/figure3a_variability_saver.pdf',
    bbox_inches = 'tight')
#%%
# Figure 3d
raw_cts = pd.read_csv(
    '/project/shared/xiao_wang/projects//MOCCA/data/Sprod_ready_data/visium/LN/Counts.txt',
    sep='\t', index_col=0)
spot_meta = pd.read_csv(
    '/project/shared/xiao_wang/projects//MOCCA/data/Sprod_ready_data/visium/LN/Spot_metadata.csv',
    index_col=0)
img_features = pd.read_csv(
    '/project/shared/xiao_wang/projects//MOCCA/data/Sprod_ready_data/visium/LN/Spot_level_haralick_features.csv',
    index_col=0)
weights = pd.read_csv(
    '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/visium/LN/denoised/W.txt',
    sep=' ', header=None)
img_fn = '/project/shared/xiao_wang/projects/MOCCA/data/Visium/LN/V1_Human_Lymph_Node_image.tif'
img = io.imread(img_fn)
ln_cpm = np.log2(1+raw_cts.apply(lambda x: 1e4*x/x.sum(),axis=1).fillna(0))
denoised_cpm = np.matmul(weights.values,ln_cpm.values)
denoised_cpm = pd.DataFrame(denoised_cpm, index = ln_cpm.index, columns = ln_cpm.columns)
#%%
for gene in ['IGHD','MS4A1','CD1C', 'CD3D']:
    for prefix in ['CPM','Sprod']:
        if prefix == 'CPM':
            cts = ln_cpm[[gene]].copy()
        else:
            cts = denoised_cpm[[gene]].copy()
        _ = marker_overlay(
            gene, spot_meta, cts, img,
            output = '/home2/s190548/temp/Visium_LN_Tao_{}_{}.pdf'.format(gene, prefix))
#%%
for gene in ['IGHD','MS4A1','CD1C', 'CD3D']:
    for prefix in ['raw', 'CPM','Sprod','Scramble','Saver','scI']:
        cts = ln_cts[[prefix + '.' + gene]].copy()
        if prefix in ['CPM','Sprod','Scramble']:
            cts = np.log2(1+(np.exp(cts)-1)/100)
        _ = marker_overlay(
            prefix + '.' + gene, spot_meta, cts, img,
            output = output_folder + 'Visium_LN_{}_{}.pdf'.format(gene, prefix))
#%%
for f in img_features.columns:
    cts = img_features.loc[ln_cts.index,[f]]
    _ = marker_overlay(
        f, spot_meta, cts, img,no_raw_img=True,
        output = output_folder + '/featuremap/Visium_LN_{}.pdf'.format(f))
#%%
for gene in ['IGHD','MS4A1','CD1C', 'CD3D']:
    for prefix in ['raw', 'CPM','Sprod','Scramble','Saver','scI']:
        cts = ln_cts[[prefix + '.' + gene]].copy()
        if prefix in ['CPM','Sprod','Scramble']:
            cts = np.log2(1+(np.exp(cts)-1)/100)
        _ = marker_overlay(
            prefix + '.' + gene, spot_meta, cts, img, no_raw_img=True,
            output = output_folder + 'Visium_LN_marker_only_{}_{}.pdf'.format(gene, prefix))
# %%
# ==============================================================================
# ==============================================================================
# ==============================================================================
# some random codes
cpm_single = raw_cpm_single.append(denoised_cpm_single)
cpm_single['dataset'] = ['Raw']*len(raw_cpm_single) + ['Denoised']*len(raw_cpm_single)
cpm_single.to_csv(
    '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/cpm_single.csv')
# %%
patch_folder = '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/patches'
spot_metas = [x for x in os.listdir(patch_folder) if 'meta' in x]
for fn in spot_metas:
    fn = os.path.join(patch_folder, fn)
    meta = pd.read_csv(fn, index_col=0)
    ori_idx = meta.index.copy()
    rand_idx = meta.index.values.copy()
    np.random.shuffle(rand_idx)
    meta.index = rand_idx
    meta = meta.loc[ori_idx]
    meta.to_csv(fn.replace('.csv','.coord_shuffled.csv'))
# %%
from scipy.stats import entropy
e_raw = entropy(raw_cpm, axis=1)
e_saver = entropy(saver_cpm, axis=1)
e_sh = entropy(sh_cpm, axis=1)
e_denoised = entropy(denoised_cpm, axis=1)

for v,l in zip([e_raw, e_saver,e_sh,e_denoised],['Raw','Saver','Scrambled','Sprod']):
    _ = plt.hist(v, bins=500, alpha=0.5, label=l)
plt.legend()
# %%
for v,l in zip([e_sh,e_denoised],['Scrambled','Sprod']):
    _ = plt.hist(v, bins=500, alpha=0.5, label=l)
plt.legend()
# %%
