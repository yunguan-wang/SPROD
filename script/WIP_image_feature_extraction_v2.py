'''
This script will extract two intensity features, median and std, as well as six haralick features from image.
It is tested to work with Visium he and if images.
A conda env is available in the project folder. 

# input folder must have the raw tif image and the spot_metadata csv file.
# The spot metadata file MUST have 'X', 'Y', and 'Spor_radius' columns

# example usage:
conda activate /project/shared/xiao_wang/projects/MOCCA/conda_env
python /project/shared/xiao_wang/projects/MOCCA/code/feature_extraction.py \
    /project/shared/xiao_wang/projects/MOCCA/data/Visium/ovarian_cancer_immune \
    if /project/shared/xiao_wang/projects/MOCCA/data/Visium/ovarian_cancer_immune/intermediate
    
'''
import numpy as np
import pandas as pd
import os
from skimage import io, img_as_float32, morphology, exposure
from skimage.feature import greycomatrix, greycoprops
from itertools import product
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale, minmax_scale
from sklearn.cluster import k_means
from skimage.color import (separate_stains, hdx_from_rgb)

def extract_img_features(
    input_path, input_type, output_path, img = None,img_meta = None,
    feature_mask_shape = 'spot'):
    if img_meta is None:
        img_meta = pd.read_csv('Spot_metadata.csv', index_col=0)
    if img is None:
        img_tif = [x for x in os.listdir(input_path) if 'tif' in x][0]
        if input_type == 'if':
            # the indexing is a workaround for the strange Visium if image channels.
            img = io.imread(img_tif)[:,:,[3,1,0]]
        else:
            img = io.imread(img_tif)
        img = img_as_float32(img)
        img = (255*img).astype('uint8')
    # normalize image with color deconv and histogram equalization
    img = separate_stains(img, hdx_from_rgb)
    img = minmax_scale(
    img.reshape(-1,3)).reshape(img.shape)
    img = np.clip(img, 0 ,1)
    img = exposure.equalize_adapthist(img,clip_limit = 0.01)
    img = (255 * img).astype('uint8')
    # Hard coded type of Haralick features and Angles for searching for neighboring pixels
    # hard coded number of angles to be 4, meaning horizontal, vertical and two diagonal directions.
    # extracting block shaped features
    if feature_mask_shape == 'block':
        tmp = img_meta.sort_values(['Row','Col'])
        block_y = int(np.median(tmp.Y.values[2:-1] - tmp.Y.values[1:-2]) // 2)
        tmp = img_meta.sort_values(['Col','Row'])
        block_x = int(np.median(tmp.X.values[2:-1] - tmp.X.values[1:-2]) // 2)
        # block_r = max(block_x,block_y)
        # block_x = block_y = 10
    print('Prossessing {}'.format(input_path))
    feature_set = [
        'contrast','dissimilarity','homogeneity','ASM','energy','correlation']
    text_features = []
    intensity_features = []
    for i in range(img_meta.shape[0]):
        if (i+1)%100 == 0:
            print('Processing {} spot out of {} spots'.format(i+1, img_meta.shape[0]))
        row = img_meta.iloc[i]
        x, y, r = row[['X','Y','Spot_radius']].astype(int)
        if feature_mask_shape == 'spot':
            spot_img = img[x-r:x+r+1,y-r:y+r+1]
            spot_mask = morphology.disk(r)
            # only use the spot, not the bbox
            spot_img = np.einsum('ij,ijk->ijk',spot_mask,spot_img)
        else:
            spot_img = img[x-block_x:x+block_x+1,y-block_y:y+block_y+1]
            spot_mask = np.ones_like(spot_img[:,:,0],dtype='bool')

        # extract texture features
        ith_texture_f = []
        for c in range(img.shape[2]):
            glcm = greycomatrix(
                spot_img[:,:,c], distances=[1], 
                # Angles are arranged in a counter clockwise manner, in radian.
                angles=[0, np.pi/4, np.pi/2, 3*np.pi/4], levels=256,
                                symmetric=True, normed=False)
            glcm = glcm[1:,1:]
            glcm = glcm / np.sum(glcm,axis=(0,1))
            for feature_name in feature_set:
                ith_texture_f += greycoprops(glcm, feature_name)[0].tolist()
        # The first 6 features are intensity features, and the rest are Haralicks.
        text_features.append(ith_texture_f)

        # extract intensity features
        int_low = 0.2
        int_high = 0.8
        int_step = 0.1
        q_bins = np.arange(int_low, int_high, int_step)
        ith_int_f = []
        for c in range(img.shape[2]):
            for t in q_bins:
                ith_int_f.append(np.quantile(spot_img[:,:,c][spot_mask==True],t))
        intensity_features.append(ith_int_f)

    # Naming the features. f stands for channels, A stands for angles.
    # construct texture feature table
    channels = ['f' + str(i) for i in range(img.shape[2])]
    col_names = product(
        channels,
        feature_set, 
        ['A1','A2','A3','A4'])
    col_names = ['_'.join(x) for x in col_names]
    text_features = pd.DataFrame(text_features,index=img_meta.index,columns=col_names)
    # construct intensity feature table
    intensity_features = pd.DataFrame(
        intensity_features,
        index=img_meta.index,
        columns= [
            '_'.join(x) for x in product(
                channels,
                ['{:.1f}'.format(x) for x in q_bins])])
    # spot_features.to_csv(output_path + '/Spot_level_haralick_features.csv')
    return text_features,intensity_features

img_fn = '/project/shared/xiao_wang/projects/MOCCA/data/Visium/LN/V1_Human_Lymph_Node_image.tif'
# img_fn_norm = '/project/shared/xiao_wang/projects/MOCCA/data/Visium/LN/V1_Human_Lymph_Node_image_adaptiveHE.tif'
img_meta = '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/visium/LN/Spot_metadata.csv'

img = io.imread(img_fn)
# img_norm = io.imread(img_fn_norm)
img_meta = pd.read_csv(img_meta, index_col=0)

# f = extract_img_features('','HE','',img, img_meta)
# f_norm = extract_img_features('','HE','',img_norm, img_meta)
# f_pca = PCA(n_components=10).fit_transform(scale(f))
# f_pca_norm = PCA(n_components=10).fit_transform(scale(f_norm))

# sns.displot(f_pca[:,0])
# sns.displot(f_pca_norm[:,0])
# sns.scatterplot(x = f_pca[:,0], y= f_pca[:,1])
# sns.scatterplot(x = f_pca_norm[:,0], y= f_pca_norm[:,1])

# clusters = [str(x) for x in k_means(f_pca,5)[1]]
# clusters_norm = [str(x) for x in k_means(f_pca_norm,5)[1]]
# _ = sns.mpl.pyplot.figure(figsize = (5,5))
# sns.scatterplot(
#     y = img_meta.X.max() - img_meta.X, x = img_meta.Y,
#     hue=clusters, hue_order = sorted(set(clusters)), s=10, markers = ['h'],
#     linewidth=0)

# extracting block features
# img_norm = separate_stains(img, hdx_from_rgb)
# img_norm = minmax_scale(
#     img_norm.reshape(-1,3)).reshape(img_norm.shape)
# img_norm = np.clip(img_norm, 0 ,1)
# img_norm = exposure.equalize_adapthist(img_norm,clip_limit = 0.01)
# img_norm = (255 * img_norm).astype('uint8')

f_t, f_i = extract_img_features(
    '','HE','',img_norm, img_meta,feature_mask_shape='spot')
# valid_cols = [x for x in f_block.columns if 'homo' not in x]
# valid_cols = [x for x in valid_cols if 'corr' not in x]
f_pca = PCA(n_components=10).fit_transform(scale(f_i))
clusters = [str(x) for x in k_means(f_pca,4)[1]]

# _ = sns.mpl.pyplot.figure(figsize = (10,10))
# ax = sns.scatterplot(
#     y = img_meta.X.max() - img_meta.X, x = img_meta.Y,
#     hue=clusters, hue_order = sorted(set(clusters)), markers = ['h'],
#     linewidth=0)
# ax.set_facecolor('grey')

_ = sns.mpl.pyplot.figure(figsize = (32,32))
ax = sns.scatterplot(
    y = img_meta.X, x = img_meta.Y,
    hue=clusters, hue_order = sorted(set(clusters)), markers = ['h'],
    linewidth=0, alpha=0.25, s=250)
ax.set_facecolor('grey')
io.imshow(img, alpha=0.9)
sns.mpl.pyplot.savefig('/home2/s190548/temp/feature_improvement_test.pdf')
sns.mpl.pyplot.close()

f_i.to_csv('/home2/s190548/temp/LN_block_feature.csv')


# center_patch = img[3000:8000,3000:8000]
# center_patch_meta = img_meta[
#     (img_meta.X>=3200)&
#     (img_meta.X<7800)&
#     (img_meta.Y>=3200)&
#     (img_meta.Y<7800)
#     ]
# center_patch_meta.X-=3000
# center_patch_meta.Y-=3000

# center_patch_sep = separate_stains(center_patch, hdx_from_rgb)
# center_patch_sep = minmax_scale(
#     center_patch_sep.reshape(-1,3)).reshape(center_patch_sep.shape)
# center_patch_sep = np.clip(center_patch_sep, 0 ,1)
# center_patch_sep = exposure.equalize_adapthist(center_patch_sep,clip_limit = 0.01)
# center_patch_sep = (255 * center_patch_sep).astype('uint8')

# f_t, f_i = extract_img_features(
#     '','HE','',img = center_patch_sep, img_meta = center_patch_meta,
#     feature_mask_shape='spot')

# # valid_cols = [x for x in f_block.columns if 'homo' not in x]
# # valid_cols = [x for x in valid_cols if 'corr' not in x]
# # f_pca_b = PCA(n_components=10).fit_transform(scale(f_block[valid_cols]))
# # clusters_b = [str(x) for x in k_means(f_pca_b,5)[1]]

# f_pca_b = PCA(n_components=10).fit_transform(scale(f_i))
# clusters_b = [str(x) for x in k_means(f_pca_b,3)[1]]

# _ = sns.mpl.pyplot.figure(figsize = (16,16))
# ax = sns.scatterplot(
#     y = center_patch_meta.X, x = center_patch_meta.Y,
#     hue=clusters_b, hue_order = sorted(set(clusters_b)), markers = ['h'],
#     linewidth=0, alpha=0.25, s=250)
# ax.set_facecolor('grey')
# io.imshow(center_patch[:,:,0], alpha=0.9)

# io.imshow(center_patch_sep[:,:,0],cmap='gray')

