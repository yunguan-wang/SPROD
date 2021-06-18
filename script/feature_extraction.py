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
import sys
from skimage import io, img_as_float32, morphology
from skimage.feature import greycomatrix, greycoprops
from itertools import product

def extract_img_features(input_path, input_type, output_path):
    img_tif = [x for x in os.listdir(input_path) if 'tif' in x][0]
    img_meta = pd.read_csv('Spot_metadata.csv', index_col=0)
    if input_type == 'if':
        # the indexing is a workaround for the strange Visium if image channels.
        img = io.imread(img_tif)[:,:,[3,1,0]]
    else:
        img = io.imread(img_tif)
    img = img_as_float32(img)
    img = (255*img).astype('uint8')

    # Hard coded type of Haralick features and Angles for searching for neighboring pixels
    feature_set = ['contrast','dissimilarity','homogeneity','ASM','energy','correlation']
    # hard coded number of angles to be 4, meaning horizontal, vertical and two diagonal directions.
    spot_features = np.zeros([img_meta.shape[0],img.shape[2]*(2+4*len(feature_set))])
    print('Prossessing {}'.format(input_path))
    for i in range(img_meta.shape[0]):
        if (i+1)%100 == 0:
            print('Processing {} spot out of {} spots'.format(i+1, img_meta.shape[0]))
        row = img_meta.iloc[i]
        x, y, r = row[['X','Y','Spot_radius']].astype(int)
        spot_img = img[x-r:x+r+1,y-r:y+r+1]
        spot_mask = morphology.disk(r)
        # only use the spot, not the bbox
        spot_img = np.einsum('ij,ijk->ijk',spot_mask,spot_img)
        haralick_features = []
        for c in range(img.shape[2]):
            glcm = greycomatrix(
                spot_img[:,:,c], distances=[1], 
                # Angles are arranged in a counter clockwise manner, in radian.
                angles=[0, np.pi/4, np.pi/2, 3*np.pi/4], levels=256,
                                symmetric=True, normed=False)
            glcm = glcm[1:,1:]
            glcm = glcm / np.sum(glcm,axis=(0,1))
            spot_features[i, 2*c] = np.median(spot_img[:,:,c][spot_mask==True]) 
            spot_features[i, 2*c + 1] = np.std(spot_img[:,:,c][spot_mask==True])
            for feature_name in feature_set:
                haralick_features += greycoprops(glcm, feature_name)[0].tolist()
        # The first 6 features are intensity features, and the rest are Haralicks.
        spot_features[i,6:] = haralick_features

    # Naming the features. f stands for channels, A stands for angles.
    col_names = product(['f0','f1','f2'],feature_set, ['A1','A2','A3','A4'])
    col_names = ['_'.join(x) for x in col_names]
    col_names = ['_'.join(x) for x in product(['f0','f1','f2'],['median','std'])] + col_names
    spot_features = pd.DataFrame(spot_features,index=img_meta.index,columns=col_names)
    spot_features.to_csv(output_path + '/Spot_level_haralick_features.csv')

if __name__ == '__main__':
    # Todo: decide if checking metadata function should be added, or force the user to provide the correct format.
    input_path = sys.argv[1]
    input_type = sys.argv[2] # {'he' or 'if'}

    try:
        output_path = sys.argv[3]
    except IndexError:
        output_path = './intermediate'

    os.chdir(input_path)
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    print('Processing {}'.format(input_path))
    if os.path.exists(output_path + '/Spot_level_haralick_features.csv'):
        print('Image feature already extracted. Stop operation.')
    else:
        extract_img_features(input_path, input_type, output_path)