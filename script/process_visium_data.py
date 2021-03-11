'''
Process Visium expression data into cell by gene matrix txt file
Need p3s env, or the imported packages.
'''
import pandas as pd
from skimage import io, morphology,color, img_as_float32
import matplotlib.pyplot as plt
from scipy.io import mmread
import os
import json
import numpy as np

fn_dir = '/project/shared/xiao_wang/projects/MOCCA/data/Visium/ovarian_cancer_immune/'
os.chdir(fn_dir)
img_type = {'HE':'HE','IF':'IF'}['IF']

mtx_fn = [x for x in os.listdir() if ('filtered' in x) and ('.gz') in x][0]
os.system("tar -xzvf {}".format(mtx_fn))
mtx = mmread(
    os.path.join('filtered_feature_bc_matrix', 'matrix.mtx.gz')
)
genes = pd.read_csv(
    os.path.join('filtered_feature_bc_matrix', 'features.tsv.gz'), 
    sep='\t', header=None,
    ).iloc[:,1].values
barcodes = pd.read_csv(
    os.path.join('filtered_feature_bc_matrix', 'barcodes.tsv.gz'),
    sep='\t', header=None).iloc[:,0].values
mtx = mtx.toarray()
cts = pd.DataFrame(mtx, columns=barcodes,index=genes)
cts = cts.groupby(cts.index).mean().T

'''
process visium image data, and extract channels
'''
img_fn = [x for x in os.listdir() if '.tif' in x][0]
spatial_fn = [x for x in os.listdir() if 'spatial.tar.gz' in x][0]
os.system("tar -xzvf {}".format(spatial_fn))
img_json = [x for x in os.listdir('spatial') if '.json' in x][0]
img_meta = pd.read_csv(
    'spatial/tissue_positions_list.csv', index_col=0, header=None)
with open('spatial/'+img_json) as f:
    img_json = json.load(f)
try:
    # Need to use spot diamater here
    radius = int(img_json['fiducial'][0]['dia'] // 2)
except:
    radius = int(img_json['spot_diameter_fullres']// 2)
    

img = io.imread(img_fn)
if img_type == 'HE':
    img_he = color.rgb2hed(img)
    c1 = 255*(img_he[:,:,0] - img_he[:,:,0].min())/np.ptp(img_he[:,:,0])
    c1 = c1.astype('uint8')
    c2 = 255*(img_he[:,:,1] - img_he[:,:,1].min())/np.ptp(img_he[:,:,1])
    c2 = c2.astype('uint8')
elif img_type == 'IF':
    img = img[:,:,[3,1,0]]
    img = img_as_float32(img)
    img = (255*img).astype('uint8')
    _ = plt.figure(figsize=(64,36))
    io.imshow(img)
    plt.savefig('Original_figure_CD45_CK_DAPI.pdf')
    plt.close()
    c1 = img[:,:,0]
    c2 = img[:,:,2]


masks = np.zeros(img.shape[:2], dtype='bool')
spot_mask = morphology.disk(radius, dtype='bool')
img_meta.columns = ['Unknown_column','Row','Col','X','Y']
img_meta['C1_median'] = 0
img_meta['C2_median'] = 0
img_meta['Spot_radius'] = radius
k = 1
for i in range(img_meta.shape[0]):
    x,y = img_meta.iloc[i,3:5].astype(int)
    masks[x,y] = 1
    x0 = x-radius
    x1 = x+radius + 1
    y0 = y-radius
    y1 = y+radius + 1
    c1_median = np.median(c1[x0:x1,y0:y1][spot_mask])
    c2_median = np.median(c2[x0:x1,y0:y1][spot_mask])
    img_meta.iloc[i,5:7] = [c1_median,c2_median]
    k+=1
    if k%1000 == 0:
        print('Processed {} samples'.format(k))

cm_cells = [x for x in img_meta.index if x in cts.index]
cts.loc[cm_cells].to_csv('Counts.txt',sep='\t')
img_meta.loc[cm_cells].to_csv('Spot_metadata.csv')