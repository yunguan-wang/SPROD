import pandas as pd
from skimage import io
import os
import seaborn as sns
import numpy as np

cts = pd.read_csv(
    '/project/shared/xiao_wang/projects/MOCCA/data/simulation/Es_0.5n_5k.txt',
    sep='\t',index_col=0)
f = pd.read_csv(
    '/project/shared/xiao_wang/projects/MOCCA/data/simulation/IF_NO_5k.csv',
    index_col=0)
meta = pd.read_csv(
    '/project/shared/xiao_wang/projects/MOCCA/data/simulation/C_simu_5k.csv',
    index_col=0)


f_img = pd.concat((f.iloc[:,[2,5]], f.iloc[:,[0,1,3,4]].sum(axis=1)), axis=1)
f_img.columns = ['c1','c2','c3']
# sns.scatterplot(x = meta.X, y = meta.Y, hue = f_img.iloc[:,0], s = 5)
# sns.scatterplot(x = meta.X, y = meta.Y, hue = f_img.iloc[:,1], s = 5)
# sns.scatterplot(x = meta.X, y = meta.Y, hue = f_img.iloc[:,2], s = 5)

img_array = np.zeros((50,50,3))
img_e = []
img_meta = []
for i in range(50):
    for j in range(50):
        dist = np.sqrt((meta.X - i)**2 + (meta.Y - j)**2)
        near_spots = dist[dist<1].index
        if len(near_spots) == 0:
            continue
        else:
            img_array[i,j,:] = f_img.loc[near_spots].mean(axis=0)
            img_e.append(cts.loc[near_spots].mean().values)
            img_meta.append([i,j])

img_meta = pd.DataFrame(img_meta)
img_meta.columns = ['X','Y']

img_e = pd.DataFrame(img_e)
img_e.columns = cts.columns

img_meta.to_csv('../example/input/Spot_metadata.csv')
img_e.to_csv('../example/input/Counts.txt',sep='\t')
io.imsave('../example/input/image.tif', img_array)