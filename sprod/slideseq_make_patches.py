import pandas as pd
import os
from multiprocessing import Pool
import sys
import numpy as np

def mpl_writer(patch):
    patch_name = patch[0]
    output_path = patch[1]
    patch_meta = patch[2]
    patch_cts = patch[3]
    patch_f = patch[4]
    print('Processing subsample {}'.format(patch_name))
    patch_cts.to_csv(
        os.path.join(output_path, patch_name + '_Counts.txt'),sep='\t')
    # check if there is i/o error
    _cts = pd.read_csv(
        os.path.join(output_path, patch_name + '_Counts.txt'),
        sep='\t', index_col=0)
    n = 0
    while not (patch_cts.index == _cts.index).all():
        if n > 5:
            print('Persistant I/O error in subsample {}'.format(patch_name))
            break
        patch_cts.to_csv(
            os.path.join(output_path, patch_name + '_Counts.txt'),sep='\t')
        _cts = pd.read_csv(
            os.path.join(output_path, patch_name + '_Counts.txt'),
            sep='\t', index_col=0)
        n += 1
    else:
        patch_meta.to_csv(
            os.path.join(output_path, patch_name + '_Spot_metadata.csv'))
        patch_f.to_csv(
            os.path.join(output_path, patch_name + '_F.csv'))

def make_patches(input_path, margin, patch_scale,output_path):
    spot_meta = pd.read_csv(
        os.path.join(input_path,'Spot_metadata.csv'), index_col=0
    )
    cts = pd.read_csv(
        os.path.join(input_path,'Counts.txt'), index_col=0, sep='\t'
    )
    # Use cpm (1e4) now
    cts = 1e4 * cts.apply(lambda x: x/x.sum(), axis=1)
    features = pd.read_csv(
        os.path.join(input_path, 'pseudo_image_features.csv'), index_col=0
    )

    # margin should be 1/scale * denoise neighborhood radius. 
    x_min, x_max = spot_meta.X.describe()[['min','max']]
    y_min, y_max = spot_meta.Y.describe()[['min','max']]

    delta_x = 1/patch_scale * (x_max-x_min)
    delta_y = 1/patch_scale * (y_max-y_min)
    margin_x = margin * (x_max-x_min)
    margin_y = margin * (y_max-y_min)

    # Determine patch core samples
    x_bins = pd.cut(
        spot_meta.X,patch_scale, labels=['x_bin_' + str(x+1) for x in range(patch_scale)])
    y_bins = pd.cut(
        spot_meta.Y,patch_scale, labels=['y_bin_' + str(x+1) for x in range(patch_scale)])
    spot_meta['patch_core'] = x_bins.astype(str) + '_' + y_bins.astype(str)

    # Output patches
    patches = []
    for i in range(patch_scale):
        x_left = x_min + i*delta_x - margin_x
        x_right = x_min + (i+1)*delta_x + margin_x
        for j in range(patch_scale):
            y_left = y_min + j*delta_y - margin_y
            y_right = y_min + (j+1)*delta_y + margin_y
            _tmp = spot_meta[
                (x_left < spot_meta.X) & (spot_meta.X <= x_right) & 
                (y_left < spot_meta.Y) & (spot_meta.Y <= y_right)]
            if _tmp.shape[0] > 10:
                patch_meta = _tmp
                patch_name = patch_meta.patch_core.value_counts().index[0]
                patches.append(
                    [
                    patch_name, output_path, patch_meta, 
                    cts.loc[patch_meta.index], features.loc[patch_meta.index]
                    ])
    # use mp to speed up IO.
    with Pool(16) as p:
        _ = p.map(mpl_writer, patches)

def subsample_patches(input_path, output_path, feature_fn = None, n_patches = 10, number_batches = 10):
    """
    sample counts 5000-spot subsamples, which is called as a patch. By default this is done number_batches times, \
    making n_patches*number_batches patches where each cell is represented number_batches times.
    """
    spot_meta = pd.read_csv(
        os.path.join(input_path,'Spot_metadata.csv'), index_col=0
    )
    cts = pd.read_csv(
        os.path.join(input_path,'Counts.txt'), index_col=0, sep='\t'
    )
    # Use cpm (1e4) now
    # cts = 1e4 * cts.apply(lambda x: x/x.sum(), axis=1)
    if feature_fn is None:
        features = pd.read_csv(
            os.path.join(input_path, 'pseudo_image_features.csv'), index_col=0
        )
    elif feature_fn[-3:] == 'csv':
        features = pd.read_csv(feature_fn, index_col=0)
    elif feature_fn[-3:] == 'txt':
        features = pd.read_csv(feature_fn, index_col=0, sep='\t')
    patches = []
    for i in range(number_batches):
        patch_ids = np.random.randint(0,n_patches,size = cts.shape[0])
        for patch_id in range(n_patches):
            mask = (patch_ids == patch_id)
            patch_name = 'Batch_{}_patch_{}'.format(i, patch_id)
            patch_meta = spot_meta[mask]
            patch_cts = cts[mask]
            patch_features = features[mask]
            patches.append([
                patch_name,output_path, patch_meta,patch_cts,patch_features
            ])
    # use mp to speed up IO.
    with Pool(16) as p:
        _ = p.map(mpl_writer, patches)

if __name__ == '__main__':
    # cmd running mode will process all slideseq data.
    # input_path = '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08'
    try:
        input_path = sys.argv[1]
    # if no input_path are defined, using default cmd running mode.
    except IndexError:
        input_path = '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/'
    try:
        margin = float(sys.argv[2])
    except IndexError:
        margin = 0.01 # This is equivalent to R = 0.1

    try:
        patch_scale = sys.argv[3]
    except IndexError:
        patch_scale = 4
    
    for slideseq_path in os.listdir(input_path):
        slideseq_path = os.path.abspath(os.path.join(input_path,slideseq_path))
        # output_path = os.path.join(slideseq_path,'patches')
        output_path = os.path.join(slideseq_path,'subsample_patches')
        if os.path.exists(output_path):
            os.system('rm -r {}'.format(output_path))
        os.makedirs(output_path)
        print('Processing {}'.format(slideseq_path))
        # make_patches(slideseq_path, margin, patch_scale,output_path)
        subsample_patches(slideseq_path, output_path)
