import pandas as pd
import os
import sys
import numpy as np
import subprocess


def stiching_denoised_patches(input_path, output_fn = None):
    cts_files = sorted([x for x in os.listdir(input_path) if 'Counts.txt' in x])
    denoised_fns = sorted([x for x in os.listdir(input_path) if 'denoised' in x])
    assert(len(cts_files) == len(denoised_fns)), 'Slideseq subsampled data not properly denoised'

    denoised_mtx = pd.DataFrame()
    for cts_fn in cts_files:
        abs_cts_fn = os.path.join(input_path,cts_fn)
        core_name = cts_fn.replace('_Counts.txt','')
        metadata_fn = abs_cts_fn.replace('Counts.txt','Spot_metadata.csv')
        denoised_fn = abs_cts_fn.replace('Counts.txt', 'denoised/Denoised_matrix.txt')
        metadata = pd.read_csv(metadata_fn, index_col=0)
        denoised_cts = pd.read_csv(denoised_fn, index_col=0, sep='\t')
        core_idx = metadata[metadata.patch_core == core_name].index
        denoised_mtx = denoised_mtx.append(denoised_cts.loc[core_idx])
    
    if output_fn is None:
        output_fn = input_path.replace('/patches','/denoised_cts.h5df')

    denoised_mtx.to_hdf(output_fn, key='denoised')

def stiching_subsampled_patches(input_path, output_fn):
    
    def get_files_len(filenames):
        total_n = 0
        for fn in filenames:
            with open(fn) as f:
                for i, _ in enumerate(f):
                    pass
                total_n += i # ignore first line
        return total_n

    cts_files = sorted([x for x in os.listdir(input_path) if 'Denoised' in x])
    batches = list(set(['_'.join(x.split('_')[:2]) for x in cts_files]))
    n_patches = [
        len([y for y in cts_files if x.split('_')[1]  == y.split('_')[1]]) for x in batches
        ]
    n_patches = max(n_patches)
    
    # Make sure all batches are good.
    batch_dict = {}
    for _, batch in enumerate(batches):
        patch_cts = [x for x in cts_files if batch in x]
        if len(patch_cts) != n_patches:
            print(batch, 'did not finished properly and will be skipped.')
            continue
        n_barcodes = get_files_len([os.path.join(input_path, x) for x in patch_cts])
        batch_dict[batch] = n_barcodes

    barcodes_counts = np.unique([x for x in batch_dict.values()], return_counts=True)
    n_total_cells = barcodes_counts[0][np.argmax(barcodes_counts[1])]
    good_batches = {
        key: val for key, val in batch_dict.items() if val == n_total_cells}
    good_batches = good_batches.keys()
    
    # Pooling data.
    tmp = pd.read_csv(
        os.path.join(input_path, cts_files[0]),
        index_col=0, sep='\t',nrows=5)
    n_genes = tmp.shape[1]
    gene_names = tmp.columns
    pooled_barcodes = []
    for batch in good_batches:
        patch_cts = [x for x in cts_files if batch in x]
        patch_cts = [os.path.join(input_path, x) for x in patch_cts] 
        cts_array = np.zeros((n_total_cells, n_genes))
        top = 0
        barcode_list = []
        print('Processing {}'.format(batch))
        for _, cts_fn in enumerate(patch_cts):
            cts = pd.read_csv(cts_fn, index_col=0, sep='\t')
            barcode_list += cts.index.tolist()
            delta = cts.shape[0]
            cts_array[top:top+delta,:] = cts.values
            top = top + delta

        barcode_list = np.array(barcode_list)
        barcode_order = np.argsort(barcode_list)
        cts_array = cts_array[barcode_order,:]
        barcode_list = barcode_list[barcode_order]

        if pooled_barcodes == []:
            pooled_barcodes = barcode_list
            pooled_cts = cts_array
        else:
            assert((barcode_list == pooled_barcodes).all()), 'Barcodes in batches does not match!'
            pooled_cts += cts_array

    pooled_cts = pooled_cts / len(good_batches)
    pooled_cts = pd.DataFrame(pooled_cts, index = pooled_barcodes, columns=gene_names)
    pooled_cts.round(2).to_csv(output_fn, sep = '\t')

if __name__ == '__main__':
    input_path = sys.argv[1]
    # input_path = '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/subsample_patches/denoised'
    try:
        output_fn = sys.argv[2]
        # output_fn = '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/subsample_patches/denoised_counts.hdf'
    except IndexError:
        output_fn = None
    # stiching_denoised_patches(input_path, output_fn)
    stiching_subsampled_patches(input_path, output_fn)




