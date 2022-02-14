'''
Process Visium expression data into cell by gene matrix txt file
'''
import pandas as pd
from scipy.io import mmread
import os
import argparse
import json

# parsing arguments

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Process Visium or Slide-seq data into formats accepted by downstream analysis.',
)

parser.add_argument(
    "--input_type", '-t',
    type=str,
    help="Input dataset type. Must be one from ['visium', 'slideseq']"
)

parser.add_argument(
    "-v",
    "--visium_folder",
    type=str,
    help="Visium dataset folder path. Only works when --type is visium.",
)

parser.add_argument(
    "-sc",
    "--slideseq_counts",
    type=str,
    help="File name for slideseq counts file. Only works when --type is slideseq.",
)

parser.add_argument(
    "-ss",
    "--slideseq_spots",
    type=str,
    help="File name for slideseq spots location file. Only works when --type is slideseq.",
)

parser.add_argument(
    "-o",
    "--output",
    type=str,
    nargs="?",
    default='./intermediate',
    help="Output path for returned counts file and spotmetadata.",
)

def data_preprocessing(input_type, output_path, cts_fn = None, spots_fn = None):
    if input_type == 'visium':
        fn_dir = args.visium_folder
        # Process counts files.
        os.chdir(fn_dir)
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
        extract image metadata, and figure out spot radius.
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
        
        # Reformat image metadata.
        img_meta = img_meta.iloc[:,1:]
        img_meta.columns = ['Row','Col','X','Y']
        img_meta['Spot_radius'] = radius
        os.system('cp {} {}'.format(os.path.join(fn_dir, img_fn), output_path))

    elif input_type == 'slideseq':
        chunks = pd.read_csv(cts_fn, index_col=0, sep='\t', comment = '#', chunksize=1000)
        cts = pd.DataFrame()
        for chunk in chunks:
            cts = pd.concat([cts,chunk.transpose()], axis=1)

        img_meta = pd.read_csv(spots_fn, index_col=0)
        img_meta.columns = ['X','Y']

    else:
        raise ValueError('Input dataset type must be either visium or slideseq')

    # Align data and output
    cm_cells = [x for x in img_meta.index if x in cts.index]
    cts = cts.loc[cm_cells]
    cts.to_csv(os.path.join(output_path,'Counts.txt'),sep='\t')
    img_meta = img_meta.loc[cm_cells]
    img_meta.to_csv(os.path.join(output_path,'Spot_metadata.csv'))

if __name__ == '__main__':
    args = parser.parse_args()
    input_type = args.input_type
    output_path = args.output
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    if os.path.exists(os.path.join(output_path,'Counts.txt')):
        print("Dataset is already processed in {}".format(output_path))
    else:
        if input_type == 'visium':
            data_preprocessing(input_type, output_path)
        else:
            cts_fn = args.slideseq_counts
            spots_fn = args.slideseq_spots 
            data_preprocessing(input_type, output_path, cts_fn, spots_fn)
    # if input_type == 'visium':
    #     fn_dir = args.visium_folder
    #     # Process counts files.
    #     os.chdir(fn_dir)
    #     mtx_fn = [x for x in os.listdir() if ('filtered' in x) and ('.gz') in x][0]
    #     os.system("tar -xzvf {}".format(mtx_fn))
    #     mtx = mmread(
    #         os.path.join('filtered_feature_bc_matrix', 'matrix.mtx.gz')
    #     )
    #     genes = pd.read_csv(
    #         os.path.join('filtered_feature_bc_matrix', 'features.tsv.gz'), 
    #         sep='\t', header=None,
    #         ).iloc[:,1].values
    #     barcodes = pd.read_csv(
    #         os.path.join('filtered_feature_bc_matrix', 'barcodes.tsv.gz'),
    #         sep='\t', header=None).iloc[:,0].values
    #     mtx = mtx.toarray()
    #     cts = pd.DataFrame(mtx, columns=barcodes,index=genes)
    #     cts = cts.groupby(cts.index).mean().T

    #     '''
    #     extract image metadata, and figure out spot radius.
    #     '''
    #     img_fn = [x for x in os.listdir() if '.tif' in x][0]
    #     spatial_fn = [x for x in os.listdir() if 'spatial.tar.gz' in x][0]
    #     os.system("tar -xzvf {}".format(spatial_fn))
    #     img_json = [x for x in os.listdir('spatial') if '.json' in x][0]
    #     img_meta = pd.read_csv(
    #         'spatial/tissue_positions_list.csv', index_col=0, header=None)
    #     with open('spatial/'+img_json) as f:
    #         img_json = json.load(f)
    #     try:
    #         # Need to use spot diamater here
    #         radius = int(img_json['fiducial'][0]['dia'] // 2)
    #     except:
    #         radius = int(img_json['spot_diameter_fullres']// 2)
        
    #     # Reformat image metadata.
    #     img_meta = img_meta.iloc[:,1:]
    #     img_meta.columns = ['Row','Col','X','Y']
    #     img_meta['Spot_radius'] = radius
    #     os.system('cp {} {}'.format(os.path.join(fn_dir, img_fn), output_path))
    #     # Align data and output
    # elif input_type == 'slideseq':
    #     cts_fn = args.slideseq_counts
    #     spots_fn = args.slideseq_spots
    #     chunks = pd.read_csv(cts_fn, index_col=0, sep='\t', comment = '#', chunksize=1000)
    #     cts = pd.DataFrame()
    #     for chunk in chunks:
    #         cts = pd.concat([cts,chunk.transpose()], axis=1)

    #     img_meta = pd.read_csv(spots_fn, index_col=0)
    #     img_meta.columns = ['X','Y']

    # else:
    #     raise ValueError('Input dataset type must be either visium or slideseq')

    # cm_cells = [x for x in img_meta.index if x in cts.index]
    # cts = cts.loc[cm_cells]
    # cts.to_csv(os.path.join(output_path,'Counts.txt'),sep='\t')
    # img_meta = img_meta.loc[cm_cells]
    # img_meta.to_csv(os.path.join(output_path,'Spot_metadata.csv'))

