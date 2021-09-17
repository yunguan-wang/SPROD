import os
from posixpath import abspath
import sys
import argparse
import logging
from multiprocessing import Pool
from feature_extraction import extract_img_features
from pseudo_image_gen import make_pseudo_img
from slide_seq_stiching import stiching_subsampled_patches
from slideseq_make_patches import subsample_patches
import subprocess

'''
Dependency
----------
R >= 4.0.2
Python 3.7

positional arguments:
  input_path            Input folder containing all necessary files.
  output_path           Output path

optional arguments:
  -h, --help            show this help message and exit
  --input_type INPUT_TYPE, -y INPUT_TYPE
                        Input image type, select from {'single','patches'}.
                        (default: single)
  --output_prefix OUTPUT_PREFIX, -p OUTPUT_PREFIX
                        Project name, First part of names in the output.
                        (default: sprod)
  --sprod_npc SPROD_NPC, -sn SPROD_NPC
                        Number of PCs to use. positive integers. set to -1 to
                        use all PCs, and use the original IF matrix (default:
                        -1)
  --sprod_umap, -su     Toggle to use UMAP on top of PCA to represent
                        features. (default: False)
  --sprod_R SPROD_R, -r SPROD_R
                        Spot neighborhood radius ratio, 0-1,
                        radius=R*min(xmax-xmin,ymax-ymin). (default: 0.08)
  --spord_perplexity SPORD_PERPLEXITY, -u SPORD_PERPLEXITY
                        Perplexity, used in Sprod to find the proper Gaussian
                        kernal for distant representation. (default: 250)
  --sprod_margin SPROD_MARGIN, -g SPROD_MARGIN
                        Margin for bisection search, used in Sprod to find the
                        proper Gaussian kernal. smaller = slower => accurate u
                        (default: 0.001)
  --sprod_latent_dim SPROD_LATENT_DIM, -k SPROD_LATENT_DIM
                        Dimension of the latent space used in sprod to
                        represent spots. (default: 10)
  --sprod_graph_reg SPROD_GRAPH_REG, -l SPROD_GRAPH_REG
                        Regularizer for spot graph contructed in sprod.
                        (default: 1)
  --sprod_weight_reg SPROD_WEIGHT_REG, -w SPROD_WEIGHT_REG
                        regularizer for the weights used to normalize
                        expression matrix. (default: 0.625)
  --sprod_diag, -d      Toggle to force graph weights to be diagnoal, useful
                        in reducing over smoothing (default: False)
  --image_feature_type IMAGE_FEATURE_TYPE, -i IMAGE_FEATURE_TYPE
                        Type of feature extracted. combination from {'spot',
                        'block'} and {'intensity', 'texture'} (default:
                        spot_intensity)
  --warm_start, -ws     Toggle for warm start, meaning the folder will have
                        all necessary files for sprod. (default: False)
  -ci IMG_TYPE, --img_type IMG_TYPE
                        Cold start option. File name for patches spots
                        location file. Only works when --type is single
                        (default: he)
'''

# Rscript /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/script/denoise_bing.R \
# -e /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/data/Counts.txt \
# -c /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/data/Spot_metadata.csv \
# -f /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/data/Spot_level_haralick_features.csv \
# -n 3 -x -r 0.08 -u 250 -k 10 -l 0.5 -t 0.625 -m 0.01 -d \
# -s /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise \
# -o /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/data/ \
# -p project_ID

def sprod_worker(cmd):
    '''
    (SPROD_PATH, CTS_FN, METADATA_FN, FEATURES_FN, N_PC, R, U, K, Lambda, L_E, M, UTIL_PATH, output)

    make_option(c("-e", "--Exp"), help="Expression matrix"),
    make_option(c("-c", "--Cspot"), help="Spot metadata, contains X, Y coordinates."),
    make_option(c("-f", "--ImageFeature"), help="Extracted image features"),
    make_option(c("-n", "--numberPC"), help="# Number of PCs to use"),
    make_option(c("-r", "--Rratio"),help="Spot neighborhood radius ratio, 0-1, radius=R*min(xmax-xmin,ymax-ymin)"), 
    make_option(c("-u", "--U"),help="# perplexity, Tunable"),
    make_option(c("-k", "--K"),help="latent space dimension, Tunable"), 
    make_option(c("-l", "--lambda"),help="# regularizer, tunable"),
    make_option(c("-t", "--L_E"),help="regularizer for denoising"), 
    make_option(c("-m", "--margin"), help="Margin for bisection search, smaller = slower => accuracy u"),
    make_option(c("-d","--diagNeighbor"),help = "To prevent over smooth."),
    make_option(c("-s", "--scriptPath"), help="path to the software folder"),
    make_option(c("-o", "--outputPath"), help="Output path"),
    make_option(c("-p", "--projectID"), help="# project name, First part of names in the output"),
    '''
    log_fn = cmd[-1]
    cmd = cmd[:-1]
    print(cmd)
    with open(log_fn, 'a') as outfile:
        subprocess.run(cmd, stdout=outfile, stderr = outfile)
    # os.system('rm -r output_fn') # Get rid of output.
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Main function for running sprod.",
    )

    parser.add_argument(
        "input_path", type=str, help="Input folder containing all necessary files."
    )

    parser.add_argument(
        "output_path", type=str, help="Output path"
    )

    parser.add_argument(
        "--input_type",
        "-y",
        default='single',
        type=str,
        help="Input image type, select from {'single','patches'}.",
    )

    parser.add_argument(
        "--output_prefix",
        "-p",
        default="sprod",
        type=str,
        help="Project name, First part of names in the output.",
    )
    parser.add_argument(
        "--sprod_npc",
        "-sn",
        type=int,
        default=-1,
        help="Number of PCs to use. positive integers. \
            set to -1 to use all PCs, and use the original IF matrix",
    )
    parser.add_argument(
        "--sprod_umap",
        "-su",
        default=False,
        action="store_true",
        help="Toggle to use UMAP on top of PCA to represent features.",
    )
    parser.add_argument(
        "--sprod_R",
        "-r",
        default=0.08,
        type=float,
        help="Spot neighborhood radius ratio, 0-1, \
            radius=R*min(xmax-xmin,ymax-ymin).",
    )
    parser.add_argument(
        "--spord_perplexity",
        "-u",
        default=250,
        type=int,
        help="Perplexity, used in Sprod to find the proper Gaussian kernal \
            for distant representation.",
    )
    parser.add_argument(
        "--sprod_margin",
        "-g",
        default=0.001,
        type=float,
        help="Margin for bisection search, used in Sprod to find the proper \
            Gaussian kernal. smaller = slower => accurate u",
    )
    parser.add_argument(
        "--sprod_latent_dim",
        "-k",
        default=10,
        type=int,
        help="Dimension of the latent space used in sprod to represent spots.",
    )
    parser.add_argument(
        "--sprod_graph_reg",
        "-l",
        default=1,
        type=float,
        help="Regularizer for spot graph contructed in sprod.",
    ),
    parser.add_argument(
        "--sprod_weight_reg",
        "-w",
        default=0.625,
        type=float,
        help="regularizer for the weights used to normalize expression matrix.",
    )
    parser.add_argument(
        "--sprod_diag",
        "-d",
        action="store_true",
        default=False,
        help="Toggle to force graph weights to be diagnoal, useful in reducing \
            over smoothing",
    )

    parser.add_argument(
        "--image_feature_type",
        "-i",
        type=str,
        default="spot_intensity",
        help="Type of feature extracted. combination from {'spot', 'block'} and \
            {'intensity', 'texture'}",
    )

    parser.add_argument(
        "--warm_start",
        "-ws",
        default=False,
        action="store_true",
        help="Toggle for warm start, meaning the folder will have all \
            necessary files for sprod.",
    )

    parser.add_argument(
        "--num_of_patches",
        "-pn",
        type=int,
        default=10,
        help="Number of subsampled patches. Only works when --type is patches.",
    )

    parser.add_argument(
        "--num_of_batches",
        "-pb",
        type=int,
        default=10,
        help="How many times subsampling is ran. Only works when --type is patches.",
    )

    # parser.add_argument(
    #     "-css",
    #     "--patches_spots",
    #     type=str,
    #     help="Cold start option. File name for patches spots location file. Only works when --type is patches.",
    # )

    parser.add_argument(
        "-ci",
        "--img_type",
        type=str,
        default='he',
        help="Cold start option. File name for patches spots location file. " +
         "Only works when --type is single",
    )

    args = parser.parse_args()
    # args = parser.parse_args(
    #     ['/home2/s190548/work_xiao/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/subsample_patches',
    #     '/home2/s190548/work_xiao/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/subsample_patches/denoised'])
    # Todo: decide if checking metadata function should be added, or force the user to provide the correct format.
    input_path = args.input_path
    output_path = args.output_path
    input_type = args.input_type
    warm_start = args.warm_start
    output_prefix = args.output_prefix
    sprod_npc = args.sprod_npc
    sprod_umap = args.sprod_umap
    sprod_R = args.sprod_R
    spord_perplexity = args.spord_perplexity
    sprod_margin = args.sprod_margin
    sprod_latent_dim = args.sprod_latent_dim
    sprod_graph_reg = args.sprod_graph_reg
    sprod_weight_reg = args.sprod_weight_reg
    sprod_diag = args.sprod_diag
    image_feature_type = args.image_feature_type
    sprod_path = os.path.abspath(__file__)
    sprod_path = '/'.join(sprod_path.split('/')[:-1])
    pn = args.num_of_patches
    pb = args.num_of_batches

    if not os.path.exists(output_path):
        os.makedirs(output_path)

    log_fn = os.path.join(output_path, 'sprod_log.txt')
    if os.path.exists(log_fn):
        os.remove(log_fn)
    logging.basicConfig(
        filename=log_fn,
        format='%(asctime)s,%(levelname)s:::%(message)s',
        datefmt='%H:%M:%S',
        level='DEBUG')

    class StreamToLogger(object):
        """
        Fake file-like stream object that redirects writes to a logger instance.
        """
        def __init__(self, logger, log_level=logging.INFO):
            self.logger = logger
            self.log_level = log_level
            self.linebuf = ''

        def write(self, buf):
            for line in buf.rstrip().splitlines():
                self.logger.log(self.log_level, line.rstrip())
        def flush(self):
            pass

    stdout_logger = logging.getLogger('STDOUT')
    sl = StreamToLogger(stdout_logger, logging.INFO)
    sys.stdout = sl

    stderr_logger = logging.getLogger('STDERR')
    sl = StreamToLogger(stderr_logger, logging.ERROR)
    sys.stderr = sl

    sprod_script = os.path.join(sprod_path, 'denoise.R')
    if not warm_start:
        logging.info('Cold start, processing from counts and image data.')
        # check if input path contains necessary data
        if not (
            os.path.exists(os.path.join(input_path, 'Counts.txt')) == 
            os.path.exists(os.path.join(input_path, 'Spot_metadata.csv')) == True
            ):
            logging.error(
                'Counts.txt and/or Spotmeta.csv is missing from {}'.format(input_path)
                )
            raise ValueError(
                'Counts.txt and/or Spotmeta.csv is missing from {}'.format(input_path))
        # check if the input path have an tif image, if not, make a pseudo image.
        num_tifs = len([x for x in os.listdir(input_path) if x[-4:] == '.tif'])
        if num_tifs == 1:
            logging.info('Extracting intensity and texture features from matching image.')
            img_type = args.img_type
            _ = extract_img_features(
                input_path, img_type, input_path)
        elif num_tifs == 0:
            logging.info('Use spot cluster probability as pseudo image features')
            cts_fn = os.path.join(input_path, 'Counts.txt')
            spots_fn = os.path.join(input_path, 'Spot_metadata.csv')
            dp_script_path = os.path.join(
                sprod_path, 'dirichlet_process_clustering.R')
            make_pseudo_img(cts_fn, spots_fn, input_path, 'dp', dp_script_path)
        else:
            logging.error(
                'More than one images are present. Please remove the unwanted ones.'
                )
            raise ValueError('More than one images are present. Please remove the unwanted ones.')
        
        if input_type == 'patches':
            logging.info('Making subsample patches from counts data.')
            intermediate_path = os.path.join(output_path, 'intermediate')
            feature_fn = os.path.join(
                input_path, 
                "{}_level_{}_features.csv".format(*image_feature_type.split('_'))
                )
            if not os.path.exists(intermediate_path):
                os.makedirs(intermediate_path)
            subsample_patches(input_path, intermediate_path, feature_fn, pn, pb)
            input_path = intermediate_path
        elif input_type == 'single':
            pass
        else:
            logging.error('Input type must be single or patches')
            raise ValueError()

    logging.info('Proceed to Sprod denoising.')
    if input_type == 'patches':
        logging.info('Sprod-ready data format is patches')
        counts = [x for x in os.listdir(input_path) if 'Counts' in x]
        patches = [x.replace('_Counts.txt','') for x in counts]
        inputs = []
        for j, patch in enumerate(sorted(patches)):
            cts_fn = os.path.join(input_path, patch + '_Counts.txt')
            meta_fn = os.path.join(input_path, patch + '_Spot_metadata.csv')
            feature_fn = os.path.join(input_path, patch + '_F.csv')
            if not (
                os.path.exists(cts_fn) == 
                os.path.exists(meta_fn) == 
                os.path.exists(feature_fn) == True
                ):
                continue
            inputs.append([cts_fn, meta_fn, feature_fn])
    else:
        logging.info('Sprod-ready data format is single')
        cts_fn = os.path.join(input_path, 'Counts.txt')
        meta_fn = os.path.join(input_path, 'Spot_metadata.csv')
        feature_fn = os.path.join(
            input_path, "{}_level_{}_features.csv".format(*image_feature_type.split('_')))
        if not (
                os.path.exists(cts_fn) == 
                os.path.exists(meta_fn) == 
                os.path.exists(feature_fn) == True
                ):
                logging.error('Counts.txt and/or Spotmeta.csv is missing from {}'.format(input_path)
                )
                raise ValueError('Required data is missing!')
        inputs = [[cts_fn, meta_fn, feature_fn]]

    logging.info('Total number of sprod jobs : {}'.format(len(inputs)))
    # SPROD_PATH, CTS_FN, METADATA_FN, FEATURES_FN, N_PC, R, U, K, Lambda, L_E, M, UTIL_PATH, output)
    list_sprod_cmds = []
    for sprod_job in inputs:
        cts_fn, meta_fn, feature_fn = [os.path.abspath(x) for x in sprod_job]
        if input_type == 'patches':
            patch_name = cts_fn.split('/')[-1].replace('_Counts.txt','')
        else:
            patch_name = 'sprod'
        # sprod_cmd = 'Rscript ' + sprod_script + ' '
        sprod_cmd = ['Rscript',sprod_script]
        for sprod_key, param in zip(
            [
                '-e', '-c', '-f', 
                '-n','-r','-u',
                '-k', '-l', '-t',
                '-m','-x','-d',
                '-s', '-o', '-p'],
            [
                cts_fn, meta_fn, feature_fn,
                sprod_npc, sprod_R, spord_perplexity,
                sprod_latent_dim, sprod_graph_reg, sprod_weight_reg,
                sprod_margin, sprod_umap, sprod_diag,
                sprod_path.replace('/script',''), os.path.abspath(output_path),
                patch_name]
                ):
            if isinstance(param,bool):
                if param:
                    # sprod_cmd = sprod_cmd + str(sprod_key) + ' '
                    sprod_cmd.append(str(sprod_key))
            else:
                # sprod_cmd = sprod_cmd + sprod_key + ' ' + str(param) + ' '
                sprod_cmd.append(sprod_key)
                sprod_cmd.append(str(param))
        # log file name is packed into cmd
        list_sprod_cmds.append(sprod_cmd + [log_fn])
    if len(list_sprod_cmds) == 1:
        sprod_worker(list_sprod_cmds[0])
    elif len(list_sprod_cmds) >1:
        with Pool(16) as p:
            _ = p.map(sprod_worker, list_sprod_cmds)
    if input_type == 'patches':
        stiching_script_path = os.path.join(sprod_path,'slide_seq_stiching.py')
        stiching_subsampled_patches(output_path, os.path.join(output_path,'denoised_stiched.hdf'))
