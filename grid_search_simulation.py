'''
A wrapper job submitting script for hyperparameter optimization based on grid search.
Need to load R/4.0 after module purge and module load shared.
Need to load the mocca env
'''
import os
from multiprocessing import Pool
from itertools import product
import argparse
import pandas as pd

def sprod_worker(params):
    '''
    (SPROD_PATH, CTS_FN, METADATA_FN, FEATURES_FN) + [N_PC,R,U,K,Lambda,L_E,M] + (UTIL_PATH, output)

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
    sprod_params = params[:-1]
    testing = params[-1]
    if testing:
        print('Rscript {} -e {} -c {} -f {} -n {} -r {} -u {} -k {} -l {} -t {} -m {} -s {} -o {}'.format(*sprod_params))
    else:
        os.system('Rscript {} -e {} -c {} -f {} -n {} -r {} -u {} -k {} -l {} -t {} -m {} -s {} -o {}'.format(*sprod_params))
    # os.system('rm -r output_fn') # Get rid of output.
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Grid search script for denoise.R. Parameters are passed as initial value, step multiplier and number of steps.\
            Seperated by "-".'
    )

    parser.add_argument(
        "--N_PC", '-P',
        type = str,
        default = '3-1-1',
        help="Number of principle components."
    )

    parser.add_argument(
        "--N_latent_space",
        "-K",
        type = str,
        default = '3-1-1',
        help = "Latent space dimension.",
    )

    parser.add_argument(
        "--neighborhood_radius_ratio",
        "-R",
        default = '0.05-2-2',
        type=str,
        help="Size of neighborhood radius to search for proximity relation.",
    )

    parser.add_argument(
        "--perplexity",
        "-U",
        type=str,
        default='1000-0.5-3',
        help="Perplexity for gaussian kernal during distance calculation.",
    )

    parser.add_argument(
        "--max_edge_weight",
        "-S",
        type=str,
        default='1-1-1',
        help="Upper bound of graph edge weights.",
    )

    parser.add_argument(
        "--graph_regularizer",
        "-L",
        type=str,
        default='10-0.5-4',
        help="Regularization constant for graph.",
    )

    parser.add_argument(
        "--denoising_regularizer",
        "-L_E",
        type=str,
        default='1-0.5-2',
        help="Regularization constant for graph.",
    )

    parser.add_argument(
        "--bisection_margin",
        "-M",
        type=str,
        default='0.001-1-1',
        help="Margin for bisection search.",
    )

    parser.add_argument(
        "--result_path",
        "-r",
        type=str,
        default='Grid_search_results.csv',
        help="Output filename.",
    )

    parser.add_argument(
        "--cts_fn",
        "-c",
        type=str,
        default='/project/shared/xiao_wang/projects/MOCCA/data/simulation/Es_0.1n_5k.txt',
        help="Input expression data with noise.",
    )

    parser.add_argument(
        "--features_fn",
        "-f",
        type=str,
        default='/project/shared/xiao_wang/projects/MOCCA/data/simulation/IF_0.1n_5k.csv',
        help="Input (pseoudo) image features.",
    )

    parser.add_argument(
        "--cts_true_fn",
        "-ct",
        type=str,
        default='/project/shared/xiao_wang/projects/MOCCA/data/simulation/Es_5k.txt',
        help="Ground truth gene expression data.",
    )

    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default='/project/shared/xiao_wang/projects/MOCCA/data/grid_search',
        help="Ground truth gene expression data.",
    )

    parser.add_argument(
        "--step_method",
        "-m",
        type=str,
        default='multiply',
        help="Parameter search method. {'add','multiply'}",
    )

    parser.add_argument(
        "--testing_mode",
        "-t",
        default = False,
        action = 'store_true',
        help="Toggle testing mode",
    )
    args = parser.parse_args()
    # args = parser.parse_args([])

    N_PC = args.N_PC
    R = args.neighborhood_radius_ratio
    U = args.perplexity
    K = args.N_latent_space
    S = args.max_edge_weight
    Lambda = args.graph_regularizer
    L_E = args.denoising_regularizer
    M = args.bisection_margin
    result_fn = args.result_path
    CTS_FN = args.cts_fn
    FEATURES_FN = args.features_fn
    E_ORIGINAL = args.cts_true_fn
    OUTPUT_PATH = args.output
    STEP_METHOD = args.step_method
    testing = args.testing_mode

    # Hard coded script path and input file locations.
    UTIL_PATH = '/project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise'
    SPROD_PATH = '/project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/script/denoise.R'
    # CTS_FN = '/project/shared/xiao_wang/projects/MOCCA/data/simulation/Es_0.1n_5k.txt'
    METADATA_FN = '/project/shared/xiao_wang/projects/MOCCA/data/simulation/C_simu_5k.csv'
    # FEATURES_FN = '/project/shared/xiao_wang/projects/MOCCA/data/simulation/IF_0.1n_5k.csv'
    # OUTPUT_PATH = '/project/shared/xiao_wang/projects/MOCCA/data/grid_search'
    # E_ORIGINAL = '/project/shared/xiao_wang/projects/MOCCA/data/simulation/Es_5k.txt'
    # Exact all paramerters from argparse
    PARAMS = []
    for p in [N_PC,R,U,K,Lambda,L_E,M]:
        if ',' in p:
            PARAMS.append([float(x) for x in p.split(',')])
        else:
            val, step, n_steps = [float(x) for x in p.split('-')]
            if STEP_METHOD == 'multiply':
                PARAMS.append([val*step**i for i in range(int(n_steps))])
            elif STEP_METHOD == 'add':
                PARAMS.append([round(val + step*(i),3) for i in range(int(n_steps))])
    # Create job parameters for Pool
    job_params = []
    output_folders = []
    for p in product(*PARAMS):
        output = OUTPUT_PATH + '/' + '_'.join([str(x) for x in p])
        output_folders.append(output)
        if not os.path.exists(output):
            os.makedirs(output)
        else:
            if os.path.exists(os.path.join(output, 'Denoised_matrix.txt')):
                print('Job is previously done in {}.'.format(output))
                continue
        p = (SPROD_PATH, CTS_FN, METADATA_FN, FEATURES_FN) + p + (UTIL_PATH, output, testing)
        job_params.append(p)
    # multiprocessing job
    if len(job_params) > 0:
        with Pool(16) as p:
            _ = p.map(sprod_worker, job_params)
    if not testing:        
        # Save grid search results.
        res = pd.Series(0,index = [x.split('/')[-1] for x in output_folders], dtype=float)
        E = pd.read_csv(E_ORIGINAL, sep='\t')
        for output_fn in output_folders:
            denoised_fn = os.path.join(output_fn, 'Denoised_matrix.txt')
            k = output_fn.split('/')[-1]
            if os.path.exists(denoised_fn):
                E_denoised = pd.read_csv(os.path.join(output_fn, 'Denoised_matrix.txt'), sep='\t')
            else:
                print('Denoised job {} failed, skip this job'.format(k))
                continue
            sae = abs(E.values - E_denoised.values).sum()
            res[k] = sae
        res_params_array = [x.split('_') for x in res.index]
        res = pd.concat(
            [pd.DataFrame(res_params_array, index=res.index),pd.DataFrame(res)], axis=1)
        res.columns = ['N_PC', 'R', 'U', 'K', 'Lambda', 'L_E', 'Margin', 'SAE_vs_True_E']
        res.to_csv(result_fn)

# Rscript /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/script/denoise.R \
# /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/data/Counts.txt \
# /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/data/Spot_metadata.csv \
# /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/data/Spot_level_haralick_features.csv \
# 5 0.2 1000 5 5 10 1 0.001 \
# /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise \
# /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/data/