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
    os.system('Rscript {} {} {} {} {} {} {} {} {} {} {} {} {} {}'.format(*params[:-1]))
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
    for p in [N_PC,R,U,K,S,Lambda,L_E,M]:
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
        p = (SPROD_PATH, CTS_FN, METADATA_FN, FEATURES_FN) + p + (UTIL_PATH, output, E_ORIGINAL)
        job_params.append(p)
    # multiprocessing job
    if len(job_params) > 0:
        with Pool(16) as p:
            _ = p.map(sprod_worker, job_params)
            
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
    res.columns = ['N_PC', 'R', 'U', 'K', 'S', 'Lambda', 'L_E', 'Margin', 'SAE_vs_True_E']
    res.to_csv(result_fn)

# Rscript /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/script/denoise.R \
# /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/data/Counts.txt \
# /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/data/Spot_metadata.csv \
# /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/data/Spot_level_haralick_features.csv \
# 5 0.2 1000 5 5 10 1 0.001 \
# /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise \
# /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/data/