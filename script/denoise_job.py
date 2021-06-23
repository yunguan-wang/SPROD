import os
import sys
import argparse

# Rscript /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/script/denoise_bing.R \
# -e /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/data/Counts.txt \
# -c /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/data/Spot_metadata.csv \
# -f /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/data/Spot_level_haralick_features.csv \
# -n 3 -x -r 0.08 -u 250 -k 10 -l 0.5 -t 0.625 -m 0.01 -d \
# -s /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise \
# -o /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/data/ \
# -p project_ID

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Main function for running sprod.",
    )
    parser.add_argument(
        "input_path", type=str, help="Input folder containing all necessary files."
    )
    parser.add_argument("output_path", type=str, help="Output path")
    parser.add_argument(
        "--input_type",
        "-y",
        default='slideseq',
        type=str,
        help="Input image type, select from {'visium','slideseq'}.",
    )
    parser.add_argument(
        "--warm_start",
        "-wr",
        default=False,
        action="store_true",
        help="Toggle for warm start, meaning the folder will have all \
            necessary files for sprod.",
    )
    parser.add_argument(
        "--sprod_path",
        "-s",
        default="/project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise",
        type=str,
        help="Path to the software folder",
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
        default=0.4,
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
        "--test_mode",
        "-T",
        action="store_true",
        default=False,
        help="For testing purpose, only output script, without submitting the job.",
    )

    # args = parser.parse_args()
    args = parser.parse_args(
        ['/home2/s190548/work_xiao/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/subsample_patches',
        '/home2/s190548/work_xiao/projects/MOCCA/data/Sprod_ready_data/slideseq/Puck_200115_08/subsample_patches/denoised'])
    # Todo: decide if checking metadata function should be added, or force the user to provide the correct format.
    input_path = args.input_path
    output_path = args.output_path
    input_type = args.input_type
    warm_start = args.warm_start
    sprod_path = args.sprod_path
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
    test_mode = args.test_mode

    if not os.path.exists(output_path):
        os.makedirs(output_path)
    sprod_script = os.path.join(sprod_path, 'script/denoise.R')
    if warm_start:
        if input_type == 'slideseq':
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
                if j % 5 == 4:
                    line_sep = '\n'
                else:
                    line_sep = '&\n'
                inputs.append([cts_fn, meta_fn, feature_fn, line_sep])
        else:
            cts_fn = os.path.join(input_path, 'Counts.txt')
            meta_fn = os.path.join(input_path, 'Spot_metadata.csv')
            feature_fn = os.path.join(
                input_path, "{}_level_{}_features.csv".format(*image_feature_type.split('_')))
            if not (
                    os.path.exists(cts_fn) == 
                    os.path.exists(meta_fn) == 
                    os.path.exists(feature_fn) == True
                    ):
                    raise ValueError('Required data is missing!')
            line_sep = '\n'
            inputs = [[cts_fn, meta_fn, feature_fn, line_sep]]
        with open(output_path + "/denoise_batch_job.sh", "w") as f:
            f.write("#!/bin/bash\n")
            f.write("\n")
            f.write("#SBATCH --job-name=sprod\n")
            f.write("#SBATCH --partition=256GB\n")
            f.write("#SBATCH --nodes=1\n")
            f.write("#SBATCH --ntasks=32\n")
            f.write("#SBATCH --time=5-00:00:00\n")
            f.write("#SBATCH --output=sprod_out.log\n")
            f.write("#SBATCH --error=sprod_err.log\n")
            f.write("module purge\n")
            f.write("module load shared\n")
            f.write("module load R/4.0.2-gccmkl\n")
            for _patch in inputs:
                cts_fn, meta_fn, feature_fn, line_sep = _patch
                patch = cts_fn.split('/')[-1].replace('_Counts.txt','')
                f.write("Rscript {} ".format(sprod_script))
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
                        sprod_path, output_path, patch]
                        ):
                    if isinstance(param,bool):
                        if param:
                            f.write("{} ".format(sprod_key))
                    else:
                        f.write("{} {} ".format(sprod_key, param))
                f.write(line_sep)
    else:
        pass
    if not test_mode:
        os.system('sbatch {}'.format(output_path + "/denoise_batch_job.sh"))

