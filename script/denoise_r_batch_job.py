import os
import sys

# Rscript /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/script/denoise.R \
# /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/data/Counts.txt \
# /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/data/Spot_metadata.csv \
# /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/data/Spot_level_haralick_features.csv \
# 5 0.2 1000 5 5 10 1 0.001 \
# /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise \
# /project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/data/

try:
    sprod_folder = sys.argv[1]
except IndexError:
    sprod_folder = '/project/shared/xiao_wang/projects/MOCCA/data/Sprod_ready_data'

try:
    parames_string = sys.argv[2]
except IndexError:
    parames_string = '5,0.2,1000,5,5,10,1,0.001'

# get all Counts files
counts = os.popen('find {} -name "Counts.txt"'.format(sprod_folder)).read()
counts = counts.split('\n')
counts.pop()
input_folders = [x[:-10] for x in counts]

r_script_path = '/project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise'
with open('./denoise_batch_job.sh','w') as f:
    f.write('module purge\n')
    f.write('module load shared\n')
    f.write('module load R/4.0.2-gccmkl\n')
    for input_folder in input_folders:
        if 'visium' in input_folder:
            output_path = os.path.join(input_folder,'denoised')
            if not os.path.exists(output_path):
                os.mkdir(output_path)
            feature_fn = 'Spot_level_haralick_features.csv'
            f.write('Rscript {} {} {} {} {} {} {}\n'.format(
                '/project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/script/denoise.R',
                input_folder + 'Counts.txt',
                input_folder + 'Spot_metadata.csv',
                input_folder + feature_fn,
                parames_string.replace(',', ' '),
                r_script_path,
                output_path
                )
            )
        else:
            if not 'patches' in os.listdir(input_folder):
                print('Slideseq patches not found. Skip {}'.format(input_folder))
                continue
            # if patches exists
            input_folder = os.path.join(input_folder, 'patches')
            cts_fns = sorted(
                [x for x in os.listdir(input_folder) if 'txt' in x])
            metadata_fns = sorted(
                [x for x in os.listdir(input_folder) if 'meta' in x])
            feature_fns = sorted(
                [x for x in os.listdir(input_folder) if '_F' in x])
            n_job = 1
            for cts_fn, meta_fn, F_fn in zip(cts_fns, metadata_fns, feature_fns):
                output_path = os.path.join(
                    input_folder,cts_fn.replace('_Counts.txt','_') + 'denoised')
                if not os.path.exists(output_path):
                    os.mkdir(output_path)
                if n_job%4 == 0: # 10 jobs at a time
                    line_end = '\n'
                else:
                    line_end = ' &\n'
                n_job += 1
                f.write('Rscript {} {} {} {} {} {} {}{}'.format(
                    '/project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/script/denoise.R',
                    input_folder + '/' + cts_fn,
                    input_folder + '/' + meta_fn,
                    input_folder + '/' + F_fn,
                    parames_string.replace(',', ' '),
                    r_script_path,
                    output_path,
                    line_end
                )
                )
            f.write('python {} {}\n'.format('/project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/script/slide_seq_stiching.py', input_folder))
