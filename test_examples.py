# Need to install sprod first
# Need to do module purge so that the R.cpp will work
# Need R/4.0.2

import pandas as pd
import os


sprod_path = os.path.abspath(__file__)
denoise_jb = sprod_path.replace('test_examples.py', 'sprod/denoise_job.py')
input_path = sprod_path.replace('test_examples.py', 'example/input')
test_path = sprod_path.replace('test_examples.py', 'example/test')
output_path = sprod_path.replace('test_examples.py', 'example/output')

os.system(
    'python {} -y patches -pn 2 -pb 2 {} {}/batch_with_img'.format(denoise_jb, input_path, test_path))
try:
    test_denoised = pd.read_hdf(os.path.join(test_path,'batch_with_img/denoised_stiched.hdf'))
    ref_denoised = pd.read_hdf(os.path.join(output_path,'batch_with_img/denoised_stiched.hdf'))
    if (test_denoised - ref_denoised).sum().sum() <= 100:
        print('Sprod test patches mode succeeded!')
    else:
        print('Sprod test patches mode was able to run, but the results were not expected. Please check if the input files match those in the repo.')
except FileNotFoundError:
    print('Sprod denoise job failed: no output found, please check sprod_log.txt.')

os.system(
    'python {} {} {}/single_with_img'.format(denoise_jb, input_path, test_path))
try:
    test_denoised = pd.read_csv(
        os.path.join(test_path,'single_with_img/sprod_Denoised_matrix.txt'), 
        sep='\t')
    ref_denoised = pd.read_csv(
        os.path.join(output_path,'single_with_img/sprod_Denoised_matrix.txt'), 
        sep='\t')
    if (test_denoised - ref_denoised).sum().sum() <= 100:
        print('Sprod test single mode succeeded!')
    else:
        print('Sprod test single mode was able to run, but the results were not expected. Please check if the input files match those in the repo.')
except FileNotFoundError:
    print('Sprod denoise job failed: no output found, please check sprod_log.txt.')