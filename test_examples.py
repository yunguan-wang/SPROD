# Need to install sprod first
# Need to do module purge so that the R.cpp will work
# Need R/4.0.2

import pandas as pd
import os
import platform

os_type = platform.system()
sprod_path = os.path.abspath(__file__)
denoise_jb = sprod_path.replace('test_examples.py', 'sprod.py')
input_path = sprod_path.replace('test_examples.py', 'test_example/input')
test_path = sprod_path.replace('test_examples.py', 'test_example/test_output')
output_path = sprod_path.replace('test_examples.py', 'test_example/expected_output')
sprod_path = sprod_path + '/sprod'

print('Removing previous testing results...', end='')
if os_type == 'Windows':
    os.system('del -r ./{}'.format(test_path))
else:
    os.system('rm -rf {}/'.format(test_path))
print('Done!')

print('Testing Sprod in single mode')
os.system(
    'python {} {} {}/single_with_img -ci if'.format(denoise_jb, input_path, test_path))
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

print('Testing Sprod in batch mode')
os.system(
    'python {} -y batch -pb 2 -ci if {} {}/batch_with_img'.format(denoise_jb, input_path, test_path))
try:
    test_denoised = pd.read_csv(os.path.join(test_path,'batch_with_img/denoised_stitched.txt'), sep='\t')
    ref_denoised = pd.read_csv(os.path.join(output_path,'batch_with_img/denoised_stitched.txt'), sep='\t')
    if (test_denoised - ref_denoised).sum().sum() <= 100:
        print('Sprod test batch mode succeeded!')
    else:
        print('Sprod test batch mode was able to run, but the results were not expected. Please check if the input files match those in the repo.')
except FileNotFoundError:
    print('Sprod denoise job failed: no output found, please check sprod_log.txt.')


print('Testing Sprod in single, pseudoimage mode')
tif_fn = os.path.join(input_path, 'image.tif')
os.rename(tif_fn, tif_fn.replace('.tif',''))
try:
    os.system(
        'python {} {} {}/single_pseudoimage'.format(denoise_jb, input_path, test_path))
    os.rename(tif_fn.replace('.tif',''), tif_fn)
except:
    os.rename(tif_fn.replace('.tif',''), tif_fn)

try:
    test_denoised = pd.read_csv(
        os.path.join(test_path,'single_pseudoimage/sprod_Denoised_matrix.txt'), 
        sep='\t')
    ref_denoised = pd.read_csv(
        os.path.join(output_path,'single_pseudoimage/sprod_Denoised_matrix.txt'), 
        sep='\t')
    if (test_denoised - ref_denoised).sum().sum() <= 100:
        print('Sprod test single pseudoimage mode succeeded!')
    else:
        print('Sprod test single pseudoimage mode was able to run, but the results were not expected. Please check if the input files match those in the repo.')
except FileNotFoundError:
    print('Sprod denoise job failed: no output found, please check sprod_log.txt.')
    
# clean up intermediate files in the input folder.
os.system('rm -rf {}/*features.csv'.format(input_path))