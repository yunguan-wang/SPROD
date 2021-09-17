# Need to install sprod first
# Need to do module purge so that the R.cpp will work
# Need R/4.0.2

python ./sprod/denoise_job.py -y patches -pn 2 -pb 2 ./example/input ./example/output/batch_with_img
python ./sprod/denoise_job.py ./example/input ./example/output/single_with_img