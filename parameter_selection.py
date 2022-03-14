"""
This is a demo function for very basic grid search style parameter selction.
By default a 3 X 3 X 3 paramater space is searched and the Sprod performance is 
prioritized based on the qualities of the constructed latent graph. Parameter sets
that preserve the overall spot physical struction and image similarity better will 
bet ranked higher.
The running time of this small search space is a few hours on a normal PC.
You can modify this scipt to your need.

usage: parameter_selection_demo.py [-h] [--input_Rs INPUT_RS]
                                   [--input_K INPUT_K]
                                   [--input_Lambda INPUT_LAMBDA]
                                   [--num_process NUM_PROCESS]
                                   input_path

Demo function for a simple 3 x 3 x 3 pamameter grid search expecting processed
Counts.txt, features and metadata in the input folder.

positional arguments:
  input_path            Input folder containing all necessary files.

optional arguments:
  -h, --help            show this help message and exit
  --input_Rs INPUT_RS, -R INPUT_RS
                        Input Rs, a string separated by ','. (default:
                        0.05,0.1,0.15)
  --input_K INPUT_K, -K INPUT_K
                        Input Ks, a string separated by ','. (default: 3,5,10)
  --input_Lambda INPUT_LAMBDA, -L INPUT_LAMBDA
                        Input RLambdas, a string separated by ','. (default:
                        5,10,15)
  --num_process NUM_PROCESS, -p NUM_PROCESS
                        Number of sprod jobs running at the same time. (default:
                        4)
  --overwrite, -o       Whether to overwrite old outputs with the same
                        parameter combo. (default: False)

"""
import os
import argparse
import pandas as pd
import numpy as np
import random
from multiprocessing import Pool
from scipy.spatial.distance import cdist

def worker(cmd):
    print(cmd)
    os.system(cmd)
    output_path = cmd.split(' ')[3]
    # Remove sprod denoised outputs as those are not needed.
    os.system(
        'rm {}/sprod_Denoised_matrix.txt {}/sprod_Detected_graph.txt {}/sprod_Latent_space.txt'.format(
            output_path, output_path, output_path
        ))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Demo function for a simple 3 x 3 x 3 pamameter grid search \
            expecting processed Counts.txt, features and metadata in the input folder.",
    )

    parser.add_argument(
        "input_path", type=str, help="Input folder containing all necessary files."
    )

    parser.add_argument(
        "output_path", type=str, help="Output folder."
    )

    parser.add_argument(
        "--input_Rs",
        "-R",
        default="0.05,0.1,0.15",
        type=str,
        help="Input Rs, a string separated by ','.",
    )

    parser.add_argument(
        "--input_K",
        "-K",
        default="3,5,10",
        type=str,
        help="Input Ks, a string separated by ','.",
    )

    parser.add_argument(
        "--input_Lambda",
        "-L",
        default="5,10,15",
        type=str,
        help="Input RLambdas, a string separated by ','.",
    )

    parser.add_argument(
        "--num_process",
        "-p",
        default=4,
        type=int,
        help="Number of sprod jobs running at the same time. ",
    )

    parser.add_argument(
        "--overwrite",
        "-o",
        default=False,
        action="store_true",
        help="Whether to overwrite old outputs with the same parameter combo. ",
    )

    args = parser.parse_args()
    input_path = os.path.abspath(args.input_path)
    output_path = os.path.abspath(args.output_path)
    Rs = args.input_Rs.split(',')
    Ks = args.input_K.split(',')
    Ls = args.input_Lambda.split(',')
    num_process = args.num_process
    overwrite = args.overwrite

    # getting script path for supporting codes.
    np.random.seed(0)
    script_path = os.path.abspath(__file__)
    sprod_path = script_path.replace('parameter_selection_demo.py', 'sprod')
    sprodpy_path = script_path.replace('parameter_selection_demo', 'sprod')
    sprodR_path = os.path.join(sprod_path, "denoise.R")
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    os.chdir(output_path)
    # Determine if downsampling is needed
    cts_fn = os.path.join(input_path, "Counts.txt")
    feature_fn = os.path.join(input_path, "spot_level_intensity_features.csv")
    metadata_fn = os.path.join(input_path, "Spot_metadata.csv")
    # IMPORTANT!!!
    # The above three files must have the same order in the rows.
    # If the input is very big, will only run the sprod param demo on a random 
    # 5000 spots
    n_spots = sum(1 for _ in open(cts_fn))
    if n_spots > 10000:
        print('These is more than 10000 spots. Random sample only 5000.')
        subsamples = 5000
        sub_idx = random.sample(range(1,n_spots), subsamples)
        sub_idx.append(0)
        sub_idx = sorted(sub_idx)
        subsampled_input = os.path.join(input_path, 'subsampled')
        if not os.path.exists(subsampled_input):
            os.makedirs(subsampled_input)
        for fn in [cts_fn,feature_fn,metadata_fn]:
            subsampled_lines = []
            with open(fn, 'r') as f:
                line_n = 0
                for line_n, line in enumerate(f):
                    if line_n in sub_idx:
                        subsampled_lines.append(line)
            subsampled_fn = os.path.join(subsampled_input, fn.split('/')[-1])
            with open(subsampled_fn,'w') as s_f:
                for line in subsampled_lines:
                    s_f.write(line)
        input_path = subsampled_input

    # Preparing the Sprod jobs
    if not overwrite:
        print('Will not overwrite existing results.')
    cmd_list = []
    for R in Rs:
        for K in Ks:
            for L in Ls:
                output_path = os.path.join(
                    os.path.abspath('.'), 'R-{}_K-{}_L-{}'.format(R,K,L))
                if not overwrite:
                    if os.path.exists(
                        os.path.join(output_path, 'sprod_log.txt')):
                        print('Skip processing of {}.'.format(output_path))
                        continue
                cmd_list.append(
                    'python {sprodpy_path} {input_path} {output_path} -r {R} -k {K} -l {L} -ws -dg'.format(
                        sprodpy_path = sprodpy_path,
                        input_path = input_path,
                        output_path = output_path,
                        R = R,
                        K = K,
                        L = L,
                    )
                )
    
    with Pool(num_process) as p:
        p.map(worker, cmd_list)

    # processing results, merging PDFs, this process is optionl
    # Will be skipped if the required packages are not found.
    try:
        from PyPDF2 import PdfFileWriter, PdfFileReader, PdfFileMerger
        import io
        from reportlab.pdfgen import canvas
        from reportlab.lib.pagesizes import letter
        for pdf_name in ['Diagnose_Spatial.pdf','Diagnose_TSNE.pdf', 'Diagnose_UMAP.pdf']:
            pdf_list = []
            for dir in os.listdir():
                if not os.path.isdir(dir):
                    continue
                dir_files = os.listdir(dir)
                if pdf_name in dir_files:
                    pdf_list.append(os.path.join(dir, pdf_name))
            # Call the PdfFileMerger
            mergedObject = PdfFileMerger()
            # Loop through all of them and append their pages
            for pdf_fn in sorted(pdf_list):
                packet = io.BytesIO()
                can = canvas.Canvas(packet, pagesize=letter)
                can.drawString(100, 485, pdf_fn)
                can.save()
                #move to the beginning of the StringIO buffer
                packet.seek(0)
                # create a new PDF with Reportlab
                title_pdf = PdfFileReader(packet)
                # read your existing PDF
                existing_pdf = PdfFileReader(open(pdf_fn, "rb"))
                output = PdfFileWriter()
                # add the "watermark" (which is the new pdf) on the existing page
                page = existing_pdf.getPage(0)
                page.mergePage(title_pdf.getPage(0))
                output.addPage(page)
                # finally, write "output" to a real file
                outputStream = open("tmp.pdf", "wb")
                output.write(outputStream)
                outputStream.close()
                mergedObject.append(PdfFileReader("tmp.pdf",'rb'))
            # Write all the files into a file which is named as shown below
            mergedObject.write('Merged ' + pdf_name)
            os.system('rm tmp.pdf')
    except:
        print('PyPDF2 package is not found, abort pdf merging')
    
    # getting R spatial correlations
    sprod_outputs = [x for x in os.listdir('.') if os.path.isdir(x)]
    df_improvements = pd.DataFrame(index = sprod_outputs)
    for dir in sprod_outputs:
        dir_files = os.listdir(dir)
        if 'sprod_log.txt' in dir_files:
            with open(os.path.join(dir,'sprod_log.txt')) as f:
                for line in f:
                    if line.split(':')[0] == '-spatial':
                        spatial_v = float(line.strip('\n').split(':')[1])
                        df_improvements.loc[
                            dir,'Average Graph Distance'] = spatial_v
                    if line.split(':')[0] == '-image tsne':
                        tsne_v = float(line.strip('\n').split(':')[1])
                        df_improvements.loc[
                            dir,'Average Image Distance'] = tsne_v
    print(df_improvements)
    rank1 = pd.Series(
        list(range(1, 1 + df_improvements.shape[0])),
        index = df_improvements.sort_values('Average Graph Distance').index)
    rank2 = pd.Series(
        list(range(1, 1 + df_improvements.shape[0])),
        index = df_improvements.sort_values('Average Image Distance').index)
    df_improvements['Rank'] = (rank1 + rank2)/2
    df_improvements = df_improvements.sort_values('Rank')
    df_improvements.to_csv('pamameter_ranks.csv')