"""
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
                        Input RLambdas, a string separated by ','. (default:
                        4)

"""
import os
import argparse
import pandas as pd
import random
from multiprocessing import Pool
from scipy.spatial.distance import cdist

def worker(cmd):
    print(cmd)
    os.system(cmd)


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
        help="Input RLambdas, a string separated by ','.",
    )

    args = parser.parse_args()
    input_path = os.path.abspath(args.input_path)
    num_process = args.num_process
    Rs = args.input_Rs.split(',')
    Ks = args.input_K.split(',')
    Ls = args.input_Lambda.split(',')

    # getting script path for supporting codes.
    script_path = os.path.abspath(__file__)
    sprod_path = script_path.replace('parameter_selection_demo.py', 'sprod')
    sprodpy_path = script_path.replace('parameter_selection_demo', 'sprod')
    sprodR_path = os.path.join(sprod_path, "denoise.R")

    # Determine if downsampling is needed
    cts_fn = os.path.join(input_path, "Counts.txt")
    feature_fn = os.path.join(input_path, "spot_level_intensity_features.csv")
    metadata_fn = os.path.join(input_path, "Spot_metadata.csv")
    # IMPORTANT!!!
    # The above three files must have the same order in the rows.
    n_spots = sum(1 for _ in open(cts_fn))
    if n_spots > 10000:
        print('These is more than 10000 spots. Random sample only 5000.')
        subsamples = 5000 - 1
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

    cmd_list = []
    for R in Rs:
        for K in Ks:
            for L in Ls:
                output_path = os.path.join(
                    os.path.abspath('.'), 'R-{}_K-{}_L-{}'.format(R,K,L))
                if os.path.exists(
                    os.path.join(output_path, 'Sprod_denoised_matrix.txt')):
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

    # processing results
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
    
    # Processing raw counts
    cpm = pd.read_csv(
        os.path.join(input_path, 'Counts.txt'), sep='\t', index_col=0)
    # raw_cts = pd.read_csv(
    #     os.path.join(input_path, 'Counts.txt'), sep='\t', index_col=0)
    # n_counts = raw_cts.sum(axis=1)
    # cpm = raw_cts.apply(lambda x: 1e6*x/x.sum(), axis=1)
    # cpm = np.log2(1+cpm)
    # rpl_mt_genes = [x for x in cpm.columns if x[:3].lower() in ['rpl', 'rps', 'mt-']]
    # cpm = cpm.drop(rpl_mt_genes, axis=1)
    # Evaluating correlation improvement
    n_top_exp_genes = 25
    n_top_cor_genges = 100
    n_top_corrs = 500
    top_expressed_genes = cpm.median().sort_values(ascending=False)
    top_expressed_genes = top_expressed_genes.index[:n_top_exp_genes]
    corr_matrix = pd.DataFrame(
        1- cdist(
            cpm[top_expressed_genes].T, cpm.drop(top_expressed_genes, axis=1).T,
            metric='correlation'),
        index = top_expressed_genes,
        columns = [x for x in cpm.columns if x not in top_expressed_genes]
    )
    corr_matrix = corr_matrix.dropna(axis=1)
    top_corr = corr_matrix.median().sort_values(ascending=False)[:n_top_cor_genges]
    bot_corr = corr_matrix.median().sort_values(ascending=False)[-n_top_cor_genges:]
    top_corr = top_corr[top_corr>0].index.tolist()
    bot_corr = bot_corr[bot_corr<0].index.tolist()
    top_cor_genes = top_corr + bot_corr
    corr_matrix = corr_matrix[top_cor_genes].stack().reset_index()
    corr_matrix['id'] = corr_matrix.level_0 + '_' + corr_matrix.level_1
    corr_matrix = corr_matrix.set_index('id').iloc[:,2]
    corr_matrix = corr_matrix.sort_values(ascending=False)
    pos_cors = corr_matrix[corr_matrix>0][:n_top_corrs]
    neg_cors = corr_matrix[corr_matrix<0][-n_top_corrs:]

    # get list of denoised counts
    denoised_cts_list = []
    collist = top_expressed_genes.tolist() + top_cor_genes
    collist_R = [x.replace('-','.') for x in collist] # fix for R messing up colnames
    gene_name_mapping = pd.Series(collist, index = collist_R)
    for dir in os.listdir():
        if not os.path.isdir(dir):
            continue
        dir_files = os.listdir(dir)
        if 'sprod_Denoised_matrix.txt' in dir_files:
            denoised_cts_list.append(
                os.path.join(dir, 'sprod_Denoised_matrix.txt'))
    # Calculate correlations
    df_improvements = pd.DataFrame(
        index = [x.split('/')[0] for x in denoised_cts_list],
        columns = ['Pos_cor_improve', 'neg_cor_improve']
    )
    print('Estimating correlation improvement based on top 1000 correlated gene pairs.')
    for denoised_cts_fn in denoised_cts_list:
        print('Processing {}'.format(denoised_cts_fn))
        denoised_cpm = pd.read_csv(
            denoised_cts_fn, usecols=collist_R, sep='\t', index_col=0)
        # denoised_cts = pd.read_csv(
        #     denoised_cts_fn, usecols=collist_R, sep='\t', index_col=0)
        # denoised_cpm = denoised_cts.apply(lambda x: 1e6*x/n_counts[x.name], axis=1)
        # denoised_cpm = np.log2(denoised_cpm+1)
        denoised_cpm.columns = gene_name_mapping[denoised_cpm.columns].values
        denoised_corr = pd.DataFrame(
            1- cdist(
                denoised_cpm[top_expressed_genes].T, denoised_cpm[top_cor_genes].T,
                metric='correlation'),
            index = top_expressed_genes,
            columns = top_cor_genes
            ).stack().reset_index()
        denoised_corr['id'] = denoised_corr.level_0 + '_' + denoised_corr.level_1
        denoised_corr = denoised_corr.set_index('id').iloc[:,2]
        improves = [
            (denoised_corr[pos_cors.index] - pos_cors).mean(),
            (neg_cors - denoised_corr[neg_cors.index]).mean()
        ]
        df_improvements.loc[denoised_cts_fn.split('/')[0]] = improves
    df_improvements.to_csv('Correlation_improvement.csv')
    
    # Plotting
    try:
        import seaborn as sns
        df_improvements['Correlation_improvement'] = df_improvements.sum(axis=1)
        df_improvements = df_improvements.sort_values('Correlation_improvement')
        plt = sns.mpl.pyplot
        _ = plt.figure(figsize=(8,6))
        sns.lineplot(
            data = df_improvements,
            x = [x.replace('_',' ').replace('-',': ') for x in df_improvements.index],
            y='Correlation_improvement')
        _ = plt.xticks(rotation=45, ha = 'right', va = 'top')
        plt.tight_layout()
        plt.savefig('Correlation_improvement.pdf')
        plt.close()
    except:
        print('Seaborn is not found, abort plooting.')
