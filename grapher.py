# pyhton wrapper for the rust file

import argparse
import glob
import pandas as pd
from plotnine import *


# TODO put it in a function
parser = argparse.ArgumentParser(description="creates graphs from cellector data")
parser.add_argument("-d","--data_directory", required=True, help="directory cellector output files are stored")
args = parser.parse_args()

data_files = sorted(glob.glob(args.data_directory + "/iteration_?.tsv"))
threshold_files = sorted(glob.glob(args.data_directory + "/iteration_?_threshold.tsv"))

data_file_indices = [int(name.split("_")[-1][:-4]) for name in data_files]
threshold_file_indices = [int(name.split("_")[-2]) for name in threshold_files]
if(any([x != y for x, y in zip(data_file_indices, threshold_file_indices)])):
    print("iteration and threshold file indices do not match")
else:
    for data_file_name, threshold_file_name in zip(data_files, threshold_files):
        with open(data_file_name, mode='r') as data_fid, open(threshold_file_name, mode="r") as threshold_fid:
            df = pd.read_csv(data_fid, sep="\t")
            df['neg_log_likelihood_loci_normalized'] = -df['log_likelihood'] / df['num_loci_used']

            sorted_loci_norm = sorted(-df['neg_log_likelihood_loci_normalized'])
            median = sorted_loci_norm[len(sorted_loci_norm)//2]
            q1 = sorted_loci_norm[len(sorted_loci_norm)//4]
            q3 = sorted_loci_norm[int(len(sorted_loci_norm)*0.75)]
            interquartile_range = q3 - q1
            #threshold = q1 - 4*interquartile_range
            threshold = float(threshold_fid.readline())

            graph = ggplot(df)+geom_point(aes(x="cell_id",y="neg_log_likelihood_loci_normalized",size="num_loci_used",color="assignment"))+geom_abline(intercept=-threshold,slope=0)
            graphname = data_file_name[:-4]+".pdf"
            save_as_pdf_pages([graph],filename=graphname)


