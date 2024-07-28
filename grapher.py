# pyhton wrapper for the rust file

import argparse
import glob
import pandas as pd
from plotnine import *


# TODO put it in a function
parser = argparse.ArgumentParser(description="creates graphs from cellector data")
parser.add_argument("-d","--data_directory", required=True, help="directory cellector output files are stored")
args = parser.parse_args()

globed_files = glob.glob(args.data_directory + "/iteration_?.tsv")

for file_name in globed_files:
    with open(file_name, mode='r') as fid:
        df = pd.read_csv(fid, sep="\t")
        df['neg_log_likelihood_loci_normalized'] = -df['log_likelihood'] / df['num_loci_used']

        sorted_loci_norm = sorted(-df['neg_log_likelihood_loci_normalized'])
        median = sorted_loci_norm[len(sorted_loci_norm)//2]
        q1 = sorted_loci_norm[len(sorted_loci_norm)//4]
        q3 = sorted_loci_norm[int(len(sorted_loci_norm)*0.75)]
        interquartile_range = q3 - q1
        threshold = q1 - 5*interquartile_range

        graph = ggplot(df)+geom_point(aes(x="cell_id",y="neg_log_likelihood_loci_normalized",size="num_loci_used",color="assignment"))+geom_abline(intercept=-threshold,slope=0)
        graphname = file_name[:-4]+".pdf"
        save_as_pdf_pages([graph],filename=graphname)



