#!/usr/bin/env python

from scipy.stats import betabinom
import argparse
from plotnine import *
import numpy as np
import subprocess
import pandas as pd
import os

parser = argparse.ArgumentParser(description="find outlier genotype cells in single cell experiment")
parser.add_argument("-a","--alt", required=True, help="alt.mtx from vartrix")
parser.add_argument("-r","--ref", required=True, help="ref.mtx from vartrix")
parser.add_argument("--min_ref", required=False, default = 10, type=int, help="minimum ref count to use loci")
parser.add_argument("--min_alt", required=False, default = 10, type=int, help="minimum alt count to use loci")
parser.add_argument("--barcodes", required=True, help="barcodes.tsv file")
parser.add_argument("--ground_truth", required=False, help="cell hashing assignments to barcodes or similar, two columns first column barcode, second column ground truth assignment")
parser.add_argument("--output_prefix", required=True, help="output prefix")
args = parser.parse_args()

if os.path.isdir(args.output_prefix):
    print("restarting pipeline in existing directory? not implemented yet " + args.output_prefix)
else:
    subprocess.check_call(["mkdir", "-p", args.output_prefix])

barcode_assignment = {} # map from barcode to assignment
cell_id_assignment = {}
cell_id_to_barcode = {}
barcode_to_cell_id = {}
with open(args.barcodes) as infile:
    for (cell_id_minus_one, line) in enumerate(infile):
        cell_id = cell_id_minus_one + 1
        bc = line.strip()
        cell_id_to_barcode[cell_id] = bc
        barcode_to_cell_id[bc] = cell_id

if args.ground_truth:
    with open(args.ground_truth) as infile:
        for line in infile:
            toks = line.strip().split()
            bc = toks[0]
            assignment = toks[1]
            barcode_assignment[bc] = assignment
            cell_id_assignment[barcode_to_cell_id[bc]] = assignment


loci_used = set() # is this gonna be a set or list?

cell_data = {} # map from cell id to map from loci to ref,alt counts
loci_counts = {} # map from locus to ref,alt counts

with open(args.alt) as altfile:
    altfile.readline()
    altfile.readline()
    altfile.readline()
    for line in altfile:
        toks = line.strip().split()
        locus = int(toks[0])
        cell = int(toks[1])
        count = int(toks[2])
        counts = loci_counts.setdefault(locus, [0,0])
        counts[1] += count
        cell_counts = cell_data.setdefault(cell, {})
        cell_loci_counts = cell_counts.setdefault(locus,[0,0])
        cell_loci_counts[1] += count



with open(args.ref) as reffile:
    reffile.readline()
    reffile.readline()
    reffile.readline()
    for line in reffile:
        toks = line.strip().split()
        locus = int(toks[0])
        cell = int(toks[1])
        count = int(toks[2])
        counts = loci_counts.setdefault(locus, [0,0])
        counts[0] += count
        cell_counts = cell_data.setdefault(cell, {})
        cell_loci_counts = cell_counts.setdefault(locus,[0,0])
        cell_loci_counts[0] += count

for locus in loci_counts.keys():
    counts = loci_counts[locus]
    if counts[0] >= args.min_ref and counts[1] >= args.min_alt:
        loci_used.add(locus)


loci_used_sorted = sorted(list(loci_used))


cells_to_exclude = set()
iteration = 0
any_change = True
while any_change:
    any_change = False
    cell_log_likelihoods = {} # map from cell id to [log_likelihood, #loci_used, #alleles_used]
    loci_normalized = []
    alleles_normalized = []
    loci_cell_likelihoods = {} # map from locus to cell to likelihood
    beta_alpha = {}
    for locus in loci_used:
        locus_counts = loci_counts[locus]
        alpha = locus_counts[1] + 1 # bayes prior alpha 1 and beta 1
        beta = locus_counts[0] + 1
        beta_alpha[locus] = [beta, alpha]
    for cell_id in cells_to_exclude:
        cell_counts = cell_data[cell_id]
        for locus in cell_counts.keys():
            if locus in loci_used:
                cell_locus_counts = cell_counts[locus]
                num_alt = cell_locus_counts[1]
                num_ref = cell_locus_counts[0]
                beta_alpha[locus][0] -= num_ref
                beta_alpha[locus][1] -= num_alt
                 
    filename = args.output_prefix+"/iteration_"+str(iteration)+".tsv"
    with open(filename, 'w') as out:
        out.write("\t".join(["cell_idx","barcode","assignment","log_likelihood","num_loci_used","num_alleles_used","neg_log_likelihood_loci_normalized","neg_log_likelihood_alleles_normalized"])+"\n")
        for cell_id in cell_data.keys():
            log_likelihood = 0
            num_loci_used = 0
            num_alleles_used = 0
            cell_counts = cell_data[cell_id]
            for locus in cell_counts.keys():
                if locus in loci_used:
                    beta, alpha = beta_alpha[locus]
                    cell_locus_counts =  cell_counts[locus]
                    num_alt = cell_locus_counts[1]
                    num_ref = cell_locus_counts[0]
                    num_loci_used += 1
                    total_alleles = num_alt + num_ref
                    num_alleles_used += total_alleles
                    #dist = betabinom(total_alleles, alpha, beta)
                    logpmf = betabinom.logpmf(num_alt, total_alleles, alpha, beta)
                    log_likelihood += logpmf
                    cell_likelihood = loci_cell_likelihoods.setdefault(locus, {})
                    cell_likelihood[cell_id] = logpmf
            if num_loci_used > 0:
                cell_log_likelihoods[cell_id] = log_likelihood/num_loci_used
            else:
                cell_log_likelihoods[cell_id] = 0.0
            assignment = "na"
            if cell_id in cell_id_assignment:
                assignment = cell_id_assignment[cell_id]
            if num_loci_used > 0 and num_alleles_used > 0:
                loci_normalized.append(log_likelihood / num_loci_used)
                alleles_normalized.append(log_likelihood / num_alleles_used)
                out.write("\t".join([str(cell_id),str(cell_id_to_barcode[cell_id]), assignment, str(log_likelihood), str(num_loci_used), str(num_alleles_used), str(-log_likelihood/num_loci_used), str(-log_likelihood/num_alleles_used)])+"\n")
            else:
                out.write("\t".join([str(cell_id),str(cell_id_to_barcode[cell_id]), assignment, str(log_likelihood), str(num_loci_used), str(num_alleles_used), "0", "0"])+"\n")

        mean_loci_normalized = np.mean(loci_normalized)
        sd_loci_normalized = np.std(loci_normalized)

        sorted_loci_norm = sorted(loci_normalized)
        median = sorted_loci_norm[len(loci_normalized)//2]
        q1 = sorted_loci_norm[len(loci_normalized)//4]
        q3 = sorted_loci_norm[int(len(loci_normalized)*0.75)]
        interquartile_range = q3 - q1
        threshold = q1 - 4*interquartile_range

        print("loci normalized median=",median, " iqr=",interquartile_range, " q1-4*iqr=",threshold)        
        print("loci normalized mean=",mean_loci_normalized, " std=",sd_loci_normalized, " 3 sigma threshold=",mean_loci_normalized-3*sd_loci_normalized)
        df = pd.read_csv(filename,sep="\t")
        graphname = filename[:-4]+".pdf"
        print(graphname)
        graph = ggplot(df)+geom_point(aes(x="cell_idx",y="neg_log_likelihood_loci_normalized",size="num_loci_used",color="assignment"))+geom_abline(intercept=-threshold,slope=0)
        save_as_pdf_pages([graph],filename=graphname)

        #threshold = mean_loci_normalized - 3*sd_loci_normalized
        removed_cells = 0
        for (cell_id, log_likelihood) in cell_log_likelihoods.items():
            if log_likelihood < threshold:
                if not(cell_id in cells_to_exclude):
                    removed_cells += 1
                    any_change = True
                    cells_to_exclude.add(cell_id)
        print("removed ",removed_cells," in the ",iteration," iteration")
        # lets find which loci contribute the most to exclusion
        exclusion_locus_loglikes = {} #map from locus to sum of loglikelihoods from excluded cells
        exclusion_locus_cellcounts = {}
        exclusion_locus_alts = {}
        exclusion_locus_refs = {}
        majority_locus_alts = {}
        majority_locus_refs = {}
        for (locus, cell_likelihoods) in loci_cell_likelihoods.items():
            exclusion_locus_loglikes.setdefault(locus, 0)
            exclusion_locus_cellcounts.setdefault(locus, 0)
            exclusion_locus_alts.setdefault(locus,0)
            exclusion_locus_refs.setdefault(locus,0)
            majority_locus_alts.setdefault(locus,0)
            majority_locus_refs.setdefault(locus,0)
            for (cell, log_likelihood) in cell_likelihoods.items():
                cell_counts = cell_data[cell]
                if locus in cell_counts:
                    cell_locus_counts =  cell_counts[locus]
                else:
                    cell_locus_counts = [0,0]
                num_alt = cell_locus_counts[1]
                num_ref = cell_locus_counts[0]
                if cell in cells_to_exclude:
                    exclusion_locus_loglikes[locus] += log_likelihood
                    exclusion_locus_cellcounts[locus] += 1
                    exclusion_locus_alts[locus] += num_alt
                    exclusion_locus_refs[locus] += num_ref    
                else:
                    majority_locus_alts[locus] += num_alt
                    majority_locus_refs[locus] += num_ref
        sorted_locus_loglikes = sorted(exclusion_locus_loglikes.items(), key=lambda x:x[1])
        outfile = filename[:-4]+"_locus_contribution_minority.tsv"
        with open(outfile,'w') as out:
            out.write("\t".join(["locus","log_likelihood_minority","minority_cellcount", "log_likelihood_minority_per_cell","minority_alt","minority_ref","majority_alt","majority_ref","minority_af","majority_af"])+"\n")
            for (locus, loglike) in sorted_locus_loglikes:
                cellcounts = exclusion_locus_cellcounts[locus]
                loglikepercell = 0
                if cellcounts > 0:
                    loglikepercell = loglike/cellcounts
                minority_alt = exclusion_locus_alts[locus]
                minority_ref = exclusion_locus_refs[locus]
                minority_af = 0
                if minority_alt + minority_ref > 0:
                    minority_af = minority_alt/(minority_alt+minority_ref)
                majority_alt = majority_locus_alts[locus]
                majority_ref = majority_locus_refs[locus]
                majority_af = 0
                if majority_alt + majority_ref > 0:
                    majority_af = majority_alt/(majority_alt+majority_ref)
                out.write("\t".join([str(locus), str(loglike), str(exclusion_locus_cellcounts[locus]), str(loglikepercell),str(minority_alt),str(minority_ref),str(majority_alt),str(majority_ref),str(minority_af),str(majority_af)])+"\n")
        iteration += 1
with open(args.output_prefix+"/cellector_assignments.tsv",'w') as out:
    out.write("\t".join(["barcode","cellector_assignment","ground_truth_assignment"])+"\n")
    with open(args.barcodes) as infile:
        for line in infile:
            bc = line.strip()
            cell_id = barcode_to_cell_id[bc]
            assignment = "majority"
            if cell_id in cells_to_exclude:
                assignment = "minority"
            gt = "na"
            if cell_id in cell_id_assignment:
                gt = cell_id_assignment[cell_id]
        
            out.write("\t".join([bc,assignment,gt])+"\n")
