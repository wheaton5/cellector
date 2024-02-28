#!/usr/bin/env python

from scipy.stats import betabinom
import argparse
from plotnine import *
import numpy as np

parser = argparse.ArgumentParser(description="find outlier genotype cells in single cell experiment")
parser.add_argument("-a","--alt", required=True, help="alt.mtx from vartrix")
parser.add_argument("-r","--ref", required=True, help="ref.mtx from vartrix")
parser.add_argument("--min_ref", required=False, default = 10, type=int, help="minimum ref count to use loci")
parser.add_argument("--min_alt", required=False, default = 10, type=int, help="minimum alt count to use loci")
parser.add_argument("--barcodes", required=True, help="barcodes.tsv file")
parser.add_argument("--ground_truth", required=False, help="cell hashing assignments to barcodes or similar, two columns first column barcode, second column ground truth assignment")
args = parser.parse_args()

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

cell_data = {} # map from cell id to map from loci to counts
loci_counts = {}
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

cell_log_likelihoods = {} # map from cell id to [log_likelihood, #loci_used, #alleles_used]
loci_normalized = []
alleles_normalized = []
print("\t".join(["cell_idx","barcode","assignment","log_likelihood","num_loci_used","num_alleles_used","log_likelihood_loci_normalized","log_likelihood_alleles_normalized"]))
for cell_id in cell_data.keys():
    log_likelihood = 0
    num_loci_used = 0
    num_alleles_used = 0
    cell_counts = cell_data[cell_id]
    for locus in cell_counts.keys():
        if locus in loci_used:
            locus_counts = loci_counts[locus]
            alpha = locus_counts[1] + 1 # bayes prior alpha 1 and beta 1
            beta = locus_counts[0] + 1
            cell_locus_counts =  cell_counts[locus]
            num_alt = cell_locus_counts[1]
            num_ref = cell_locus_counts[0]
            num_loci_used += 1
            total_alleles = num_alt + num_ref
            num_alleles_used += total_alleles
            #dist = betabinom(total_alleles, alpha, beta)
            log_likelihood += betabinom.logpmf(num_alt, total_alleles, alpha, beta)
    cell_log_likelihoods[cell_id] = [log_likelihood, num_loci_used, num_alleles_used]
    assignment = "na"
    if cell_id in cell_id_assignment:
        assignment = cell_id_assignment[cell_id]
    if num_loci_used > 0 and num_alleles_used > 0:
        loci_normalized.append(log_likelihood / num_loci_used)
        alleles_normalized.append(log_likelihood / num_alleles_used)
        print("\t".join([str(cell_id),str(cell_id_to_barcode[cell_id]), assignment, str(log_likelihood), str(num_loci_used), str(num_alleles_used), str(log_likelihood/num_loci_used), str(log_likelihood/num_alleles_used)]))
    else:
        print("\t".join([str(cell_id),str(cell_id_to_barcode[cell_id]), assignment, str(log_likelihood), str(num_loci_used), str(num_alleles_used), "0", "0"]))

mean_loci_normalized = np.mean(loci_normalized)
sd_loci_normalized = np.std(loci_normalized)

print("loci normalized mean=",mean_loci_normalized, " std=",sd_loci_normalized, " 3 sigma threshold=",mean_loci_normalized-3*sd_loci_normalized)

mean_alleles_normalized = np.mean(alleles_normalized)
sd_alleles_normalized = np.std(alleles_normalized)

print("alleles normalized mean=",mean_alleles_normalized, " std=",sd_alleles_normalized, " 3 sigma threshold=",mean_alleles_normalized-3*sd_alleles_normalized)
