# cellector
statistical model for finding anomalous genotype cells in mixed genotype scRNAseq data

usage: cellector.py [-h] -a ALT -r REF [--min_ref MIN_REF] [--min_alt MIN_ALT] --barcodes BARCODES [--ground_truth GROUND_TRUTH] --output_prefix OUTPUT_PREFIX

find outlier genotype cells in single cell experiment

optional arguments:
  -h, --help            show this help message and exit
  -a ALT, --alt ALT     alt.mtx from vartrix
  -r REF, --ref REF     ref.mtx from vartrix
  --min_ref MIN_REF     minimum ref count to use loci
  --min_alt MIN_ALT     minimum alt count to use loci
  --barcodes BARCODES   barcodes.tsv file
  --ground_truth GROUND_TRUTH (optional)
                        cell hashing assignments to barcodes or similar, two columns first column barcode, second column ground truth assignment
  --output_prefix OUTPUT_PREFIX
                        output prefix




