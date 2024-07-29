# cellector
statistical model for finding anomalous genotype cells in mixed genotype scRNAseq data

From cellranger outputs (bam, cell barcodes.tsv) you can use cellector_pipeline.py with the following usage.
```
./cellector_pipeline.py -h
usage: cellector_pipeline.py [-h] -i BAM -b BARCODES -f FASTA -t THREADS -o OUT_DIR --common_variants COMMON_VARIANTS
                             [--min_alt MIN_ALT] [--min_ref MIN_REF] [--ignore IGNORE] [--cellector_binary CELLECTOR_BINARY]
                             [--souporcell_binary SOUPORCELL_BINARY]

single cell RNAseq foreign genotype cell detection

optional arguments:
  -h, --help            show this help message and exit
  -i BAM, --bam BAM     cellranger bam
  -b BARCODES, --barcodes BARCODES
                        barcodes.tsv from cellranger
  -f FASTA, --fasta FASTA
                        reference fasta file
  -t THREADS, --threads THREADS
                        max threads to use
  -o OUT_DIR, --out_dir OUT_DIR
                        name of directory to place cellector files
  --common_variants COMMON_VARIANTS
                        common variant loci or known variant loci vcf, must be vs same reference fasta
  --min_alt MIN_ALT     min alt to use locus, default = 4.
  --min_ref MIN_REF     min ref to use locus, default = 4.
  --ignore IGNORE       set to True to ignore data error assertions
  --cellector_binary CELLECTOR_BINARY
                        /path/to/cellector
  --souporcell_binary SOUPORCELL_BINARY
                        /path/to/souporcell
```

If you have the alt.mtx and ref.mtx you can use cellector directly with the following usage.

static binary for linux x64/x86 included in main directory
python version is now depricated
```
usage: ./cellector_linux -h
cellector 1.0.0
Haynes Heaton <whheaton@gmail.com>
genotype outlier detection for scRNAseq

USAGE:
    cellector [OPTIONS] --alt <alt> --ref <ref> --barcodes <barcodes> --output_directory <output_directory> 

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -r, --ref <ref>                                                      ref.mtx matrix from vartrix
    -a, --alt <alt>                                                      alt.mtx matrix from vartrix
    -v, --vcf <vcf>
            vcf associated with the alt.mtx and ref.mtx, only required to associated loci with variants and to get the
            genotypes of the minority and majority cells in a vcf output
    -b, --barcodes <barcodes>                                            cell barcodes
    -g, --ground_truth <ground_truth> 
            cell hashing assignments or other ground truth, first column barcodes second column assignment

        --interquartile_range_multiple <interquartile_range_multiple>
            number of interquartile range multiples away from 25th percentile to make the threshold to call an outline

        --min_alleles_posterior <min_alleles_posterior>
            minimum number of alleles for both minority distribution and majority distribution to have for a locus to be
            used for posterior probability calculation
        --min_alt <min_alt>
            minimum number of cells containing the alt allele for the variant to be used for clustering (default 4)

        --min_ref <min_ref>
            minimum number of cells containing the ref allele for the variant to be used for clustering (default 4)

        --output_directory <output_directory>
            output directory where other output files will be stored

        --posterior_threshold <posterior_threshold>
            posterior probability threshold for assignment of minority or majority (default 0.999)

```
