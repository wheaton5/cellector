# cellector
statistical model for finding anomalous genotype cells in mixed genotype scRNAseq data

static binary for linux x64/x86 included in main directory
python version is now depricated

usage: ./cellector -h
cellector 1.0.0
Haynes Heaton <whheaton@gmail.com>
genotype outlier detection for scRNAseq

USAGE:
    cellector [OPTIONS] --alt <alt> --barcodes <barcodes> --output_directory <output_directory> --ref <ref>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -a, --alt <alt>                                                      alt.mtx matrix from vartrix
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

    -r, --ref <ref>                                                      ref.mtx matrix from vartrix
    -v, --vcf <vcf>
            vcf associated with the alt.mtx and ref.mtx, only required to associated loci with variants and to get the
            genotypes of the minority and majority cells in a vcf output
