#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(
    description="single cell RNAseq foreign genotype cell detection")
parser.add_argument("-i", "--bam", required = True, help = "cellranger bam")
parser.add_argument("-b", "--barcodes", required = True, help = "barcodes.tsv from cellranger")
parser.add_argument("-f", "--fasta", required = True, help = "reference fasta file")
parser.add_argument("-t", "--threads", required = True, type = int, help = "max threads to use")
parser.add_argument("-o", "--out_dir", required = True, help = "name of directory to place cellector files")
parser.add_argument("--common_variants", required = True, default = None, 
    help = "common variant loci or known variant loci vcf, must be vs same reference fasta")
parser.add_argument("--min_alt", required = False, default = "4", help = "min alt to use locus, default = 4.")
parser.add_argument("--min_ref", required = False, default = "4", help = "min ref to use locus, default = 4.")
parser.add_argument("--ignore", required = False, default = False, type = bool, help = "set to True to ignore data error assertions")
parser.add_argument("--cellector_binary", required=False, default = "cellector_linux", help = "/path/to/cellector")
parser.add_argument("--souporcell_binary", required=False, default = "souporcell_linux", help="/path/to/souporcell")
parser.add_argument("--troublet_binary", required=False, default = "troublet_linux", help="/path/to/troublet")
parser.add_argument("--grapher_script", required=False, default = "grapher.py", help="/path/to/grapher.py")
args = parser.parse_args()

print("checking modules")
# importing all reqs to make sure things are installed
import numpy as np
import gzip
import math
import pysam
import pyfaidx
import subprocess
import time
import os
print("imports done")

open_function = lambda f: gzip.open(f,"rt") if f[-3:] == ".gz" else open(f)

UMI_TAG = "UB"
CELL_TAG = "CB"

#load each file to make sure it is legit
bc_set = set()
with open_function(args.barcodes) as barcodes:
    for (index, line) in enumerate(barcodes):
        bc = line.strip()
        bc_set.add(bc)

assert len(bc_set) > 50, "Fewer than 50 barcodes in barcodes file? We expect 1 barcode per line."

#test bam load
bam = pysam.AlignmentFile(args.bam)
num_cb = 0
num_cb_cb = 0 # num reads with barcodes from barcodes.tsv file
num_umi = 0
num_read_test = 100000
for (index,read) in enumerate(bam):
    if index >= num_read_test:
        break
    if read.has_tag(CELL_TAG):
        num_cb += 1
        if read.get_tag(CELL_TAG) in bc_set:
            num_cb_cb += 1
    if read.has_tag(UMI_TAG):
        num_umi += 1
if not args.ignore:
    assert float(num_cb) / float(num_read_test) > 0.5, "Less than 50% of first 100000 reads have cell barcode tag (CB), turn on --ignore True to ignore"
    assert float(num_umi) / float(num_read_test) > 0.5, "Less than 50% of first 100000 reads have UMI tag (UB), turn on --ignore True to ignore"
    assert float(num_cb_cb) / float(num_read_test) > 0.05, "Less than 25% of first 100000 reads have cell barcodes from barcodes file, is this the correct barcode file? turn on --ignore True to ignore"

print("checking fasta")
fasta = pyfaidx.Fasta(args.fasta, key_function = lambda key: key.split()[0])

def get_bam_regions(bamname, threads):
    bam = pysam.AlignmentFile(bamname)
    total_reference_length = 0
    for chrom in bam.references:
        total_reference_length += bam.get_reference_length(chrom)
    step_length = int(math.ceil(total_reference_length / threads))
    regions = []
    region = []
    region_so_far = 0
    chrom_so_far = 0
    for chrom in bam.references:
        chrom_length = bam.get_reference_length(chrom)
        #print(chrom+" size "+str(chrom_length)+" and step size "+str(step_length))
        while True:
            #print("\tregion so far "+str(region_so_far)+" chrom so far "+str(chrom_so_far)) 
            #print("\t"+str(chrom_length - chrom_so_far)+" <= "+str(step_length - region_so_far))
            #print("\t"+str((chrom_length - chrom_so_far) <= step_length - region_so_far))
            if (chrom_length - chrom_so_far) <= step_length - region_so_far:
                region.append((chrom, chrom_so_far, chrom_length))
                #print("\t\tending chrom\t"+chrom+":"+str(chrom_so_far)+"-"+str(chrom_length))
                region_so_far += chrom_length - chrom_so_far
                chrom_so_far = 0
                break
            else:
                region.append((chrom, chrom_so_far, chrom_so_far + step_length - region_so_far))
                #print("\t\tending region\t"+chrom+":"+str(chrom_so_far)+"-"+str(chrom_so_far + step_length - region_so_far))
                regions.append(region)
                region = []
                chrom_so_far += step_length - region_so_far
                region_so_far = 0
    if len(region) > 0:
        regions.append(region)
    
    return regions

def freebayes(args, bam, fasta):
        print("using common variants")
        
        # parallelize the samtools depth call. It takes too long
        regions = get_bam_regions(bam, int(args.threads))
        depth_files = []
        depth_procs = []
        print(len(regions))
        for (index, region) in enumerate(regions):
            region_args = []
            for (chrom, start, stop) in region:
                region_args.append(chrom+":"+str(start)+"-"+str(stop))
            depthfile = args.out_dir+"/depth_"+str(index)+".bed"
            depth_files.append(depthfile)
            min_cov = int(args.min_ref)+int(args.min_alt)
            with open(depthfile, 'w') as bed:
                with open(depthfile+".sh",'w') as depther:
                    depther.write("samtools view -hb "+bam+" "+" ".join(region_args)+ " | samtools depth - | "+
                    "awk '{ if ($3 >= "+str(min_cov)+ " && $3 < 100000) { print $1 \"\t\" $2 \"\t\" $2+1 \"\t\" $3 } }'")
                subprocess.check_call(["chmod", "777", depthfile+".sh"])
                #ps0 = subprocess.Popen(['samtools', 'view', bam]+region_args, stdout = subprocess.PIPE)
                #ps1 = subprocess.Popen(['samtools', 'depth', '-'], stdin = ps0.stdout, stdout = subprocess.PIPE)
                # awk magic 
                #ps2 = subprocess.Popen(["awk '{ if ($3 >= " + str(min_cov) + " && $3 < 100000) { print $1 \"\t\" $2 \"\t\" $2+1 \"\t\" $3 } }'"], 
                #    shell = True, stdin = ps1.stdout, stdout = bed)
                ps = subprocess.Popen([depthfile+".sh"], shell = True, stdout = bed)
                depth_procs.append(ps)

        for proc in depth_procs:
            proc.wait()
        merged_depthfiles = []
        for depth_file in depth_files:
            merged_depthfile = depth_file[:-4]+"_merged.bed"
            with open(merged_depthfile, 'w') as bed:
                subprocess.check_call(["bedtools", "merge", "-i", depth_file], stdout = bed)
            merged_depthfiles.append(merged_depthfile)
        with open(args.out_dir + "/depth_merged.bed", 'w') as merged_bed:
            subprocess.check_call(['cat']+merged_depthfiles, stdout = merged_bed)
        for tmp in depth_files: # clean up tmp bed files
            subprocess.check_call(['rm', tmp, tmp+".sh"])
        for tmp in merged_depthfiles:
            subprocess.check_call(['rm', tmp])


        with open(args.out_dir + "/common_variants_covered_tmp.vcf", 'w') as vcf:
            subprocess.check_call(["bedtools", "intersect", "-wa", "-a", args.common_variants, "-b", args.out_dir + "/depth_merged.bed"], stdout = vcf)
        with open(args.out_dir + "/common_variants_covered_tmp.vcf") as vcf:
            with open(args.common_variants) as common:
                with open(args.out_dir + "/common_variants_covered.vcf",'w') as out:
                    for line in common:
                        if line.startswith("#"):
                            out.write(line)
                        else:
                            break
                    for line in vcf:
                        out.write(line)
        with open(args.out_dir + "/variants.done", 'w') as done:
            done.write(args.out_dir + "/common_variants_covered.vcf" + "\n")
        return(args.out_dir + "/common_variants_covered.vcf")

def vartrix(args, final_vcf, final_bam):
    print("running vartrix")
    ref_mtx = args.out_dir + "/ref.mtx"
    alt_mtx = args.out_dir + "/alt.mtx"  
    barcodes = args.barcodes
    if barcodes[-3:] == ".gz":
        with open(args.out_dir + "/barcodes.tsv",'w') as bcsout:
            subprocess.check_call(['gunzip', '-c', barcodes],stdout = bcsout)
        barcodes = args.out_dir + "/barcodes.tsv"
    with open(args.out_dir + "/vartrix.err", 'w') as err:
        with open(args.out_dir + "/vartrix.out", 'w') as out:
            cmd = ["vartrix", "--mapq", "30", "-b", final_bam, "-c", barcodes, "--scoring-method", "coverage", "--threads", str(args.threads),
                "--ref-matrix", ref_mtx, "--out-matrix", alt_mtx, "-v", final_vcf, "--fasta", args.fasta]
            cmd.append("--umi")
            subprocess.check_call(cmd, stdout = out, stderr = err)
    subprocess.check_call(['touch', args.out_dir + "/vartrix.done"])
    subprocess.check_call(['rm', args.out_dir + "/vartrix.out", args.out_dir + "/vartrix.err"])
    return((ref_mtx, alt_mtx))

def read_tsv(filename):
    with open(filename,'r') as fid:
        table = []
        for line in fid.readlines():
            table.append(line.strip().split("\t"))
    return table


#### MAIN RUN SCRIPT
if os.path.isdir(args.out_dir):
    print("restarting pipeline in existing directory " + args.out_dir)
else:
    subprocess.check_call(["mkdir", "-p", args.out_dir])
bam = args.bam
if not os.path.exists(args.out_dir + "/variants.done"):
    final_vcf = freebayes(args, bam, fasta)
else:
    with open(args.out_dir + "/variants.done") as done:
        final_vcf = done.readline().strip()
if not os.path.exists(args.out_dir + "/vartrix.done"):
    vartrix(args, final_vcf, bam)
print(final_vcf)
ref_mtx = args.out_dir + "/ref.mtx"
alt_mtx = args.out_dir + "/alt.mtx"
subprocess.check_call(["cp", args.barcodes, args.out_dir+"/."])

cellector_bin = args.cellector_binary if args.cellector_binary[0] == "/" else "./"+args.cellector_binary
cellector_cmd = [cellector_bin, "-a", alt_mtx, "-r", ref_mtx, "--output_directory",args.out_dir, 
                "--min_alt",args.min_alt, "--min_ref", args.min_ref, "--barcodes", args.barcodes, "--vcf", final_vcf]

print("running cellector")
print(" ".join(cellector_cmd))
with open(args.out_dir+"/cellector.err",'w') as err:
    with open(args.out_dir+"/cellector.out",'w') as out:
        subprocess.check_call(cellector_cmd,stdout=out, stderr=err)

souporcell_bin = args.souporcell_binary if args.souporcell_binary[0] == "/" else "./"+args.souporcell_binary
souporcell_cmd = [souporcell_bin,"-a", alt_mtx, "-r", ref_mtx, "--barcodes", args.barcodes,
    "-t", str(args.threads), "-k", "2", "--min_ref", str(args.min_ref), "--min_alt", str(args.min_alt)]
print("running souporcell")
print(" ".join(souporcell_cmd))
with open(args.out_dir+"/souporcell.err",'w') as err:
    with open(args.out_dir+"/souporcell.out", 'w') as out:
        subprocess.check_call(souporcell_cmd, stdout=out, stderr=err)
troublet_bin = args.troublet_binary if args.troublet_binary[0] == "/" else "./"+args.troublet_binary
troublet_cmd = [troublet_bin, "--alts", alt_mtx, "--refs", ref_mtx, "--clusters", args.out_dir+"/souporcell.out"]
with open(args.out_dir+"/troublet.err",'w') as err:
    with open(args.out_dir+"/troublet.out", 'w') as out:
        subprocess.check_call(troublet_cmd, stdout=out, stderr=err)

with open(args.out_dir+"/cellector_assignments.tsv", "r") as cellector_assign_fid:
    cellector_log_likelihood_0 = []
    cellector_log_likelihood_1 = []
    cellector_assign_fid.readline()
    for line in cellector_assign_fid.readlines():
        tokens = line.split("\t")
        majority_likelihood = float(tokens[6])
        minority_likelihood = float(tokens[7])
        val = majority_likelihood / np.average([majority_likelihood, minority_likelihood])
        if tokens[1] == "0":
            cellector_log_likelihood_0.append(val)
        elif tokens[1] == "1":
            cellector_log_likelihood_1.append(val)
    cellector_value = abs(np.average(cellector_log_likelihood_0) - np.average(cellector_log_likelihood_1))
with open(args.out_dir+"/troublet.out", "r") as troublet_assign_fid:
    troublet_log_likelihood_0 = []
    troublet_log_likelihood_1 = []
    for line in troublet_assign_fid.readlines():
        tokens = line.split("\t")
        if tokens[1] == "singlet":
            cluster_0_likelihood = float(tokens[7])
            cluster_1_likelihood = float(tokens[8])
            val = cluster_0_likelihood / np.average([cluster_0_likelihood, cluster_1_likelihood])
            if tokens[2] == "0":
                troublet_log_likelihood_0.append(val)
            elif tokens[2] == "1":
                troublet_log_likelihood_1.append(val)
    souporcell_value = abs(np.average(troublet_log_likelihood_0) - np.average(troublet_log_likelihood_1))

preference = "cellector" if cellector_value > souporcell_value else "souporcell"
print("prefering the output of ", preference)
cellector_values  = read_tsv(args.out_dir+"/cellector_assignments.tsv")
souporcell_values = read_tsv(args.out_dir+"/troublet.out")
final_output_values = [[] for row in cellector_values]
final_output_values[0].append("barcode")
final_output_values[0].append("assignment")
if preference == "cellector":
    for i in range(1, len(cellector_values)):
        row = cellector_values[i]
        final_output_values[i].append(row[0])  # barcode
        final_output_values[i].append(row[1])  # assignment
else:
    for i in range(1, len(souporcell_values)):
        row = souporcell_values[i]
        final_output_values[i].append(row[0])  # barcode
        status = row[1]
        final_output_values[i].append(row[2] if status == "singlet" else status)  # assignment
# the rest of the cellector columns    2,3...
header = cellector_values[0]
for j in range(2, len(header)):
    final_output_values[0].append("cellector_"+header[j])

for i in range(1, len(cellector_values)):
    row = cellector_values[i]
    for j in range(2, len(row)):
        final_output_values[i].append(row[j])
# the rest of the souporcell columns   1, 3,4...
header = souporcell_values[0]
for j in range(2, len(header)):
    final_output_values[0].append("souporcell_"+header[j])
for i in range(1, len(souporcell_values)):
    row = souporcell_values[i]
    final_output_values[i].append(row[1])
    for j in range(3, len(row)):
        final_output_values[i].append(row[j])

with open(args.out_dir+"final_output.out", 'w') as final_output_fid:
    for row in range(len(final_output_values)):
        final_output_fid.write("\t".join(final_output_values[row])+"\n")

grapher_cmd = ["python", args.grapher_script, "-d", args.out_dir]
print("running grapher")
print(" ".join(grapher_cmd))
with open(args.out_dir+"/grapher.err",'w') as err:
    with open(args.out_dir+"/grapher.out", 'w') as out:
        subprocess.check_call(grapher_cmd, stdout=out, stderr=err)


