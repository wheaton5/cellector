import argparse
import math
import random
import os
import subprocess

parser = argparse.ArgumentParser(
    description="takes the output of combiner and makes a given percentage of cells into doublets")
parser.add_argument("-a", "--alt", required = True)
parser.add_argument("-r", "--ref", required = True)
parser.add_argument("-g", "--ground_truth", required = True)
parser.add_argument("--percent_doublets", required = True,
    help = "percentage of doublets present in the output")
parser.add_argument("-o", "--out_dir", required = True,
    help = "name of directory to place souporcell files")
args = parser.parse_args()

# read cells into dictionary formatted this way: {cell_id:{locus_id: count, locus_id: count, ...}, cell_id:{...}, ...}
def read_cells_mtx(file_name):
    cell_dict = dict()
    header_lines = list()
    with open(file_name, "r") as fid:
        header_lines.append(fid.readline())
        header_lines.append(fid.readline())
        header_lines.append(fid.readline())
        for line in fid:
            (locus_id, cell_id, count) = line.split("\t")
            cell_id = int(cell_id)
            locus_id = int(locus_id)
            count = int(count)
            if cell_id not in cell_dict:
                cell_dict[cell_id] = dict()
            cell_dict[cell_id][locus_id] = count
    return (cell_dict, header_lines)

# pick out cells to become doublet. each doublet is a pair of cell_idss
def pick_doublet_ids(cell_ids, percent_doublets):
    fraction_to_pick = 0 if (float(percent_doublets) == 0.0) else (1 / (1 / (float(percent_doublets) / 100) + 1))
    num_doublets = math.floor(fraction_to_pick * len(cell_ids))
    doublet_ids = random.sample(cell_ids, num_doublets * 2)
    doublet_ids_set1 = doublet_ids[:len(doublet_ids)//2]  # first half
    doublet_ids_set2 = doublet_ids[len(doublet_ids)//2:]  # last half
    return list(zip(doublet_ids_set1, doublet_ids_set2))

# make doublets
# pairs of cells in cell_dict with the keys of doublet_ids will have their alt/ref values summed
def make_doublets(cell_dict, doublet_ids):
    for pair_idx in range(0, len(doublet_ids)):
        cell_id1, cell_id2 = doublet_ids[pair_idx]
        cell1 = cell_dict[cell_id1]
        cell2 = cell_dict[cell_id2]
        for locus_id in list(cell2):
            if locus_id in cell1:
                cell1[locus_id] = cell1[locus_id] + cell2[locus_id]
            else:
                cell1[locus_id] = cell2[locus_id]
        del cell_dict[cell_id2]
    return cell_dict

def write_gt_and_bc(original_gt_filename, doublet_ids, out_bc_filename, out_gt_filename):
    removed_doublet_ids = [pair[1] for pair in doublet_ids]
    kept_doublet_ids = [pair[0] for pair in doublet_ids]
    barcodes = []
    ground_truths = []
    with open(original_gt_filename, "r") as original_gt_fid:
        for line in original_gt_fid:
            bc, gt = line.split("\t")
            barcodes.append(bc)
            ground_truths.append(gt)
    with open(out_gt_filename, "w") as out_gt_fid, open(out_bc_filename, "w") as out_bc_fid:
        for cell_id in range(1, len(barcodes)+1):
            if cell_id in removed_doublet_ids:
                continue
            bc = barcodes[cell_id - 1]
            gt = ground_truths[cell_id - 1]
            if cell_id in kept_doublet_ids:
                paired_doublet_id = removed_doublet_ids[kept_doublet_ids.index(cell_id)]
                paired_gt = ground_truths[paired_doublet_id - 1]
                if gt != paired_gt:
                    gt = "doublet\n"
                    #bc =  # currently unchanged
            out_gt_fid.write(f"{bc}\t{gt}")
            out_bc_fid.write(f"{bc}\n")

# cell id mapping  # cell ids match the line number of a barcode in the barcodes file
def remap_cell_ids(cell_dict):
    counter = 1
    for cell_id in sorted(cell_dict):
        cell_dict[counter] = cell_dict.pop(cell_id)
        counter += 1
    return cell_dict

# output new mtx files
def write_cells_to_mtx(cell_dict, header_lines, file_name):
    # organize data better for writing
    lines_info = [(locus_id, cell_id, cell_dict[cell_id][locus_id]) for cell_id in cell_dict for locus_id in cell_dict[cell_id]]  # unreadable
    # sort information by locus_id then cell_id
    lines_info = sorted(lines_info)
    with open(file_name, "w") as fid:
        fid.write(header_lines[0])
        fid.write(header_lines[1])
        num_lines = len(lines_info)
        num_cells = len(cell_dict)
        fid.write(f"{num_lines}\t{num_cells}\t0\n")  # it is unclear what the third number should be
        for line_info in lines_info:
            locus_id = line_info[0]
            cell_id = line_info[1]
            count = line_info[2]
            fid.write(f"{locus_id}\t{cell_id}\t{count}\n")


# main
alt_cell_dict, alt_header_lines = read_cells_mtx(args.alt)
ref_cell_dict, ref_header_lines = read_cells_mtx(args.ref)
doublet_ids = pick_doublet_ids(list(alt_cell_dict), args.percent_doublets)
make_doublets(alt_cell_dict, doublet_ids)
make_doublets(ref_cell_dict, doublet_ids)
remap_cell_ids(alt_cell_dict)
remap_cell_ids(ref_cell_dict)
if not os.path.isdir(args.out_dir):
    subprocess.check_call(["mkdir", "-p", args.out_dir])
write_cells_to_mtx(alt_cell_dict, alt_header_lines, f"{args.out_dir}/alt.mtx")
write_cells_to_mtx(ref_cell_dict, ref_header_lines, f"{args.out_dir}/ref.mtx")
write_gt_and_bc(args.ground_truth, doublet_ids, f"{args.out_dir}/bc.tsv", f"{args.out_dir}/gt.tsv")

