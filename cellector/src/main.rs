#[macro_use]
extern crate clap;
extern crate hashbrown;
extern crate statrs;
extern crate flate2;
extern crate itertools;

mod stats;
mod load_data;
use load_data::CellData;
use load_data::VcfLocusData;

use clap::App;
use std::fs::File;
use std::io::{BufWriter, Write};

use hashbrown::{HashMap,HashSet};
use itertools::izip;

use statrs::statistics::OrderStatistics;
use statrs::statistics::Data;

fn main() {
    let params = load_params();
    load_data::create_output_dir(&params);
    let (cell_id_to_barcode, barcode_to_cell_id) = load_data::load_barcodes(&params);
    let cell_id_to_assignment = load_data::load_ground_truth(&params, &barcode_to_cell_id);
    let (mut loci_used, locus_ids, cell_data, locus_counts, precomputed_log_binomial_coefficients) = 
        load_data::load_cell_data(&params, &cell_id_to_barcode, &cell_id_to_assignment);
    let vcf_data = load_data::load_vcf_data(&params);
    cellector(&params, &mut loci_used, &locus_ids, &cell_data, &locus_counts, &vcf_data, &precomputed_log_binomial_coefficients);
}

fn cellector(params: &Params, loci_used: &mut Vec<bool>, locus_ids: &Vec<usize>, cell_data: &Vec<CellData>, locus_counts: &Vec<[f64; 2]>, vcf_data: &Option<Vec<VcfLocusData>>, precomputed_log_binomial_coefficients: &Vec<Vec<f64>>) {
    let mut excluded_cells: HashSet<usize> = HashSet::new();
    let mut any_change;
    let mut iteration = 0;
    loop {
        (any_change, excluded_cells) = compute_new_excluded(params, loci_used, locus_ids, cell_data, locus_counts, &excluded_cells, iteration, vcf_data, precomputed_log_binomial_coefficients);
        iteration += 1;
        if !any_change { break; }
    }
    let posteriors = calculate_posteriors(params, loci_used, cell_data, locus_counts, &excluded_cells, precomputed_log_binomial_coefficients);
    output_final_assignments(params, cell_data, &posteriors, &excluded_cells);
}

fn output_final_assignments(params: &Params, cell_data: &Vec<CellData>, posteriors: &Vec<f64>, excluded_cells: &HashSet<usize>) {
    let filename = format!("{}/cellector_assignments.tsv",params.output_directory);
    let filehandle = File::create(&filename).expect(&format!("Unable to create file {}", &filename));
    let mut writer = BufWriter::new(filehandle);
    let header = format!("barcode\tposterior_assignment\tanomally_assignment\tminority_posterior\tmajority_posterior\tground_truth_assignment\n");
    writer.write_all(header.as_bytes()).expect("could not write to cellector assignment file");
    let mut assignment_gt_counts: HashMap<String, HashMap<String, usize>> = HashMap::new();
    let mut gt_counts: HashMap<String, usize> = HashMap::new();
    for cell_id in 0..cell_data.len() {
        let cell = &cell_data[cell_id];
        let mut posterior_assignment = "unassigned";
        if posteriors[cell_id] > params.posterior_threshold {
            posterior_assignment = "0";
        } else if 1.0 - posteriors[cell_id] > params.posterior_threshold {
            posterior_assignment = "1";
        }
        let counts = assignment_gt_counts.entry(posterior_assignment.to_string()).or_insert(HashMap::new());
        let count = counts.entry(cell.assignment.clone()).or_insert(0);
        *count += 1;
        let count = gt_counts.entry(cell.assignment.clone()).or_insert(0);
        *count += 1;
        
        let anomally_assignment;
        if excluded_cells.contains(&cell_id) {
            anomally_assignment = "0";
        } else { anomally_assignment = "1"; } // maybe do something here where we leave some unassigned if they are near the threshold?
        let line = format!("{}\t{}\t{}\t{:.5}\t{:.5}\t{}\n", cell.barcode, posterior_assignment, anomally_assignment, posteriors[cell_id], 1.0 - posteriors[cell_id], cell.assignment);
        writer.write_all(line.as_bytes()).expect("could not write to cellector assignment file");
    }
    pretty_print(params, assignment_gt_counts, gt_counts);
}

//TODO add detailed output for incorrectly assigned cells
fn pretty_print(params: &Params, assignment_gt_counts: HashMap<String, HashMap<String, usize>>, gt_counts: HashMap<String, usize>) {
    let mut count_vec: Vec<(&String, &usize)> = gt_counts.iter().collect();
    count_vec.sort_by(|a, b| b.1.cmp(a.1));    
    let mut string_build: String = String::new();
    let first_header = "cellector assignment   ".to_string();
    let header = "      0      1      unassigned\n".to_string();
    let header_offset = first_header.len();
    string_build.push_str(&first_header);
    string_build.push_str(&header);
    let mut xoffset = 3;
    xoffset = xoffset.max(first_header.len()+2);
    string_build.push_str("cell_hashing");
    for _i in 0..(xoffset.checked_sub(12).unwrap_or(0)) { string_build.push_str(" "); }
    string_build.push_str("|");
    for _i in 0..(header.len().checked_sub(1).unwrap_or(0)) { string_build.push_str("-"); }
    string_build.push_str("|\n");
    for (gt, count) in & count_vec {
        xoffset = xoffset.max(gt.len() + 3);
        let (mut count0, mut count1, mut unassigned) = (0,0,0);
        if assignment_gt_counts.contains_key("0") {
            let counts = assignment_gt_counts.get("0").unwrap();
            count0 = *counts.get(&gt.to_string()).unwrap_or(&0);
        }
        if assignment_gt_counts.contains_key("1") {
            let counts = assignment_gt_counts.get("1").unwrap();
            count1 = *counts.get(&gt.to_string()).unwrap_or(&0);
        }
        if assignment_gt_counts.contains_key("unassigned") {
            let counts = assignment_gt_counts.get("unassigned").unwrap();
            unassigned = *counts.get(&gt.to_string()).unwrap_or(&0);
        }
        let count0 = format!("{}",count0);
        let count1 = format!("{}",count1);
        let unassigned = format!("{}",unassigned);
        string_build.push_str(gt);
        for _i in 0..(xoffset.checked_sub(gt.len().checked_sub(1).unwrap_or(0)).unwrap_or(0)) { string_build.push_str(" "); }
        string_build.push_str(&format!(" |  {}",count0));
        for _i in 0..(4_usize.checked_sub(count0.len()).unwrap_or(0)) { string_build.push_str(" "); }
        string_build.push_str(&format!(" |  {}",count1));
        for _i in 0..(4_usize.checked_sub(count1.len()).unwrap_or(0)) { string_build.push_str(" "); }
        string_build.push_str(&format!(" |  {}",unassigned));
        for _i in 0..(12_usize.checked_sub(unassigned.len()).unwrap_or(0)) { string_build.push_str(" "); }
        string_build.push_str("|\n");
    }
    for _i in 0..xoffset { string_build.push_str(" "); }
    string_build.push_str("|");
    for _i in 0..(header.len().checked_sub(1).unwrap_or(0) ) { string_build.push_str("-"); }
    string_build.push_str("|\n");
    println!("\n\n{}",string_build);
} 

fn calculate_posteriors(params: &Params, loci_used: &Vec<bool>, cell_data: &Vec<CellData>, locus_counts: &Vec<[f64;2]>, excluded_cells: &HashSet<usize>, precomputed_log_binomial_coefficients: &Vec<Vec<f64>>) -> Vec<f64> {
    let mut posteriors: Vec<f64> = Vec::new();
    let mut included_cells: HashSet<usize> = HashSet::new();
    for cell_id in 0..cell_data.len() {
        if !excluded_cells.contains(&cell_id) {
            included_cells.insert(cell_id);
        }
    }
    let mut alpha_betas_majority_dist = init_alpha_betas(locus_counts, excluded_cells, cell_data);
    let minority_fraction = (excluded_cells.len() as f64)/(cell_data.len() as f64);
    let alpha_betas_minority_dist = init_alpha_betas(locus_counts, &included_cells, cell_data);
    for locus in 0..loci_used.len() {
        alpha_betas_majority_dist[locus].alpha = (alpha_betas_majority_dist[locus].alpha - 1.0)*minority_fraction + 1.0;
        alpha_betas_majority_dist[locus].beta = (alpha_betas_majority_dist[locus].beta - 1.0)*minority_fraction + 1.0;
    }
    let loci_used_for_posteriors: Vec<bool> = get_loci_used_for_posterior_calc(params, loci_used, cell_data, excluded_cells, locus_counts);
    let minority_dist_likelihoods: CellLogLikelihoodData = get_cell_log_likelihoods(&loci_used_for_posteriors, cell_data, &alpha_betas_minority_dist, excluded_cells, precomputed_log_binomial_coefficients);
    let majority_dist_likelihoods: CellLogLikelihoodData = get_cell_log_likelihoods(&loci_used_for_posteriors, cell_data, &alpha_betas_majority_dist, &included_cells, precomputed_log_binomial_coefficients);
    let log_prior_minority: f64 = minority_fraction.ln();
    let log_prior_majority: f64 = (1.0 - minority_fraction).ln();
    for cell_id in 0..cell_data.len() {
        let log_numerator = log_prior_minority + minority_dist_likelihoods.log_likelihoods[cell_id];
        let log_denominator = stats::logsumexp(log_numerator, log_prior_majority + majority_dist_likelihoods.log_likelihoods[cell_id]);
        let log_minority_posterior = log_numerator - log_denominator;
        let posterior = log_minority_posterior.exp();
        posteriors.push(posterior);
    }
    return posteriors;
}

fn get_loci_used_for_posterior_calc(params: &Params, loci_used: &Vec<bool>, cell_data: &Vec<CellData>, excluded_cells: &HashSet<usize>, locus_counts: &Vec<[f64; 2]>) -> Vec<bool> {
    let mut locus_counts_minority: Vec<[usize; 2]> = Vec::new();
    let mut loci_used_for_posteriors: Vec<bool> = Vec::new();
    for _i in 0..loci_used.len() {
        locus_counts_minority.push([0;2]);
    }
    for cell_id in excluded_cells {
        let cell = &cell_data[*cell_id];
        for locus in &cell.cell_loci_data {
            locus_counts_minority[locus.locus_index][0] += locus.ref_count.round() as usize;
            locus_counts_minority[locus.locus_index][1] += locus.alt_count.round() as usize;
        }
    }
    for locus_index in 0..loci_used.len() {
        let minority_alt = locus_counts_minority[locus_index][1];
        let minority_ref = locus_counts_minority[locus_index][0];
        let majority_alt = locus_counts[locus_index][1].round() as usize - locus_counts_minority[locus_index][1];
        let majority_ref = locus_counts[locus_index][0].round() as usize - locus_counts_minority[locus_index][0];
        let min_alleles = params.min_alleles_posterior;
        if loci_used[locus_index] && minority_alt + minority_ref >= min_alleles && majority_alt + majority_ref >= min_alleles {
            loci_used_for_posteriors.push(true);
        } else { loci_used_for_posteriors.push(false); }
    }
    return loci_used_for_posteriors;
}

fn compute_new_excluded(params: &Params, loci_used: &mut Vec<bool>, locus_ids: &Vec<usize>, cell_data: &Vec<CellData>, locus_counts: &Vec<[f64; 2]>, excluded_cells: &HashSet<usize>, iteration: usize, vcf_data: &Option<Vec<VcfLocusData>>, precomputed_log_binomial_coefficients: &Vec<Vec<f64>>) -> (bool, HashSet<usize>) {
    let alpha_betas = init_alpha_betas(locus_counts, excluded_cells, cell_data);
    let mut new_excluded: HashSet<usize> = HashSet::new();

    let cell_log_likelihood_data = get_cell_log_likelihoods(loci_used, cell_data, &alpha_betas, excluded_cells, precomputed_log_binomial_coefficients);
    let mut normalized_log_likelihoods: Vec<f64> = Vec::new();
    for i in 0..cell_data.len() {
        if cell_log_likelihood_data.loci_used_per_cell[i] > 0.0 {
            normalized_log_likelihoods.push(cell_log_likelihood_data.log_likelihoods[i] / cell_log_likelihood_data.loci_used_per_cell[i]);
            //normalized_log_likelihoods.push((cell_log_likelihood_data.log_likelihoods[i] - cell_log_likelihood_data.expected_log_likelihoods[i]) /
            //    cell_log_likelihood_data.expected_log_variances[i].sqrt());
        
        } else { // need to filter cells with < some number of used loci TODO
            normalized_log_likelihoods.push(f64::NEG_INFINITY);
        }
    }
    let mut normalized_tmp = Data::new(normalized_log_likelihoods.clone());
    let median = normalized_tmp.median();
    let q1 = normalized_tmp.lower_quartile();
    let q3 = normalized_tmp.upper_quartile();
    let iqr = q3 - q1;
    let threshold = q1 - params.interquartile_range_multiple * iqr;
    for (cell_id, normalized_likelihood) in normalized_log_likelihoods.iter().enumerate() {
        if *normalized_likelihood < threshold { new_excluded.insert(cell_id); }
    }
    let num_new_cells_excluded = new_excluded.difference(&excluded_cells).collect::<Vec<&usize>>().len();
    let num_cells_rescued = excluded_cells.difference(&new_excluded).collect::<Vec<&usize>>().len();
    let any_change = num_new_cells_excluded > 0 || num_cells_rescued > 0;
    
    let locus_data = get_locus_log_likelihoods(&cell_log_likelihood_data.all_pmfs, cell_data, loci_used, &new_excluded);
    locus_filter_and_output_locus_data(params, loci_used, &locus_data, locus_ids, vcf_data, iteration);
    println!("detected {} new anomylous cells and rescued {} cells to the majority in iteration {}", num_new_cells_excluded, num_cells_rescued, iteration+1);
    println!("median normalized log likelihood {} with interquartile range {}, threshold {}", median, iqr, threshold);
    output_iteration_tsv(params, cell_data, &cell_log_likelihood_data, threshold, iteration);
    return (any_change, new_excluded);
}

fn output_iteration_tsv(params: &Params, cell_data: &Vec<CellData>, cell_log_likelihood_data: &CellLogLikelihoodData, threshold: f64, iteration: usize) {
    let filename = format!("{}/iteration_{}.tsv",params.output_directory, iteration);
    let filehandle = File::create(&filename).expect(&format!("Unable to create file {}", &filename));
    let mut writer = BufWriter::new(filehandle);
    let header = format!("cell_id\tbarcode\tassignment\tlog_likelihood\texpected_log_likelihood\tnum_loci_used\n");
    writer.write_all(header.as_bytes()).expect("could not write to iteration tsv file");
    for cell_id in 0..cell_data.len() {
        let cell = &cell_data[cell_id];
        assert!(cell.cell_id == cell_id, "I did something wrong, cell_id != cell_data[cell_id].cell_id");
        let line = format!("{}\t{}\t{}\t{}\t{}\t{}\n", cell.cell_id, cell.barcode, cell.assignment, cell_log_likelihood_data.log_likelihoods[cell_id], cell_log_likelihood_data.expected_log_likelihoods[cell_id], cell_log_likelihood_data.loci_used_per_cell[cell_id]);
        writer.write_all(line.as_bytes()).expect("could not write to iteration tsv file"); 
    }
    let filename = format!("{}/iteration_{}_threshold.tsv",params.output_directory, iteration);
    let filehandle = File::create(&filename).expect(&format!("Unable to create file {}", &filename));
    let mut writer = BufWriter::new(filehandle);
    let line = format!("{}",threshold);
    writer.write_all(line.as_bytes()).expect("could not write to threshold file");
}

fn get_locus_log_likelihoods(all_pmfs: &Vec<PMFData>, cell_data: &Vec<CellData>, loci_used: &Vec<bool>, excluded_cells: &HashSet<usize>) -> LocusLogLikelihoodData {
    let mut locus_contributions_minority: Vec<f64> = Vec::new();
    let mut locus_contributions_majority: Vec<f64> = Vec::new();
    let mut locus_cells_minority: Vec<usize> = Vec::new();
    let mut locus_cells_majority: Vec<usize> = Vec::new();
    let mut locus_expected_contribution_minority: Vec<f64> = Vec::new();
    let mut locus_expected_contribution_majority: Vec<f64> = Vec::new();
    let mut locus_alleles_minority: Vec<AlleleCount> = Vec::new();
    let mut locus_alleles_majority: Vec<AlleleCount> = Vec::new();
    
    for _locus in 0..loci_used.len() {
        locus_contributions_minority.push(0.0);
        locus_contributions_majority.push(0.0);
        locus_cells_minority.push(0);
        locus_cells_majority.push(0);
        locus_expected_contribution_minority.push(0.0);
        locus_expected_contribution_majority.push(0.0);
        locus_alleles_minority.push(AlleleCount{ ref_count: 0, alt_count: 0});
        locus_alleles_majority.push(AlleleCount{ ref_count: 0, alt_count: 0});
    }

    let mut cells_to_excluded: Vec<bool> = Vec::new();
    for _ in 0..cell_data.len() { cells_to_excluded.push(false); }
    for excluded in excluded_cells { cells_to_excluded[*excluded] = true; }
    for pmf_data in all_pmfs {
        let log_pmf = pmf_data.log_pmf;
        let expected_log_pmf = pmf_data.log_pmf;
                if cells_to_excluded[pmf_data.cell_id] {
            locus_contributions_minority[pmf_data.locus_index] += log_pmf;
            locus_cells_minority[pmf_data.locus_index] += 1;
            locus_expected_contribution_minority[pmf_data.locus_index] += expected_log_pmf;
            locus_alleles_minority[pmf_data.locus_index].ref_count += pmf_data.ref_count as usize;
            locus_alleles_minority[pmf_data.locus_index].alt_count += pmf_data.alt_count as usize;
        } else {
            locus_contributions_majority[pmf_data.locus_index] += log_pmf;
            locus_cells_majority[pmf_data.locus_index] += 1;
            locus_expected_contribution_majority[pmf_data.locus_index] += expected_log_pmf;
            locus_alleles_majority[pmf_data.locus_index].ref_count += pmf_data.ref_count as usize;
            locus_alleles_majority[pmf_data.locus_index].alt_count += pmf_data.alt_count as usize;
        }
    }
                
    return LocusLogLikelihoodData{
        locus_contributions_minority: locus_contributions_minority,
        locus_contributions_majority: locus_contributions_majority,
        locus_cells_minority: locus_cells_minority,
        locus_cells_majority: locus_cells_majority,
        locus_alleles_minority: locus_alleles_minority,
        locus_alleles_majority: locus_alleles_majority,
        locus_expected_contribution_minority: locus_expected_contribution_minority,
        locus_expected_contribution_majority: locus_expected_contribution_majority,
    };
}

fn locus_filter_and_output_locus_data(params: &Params, used_loci: &mut Vec<bool>, likelihood_data: &LocusLogLikelihoodData, locus_ids: &Vec<usize>, vcf_data: &Option<Vec<VcfLocusData>>, iteration: usize) {
    let filename = format!("{}/iteration_{}_locus_contribution.tsv",params.output_directory, iteration);
    let filehandle = File::create(&filename).expect(&format!("Unable to create file {}", &filename));
    let mut writer = BufWriter::new(filehandle);
    writer.write_all(format!("locus_id\tchrom\tpos\tlog_likelihood_minority\tlog_likelihood_majority\texpected_loglike_minority\texpected_loglike_majority\tminority_cellcount\tmajority_cellcount\tlog_likelihood_minority_per_cell\tlog_likelihood_majority_per_cell\tminority_alt\tminority_ref\tmajority_alt\tmajority_ref\tminority_af\tmajority_af\n").as_bytes()).expect("could not write locus contribution file");

    let mut locus_loglike_per_cell_minority: Vec<f64> = Vec::new();
    let mut locus_loglike_per_cell_majority: Vec<f64> = Vec::new();
    let mut locus_loglike_per_cell_minority_for_threshold: Vec<f64> = Vec::new();
    for i in 0..likelihood_data.locus_contributions_minority.len() {
        if likelihood_data.locus_cells_minority[i] != 0 {
            locus_loglike_per_cell_minority.push(likelihood_data.locus_contributions_minority[i] / (likelihood_data.locus_cells_minority[i] as f64));
            locus_loglike_per_cell_minority_for_threshold.push(likelihood_data.locus_contributions_minority[i] / (likelihood_data.locus_cells_minority[i] as f64));
        } else { locus_loglike_per_cell_minority.push(0.0); }
        if likelihood_data.locus_cells_majority[i] != 0 {
            locus_loglike_per_cell_majority.push(likelihood_data.locus_contributions_majority[i] / (likelihood_data.locus_cells_majority[i] as f64));
        } else { locus_loglike_per_cell_majority.push(0.0); }
    }
    let ordering = argsort(&locus_loglike_per_cell_minority);
    // get median and filter things 100*median
    let mut data = Data::new(locus_loglike_per_cell_minority_for_threshold);
    let median = data.median();
    let threshold = median * 100.0; // TODO don't code things or choose arbitrary numbers
    for locus_index in 0..used_loci.len() {
        if locus_loglike_per_cell_minority[locus_index] < threshold {
            used_loci[locus_index] = false;
            println!("filtering locus {} locus index {} because it was contributing {} vs median {} per cell to log likelihood of minority cells", 
                locus_ids[locus_index], locus_index, locus_loglike_per_cell_minority[locus_index], median);
        }
    }
    
    for index in ordering {
        let minority_loglike = likelihood_data.locus_contributions_minority[index];
        let majority_loglike = likelihood_data.locus_contributions_majority[index];
        let minority_loglike_percell = locus_loglike_per_cell_minority[index];
        let majority_loglike_percell = locus_loglike_per_cell_majority[index];
        let locus_id = locus_ids[index];
        let minority_alt = likelihood_data.locus_alleles_minority[index].alt_count;
        let minority_ref = likelihood_data.locus_alleles_minority[index].ref_count;
        let majority_alt = likelihood_data.locus_alleles_majority[index].alt_count;    
        let majority_ref = likelihood_data.locus_alleles_majority[index].ref_count;    
        let minority_expected_loglike = likelihood_data.locus_expected_contribution_minority[index];
        let majority_expected_loglike = likelihood_data.locus_expected_contribution_majority[index];
        let minority_cells = likelihood_data.locus_cells_minority[index];
        let majority_cells = likelihood_data.locus_cells_majority[index];
        let mut minority_af = 0.0;
        let mut majority_af = 0.0;
        if minority_alt + minority_ref > 0 { minority_af = (minority_alt as f64)/((minority_alt + minority_ref) as f64); }
        if majority_alt + majority_ref > 0 { majority_af = (majority_alt as f64)/((majority_alt + majority_ref) as f64); }
        let mut chrom = "na";
        let mut pos = "na";
        if let Some(vcf_data) = vcf_data {
            chrom = &vcf_data[locus_id].chrom;
            pos = &vcf_data[locus_id].pos;
        }
        let line = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                locus_id,
                chrom,
                pos,
                minority_loglike,
                majority_loglike,
                minority_expected_loglike,
                majority_expected_loglike,
                minority_cells,
                majority_cells,
                minority_loglike_percell,
                majority_loglike_percell,
                minority_alt,
                minority_ref,
                majority_alt,
                majority_ref,
                minority_af,
                majority_af
                );
        writer.write_all(line.as_bytes()).expect("could not write to locus contributions file");
    } 
}

struct AlleleCount { ref_count: usize, alt_count: usize, }

pub fn argsort<T: PartialOrd>(data: &Vec<T>) -> Vec<usize> {
    let mut indices = (0..data.len()).collect::<Vec<_>>();
    indices.sort_by(|&i, &j| data[i].partial_cmp(&data[j]).unwrap());
    return indices;
}
    
struct CellLogLikelihoodData {
    log_likelihoods: Vec<f64>,
    loci_used_per_cell: Vec<f64>,
    expected_log_likelihoods: Vec<f64>,
    expected_log_variances: Vec<f64>,
    all_pmfs: Vec<PMFData>, // for data analysis, probably will remove later
}

struct LocusLogLikelihoodData {
    locus_contributions_minority: Vec<f64>,
    locus_contributions_majority: Vec<f64>,
    locus_cells_minority: Vec<usize>,
    locus_cells_majority: Vec<usize>,
    locus_alleles_minority: Vec<AlleleCount>,
    locus_alleles_majority: Vec<AlleleCount>,
    locus_expected_contribution_minority: Vec<f64>,
    locus_expected_contribution_majority: Vec<f64>,
}

struct PMFData { // this is overkill, but maybe will provide some insight
    cell_id: usize,
    locus_index: usize,
    locus: usize,
    excluded: bool,
    log_pmf: f64,
    alt_count: usize,
    ref_count: usize,
    alpha: f64,
    beta: f64,
    expected_log_pmf: f64,
    expected_log_variance: f64,
}

fn get_cell_log_likelihoods(loci_used: &Vec<bool>, cell_data: &Vec<CellData>, alpha_betas: &Vec<AlphaBeta>, excluded_cells: &HashSet<usize>, precomputed_log_binomial_coefficients: &Vec<Vec<f64>>) -> 
        CellLogLikelihoodData {
    let mut cell_log_likelihoods: Vec<f64> = Vec::new();
    let mut cell_loci_used: Vec<f64> = Vec::new();
    let mut cell_expected_log_likelihood: Vec<f64> = Vec::new();
    let mut cell_expected_log_variance: Vec<f64> = Vec::new();
    let mut all_pmfs: Vec<PMFData> = Vec::new();
    for cell_id in 0..cell_data.len() {
        let mut log_likelihood: f64 = 0.0;
        let mut expected_log_likelihood: f64 = 0.0;
        let mut expected_log_variance: f64 = 0.0;
        let mut loci_used_for_cell: f64 = 0.0;
        let excluded: bool = excluded_cells.contains(&cell_id);
        let cell = &cell_data[cell_id];
        for locus in &cell.cell_loci_data {
            if loci_used[locus.locus_index] {
                let log_pmf = stats::log_beta_binomial_pmf(locus.alt_count, locus.ref_count, alpha_betas[locus.locus_index].alpha, alpha_betas[locus.locus_index].beta, locus.log_binomial_coefficient);
                log_likelihood += log_pmf;
                let expected_log_pmf = stats::expected_log_beta_binomial_pmf(locus.total, alpha_betas[locus.locus_index].alpha, alpha_betas[locus.locus_index].beta, precomputed_log_binomial_coefficients);
                expected_log_likelihood += expected_log_pmf.expected_log_likelihood;
                expected_log_variance += expected_log_pmf.expected_log_variance;
                all_pmfs.push(PMFData{
                    log_pmf: log_pmf,
                    cell_id: cell_id,
                    excluded: excluded,
                    locus: locus.locus_id,
                    locus_index: locus.locus_index,
                    alt_count: locus.alt_count as usize,
                    ref_count: locus.ref_count as usize,
                    alpha: alpha_betas[locus.locus_index].alpha,
                    beta: alpha_betas[locus.locus_index].beta,
                    expected_log_pmf: expected_log_pmf.expected_log_likelihood,
                    expected_log_variance: expected_log_pmf.expected_log_variance,
                });
                loci_used_for_cell += 1.0;
            }
        }
        cell_log_likelihoods.push(log_likelihood);
        cell_loci_used.push(loci_used_for_cell);
        cell_expected_log_likelihood.push(expected_log_likelihood);
        cell_expected_log_variance.push(expected_log_variance);
    }
    return CellLogLikelihoodData{
        log_likelihoods: cell_log_likelihoods, 
        loci_used_per_cell: cell_loci_used,
        expected_log_likelihoods: cell_expected_log_likelihood,
        expected_log_variances: cell_expected_log_variance,
        all_pmfs: all_pmfs,
    };
        
}

struct AlphaBeta {
    alpha: f64,
    beta: f64,
}

fn init_alpha_betas(locus_counts: &Vec<[f64;2]>, excluded_cells: &HashSet<usize>, cell_data: &Vec<CellData>) -> Vec<AlphaBeta> {
    let mut alpha_betas: Vec<AlphaBeta> = Vec::new();
    for counts in locus_counts {
        alpha_betas.push(AlphaBeta{alpha: counts[1] + 1.0, beta: counts[0] + 1.0});
    }
    for excluded_cell in excluded_cells {
        let cell = &cell_data[*excluded_cell];
        for locus_data in &cell.cell_loci_data {
            alpha_betas[locus_data.locus_index].alpha -= locus_data.alt_count as f64;
            alpha_betas[locus_data.locus_index].beta -= locus_data.ref_count as f64;
        }
    }
    return alpha_betas;
}

pub struct Params {
    ref_mtx: String,
    alt_mtx: String,
    barcodes: String,
    min_alt: usize,
    min_ref: usize,
    ground_truth: Option<String>,
    vcf: Option<String>,
    posterior_threshold: f64,
    interquartile_range_multiple: f64,
    output_directory: String,
    min_alleles_posterior: usize,
}

fn load_params() -> Params{
    let yaml = load_yaml!("params.yml");
    let params = App::from_yaml(yaml).get_matches();
    let ref_mtx = params.value_of("ref").unwrap().to_string();
    let alt_mtx = params.value_of("alt").unwrap().to_string();
    let barcodes = params.value_of("barcodes").unwrap().to_string();
    let min_alt = params.value_of("min_alt").unwrap_or("4");
    let min_alt = min_alt.to_string().parse::<usize>().unwrap();
    let min_ref = params.value_of("min_ref").unwrap_or("4");
    let min_ref = min_ref.to_string().parse::<usize>().unwrap();
    let ground_truth: Option<String> = match params.value_of("ground_truth") {
        None => None,
        Some(x) => Some(x.to_string()),
    };
    let vcf: Option<String> = match params.value_of("vcf") {
        None => None,
        Some(x) => Some(x.to_string()),
    };
    let posterior_threshold = params.value_of("posterior_threshold").unwrap_or("0.999");
    let posterior_threshold = posterior_threshold.to_string().parse::<f64>().unwrap();
    let interquartile_range_multiple = params.value_of("interquartile_range_multiple").unwrap_or("5");
    let interquartile_range_multiple = interquartile_range_multiple.to_string().parse::<f64>().unwrap();
    let output_directory = params.value_of("output_directory").unwrap().to_string();
    let min_alleles_posterior = params.value_of("min_alleles_posterior").unwrap_or("5");
    let min_alleles_posterior = min_alleles_posterior.to_string().parse::<usize>().unwrap();

    let params = Params {
        ref_mtx: ref_mtx,
        alt_mtx: alt_mtx,
        barcodes: barcodes,
        vcf: vcf,
        min_alt: min_alt,
        min_ref: min_ref,
        ground_truth: ground_truth, 
        posterior_threshold: posterior_threshold,
        interquartile_range_multiple: interquartile_range_multiple,
        output_directory: output_directory,
        min_alleles_posterior: min_alleles_posterior,
    };
    return params;
}


