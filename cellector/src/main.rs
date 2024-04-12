#[macro_use]
extern crate clap;
extern crate hashbrown;
extern crate statrs;
extern crate flate2;
extern crate itertools;
extern crate logaddexp;

mod stats;
mod load_data;
use load_data::CellData;

use logaddexp::LogAddExp;

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
    let (mut loci_used, cell_data, locus_counts) = load_data::load_cell_data(&params, &cell_id_to_barcode, &cell_id_to_assignment);
    
    cellector(&params, &mut loci_used, &cell_data, &locus_counts);
}

fn cellector(params: &Params, loci_used: &mut Vec<bool>, cell_data: &Vec<CellData>, locus_counts: &Vec<[f64; 2]>) {
    let mut excluded_cells: HashSet<usize> = HashSet::new();
    let mut any_change;
    let mut iteration = 0;
    loop {
        (any_change, excluded_cells) = compute_new_excluded(params, loci_used, cell_data, locus_counts, &excluded_cells, iteration);
        iteration += 1;
        if !any_change { break; }
    }
    let posteriors = calculate_posteriors(params, loci_used, cell_data, locus_counts, &excluded_cells);
    output_final_assignments(params, cell_data, &posteriors, &excluded_cells);
}

fn output_final_assignments(params: &Params, cell_data: &Vec<CellData>, posteriors: &Vec<f64>, excluded_cells: &HashSet<usize>) {
    let filename = format!("{}/cellector_assignments.tsv",params.output_directory);
    let filehandle = File::create(&filename).expect(&format!("Unable to create file {}", &filename));
    let mut writer = BufWriter::new(filehandle);
    let header = format!("barcode\tposterior_assignment\tanomally_assignment\tminority_posterior\tmajority_posterior\tground_truth_assignment\n");
    writer.write_all(header.as_bytes()).expect("could not write to cellector assignment file");
    let mut assignment_gt_counts: HashMap<String, HashMap<String, usize>> = HashMap::new();
    let mut gts: HashSet<String> = HashSet::new();
    let mut gt_counts: HashMap<String, usize> = HashMap::new();
    for cell_id in 0..cell_data.len() {
        let cell = &cell_data[cell_id];
        let mut posterior_assignment = "na";
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
        gts.insert(cell.assignment.clone());
        
        let anomally_assignment;
        if excluded_cells.contains(&cell_id) {
            anomally_assignment = "0";
        } else { anomally_assignment = "1"; } // maybe do something here where we leave some unassigned if they are near the threshold?
        let line = format!("{}\t{}\t{}\t{:.5}\t{:.5}\t{}\n", cell.barcode, posterior_assignment, anomally_assignment, posteriors[cell_id], 1.0 - posteriors[cell_id], cell.assignment);
        writer.write_all(line.as_bytes()).expect("could not write to cellector assignment file");
    }
}

fn calculate_posteriors(params: &Params, loci_used: &Vec<bool>, cell_data: &Vec<CellData>, locus_counts: &Vec<[f64;2]>, excluded_cells: &HashSet<usize>) -> Vec<f64> {
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
    let minority_dist_likelihoods: CellLogLikelihoodData = get_cell_log_likelihoods(&loci_used_for_posteriors, cell_data, &alpha_betas_minority_dist);
    let majority_dist_likelihoods: CellLogLikelihoodData = get_cell_log_likelihoods(&loci_used_for_posteriors, cell_data, &alpha_betas_majority_dist);
    let log_prior_minority: f64 = minority_fraction.ln();
    let log_prior_majority: f64 = (1.0 - minority_fraction).ln();
    for cell_id in 0..cell_data.len() {
        let log_numerator = log_prior_minority + minority_dist_likelihoods.log_likelihoods[cell_id];
        let log_denominator = log_numerator.ln_add_exp(log_prior_majority + majority_dist_likelihoods.log_likelihoods[cell_id]);
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

fn compute_new_excluded(params: &Params, loci_used: &mut Vec<bool>, cell_data: &Vec<CellData>, locus_counts: &Vec<[f64; 2]>, excluded_cells: &HashSet<usize>, iteration: usize) -> (bool, HashSet<usize>) {
    let alpha_betas = init_alpha_betas(locus_counts, excluded_cells, cell_data);
    let mut new_excluded: HashSet<usize> = HashSet::new();

    let cell_log_likelihood_data: CellLogLikelihoodData = get_cell_log_likelihoods(loci_used, cell_data, &alpha_betas);
    let mut normalized_log_likelihoods: Vec<f64> = Vec::new();
    for i in 0..cell_data.len() {
        normalized_log_likelihoods.push(cell_log_likelihood_data.log_likelihoods[i] / cell_log_likelihood_data.loci_used_per_cell[i]);
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
    // TODO also compute per loci per cell log likelihoods and check for extreme outliers to remove from loci_used
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
    let filename = format!("{}.iteration_{}_threshold.tsv",params.output_directory, iteration);
    let filehandle = File::create(&filename).expect(&format!("Unable to create file {}", &filename));
    let mut writer = BufWriter::new(filehandle);
    let line = format!("{}",threshold);
    writer.write_all(line.as_bytes()).expect("could not write to threshold file");
}
    
struct CellLogLikelihoodData {
    log_likelihoods: Vec<f64>,
    loci_used_per_cell: Vec<f64>,
    expected_log_likelihoods: Vec<f64>,
}

fn get_cell_log_likelihoods(loci_used: &Vec<bool>, cell_data: &Vec<CellData>, alpha_betas: &Vec<AlphaBeta>) -> 
        CellLogLikelihoodData {
    let mut cell_log_likelihoods: Vec<f64> = Vec::new();
    let mut cell_loci_used: Vec<f64> = Vec::new();
    let mut cell_expected_log_likelihood: Vec<f64> = Vec::new();
    for cell_id in 0..cell_data.len() {
        let mut log_likelihood: f64 = 0.0;
        let mut expected_log_likelihood: f64 = 0.0;
        let mut loci_used_for_cell: f64 = 0.0;
        let cell = &cell_data[cell_id];
        for locus in &cell.cell_loci_data {
            if loci_used[locus.locus_index] {
                log_likelihood += stats::log_beta_binomial_pmf(locus.alt_count, locus.ref_count, alpha_betas[locus.locus_index].alpha, alpha_betas[locus.locus_index].beta, locus.log_binomial_coefficient);
                expected_log_likelihood += stats::expected_log_beta_binomial_pmf(locus.total, alpha_betas[locus.locus_index].alpha, alpha_betas[locus.locus_index].beta, locus.log_binomial_coefficient);
                loci_used_for_cell += 1.0;
            }
        }
        cell_log_likelihoods.push(log_likelihood);
        cell_loci_used.push(loci_used_for_cell);
        cell_expected_log_likelihood.push(expected_log_likelihood);
    }
    return CellLogLikelihoodData{
        log_likelihoods: cell_log_likelihoods, 
        loci_used_per_cell: cell_loci_used,
        expected_log_likelihoods: cell_expected_log_likelihood};
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
    let posterior_threshold = params.value_of("posterior_threshold").unwrap_or("0.999");
    let posterior_threshold = posterior_threshold.to_string().parse::<f64>().unwrap();
    let interquartile_range_multiple = params.value_of("interquartile_range_multiple").unwrap_or("4");
    let interquartile_range_multiple = interquartile_range_multiple.to_string().parse::<f64>().unwrap();
    let output_directory = params.value_of("output_directory").unwrap().to_string();
    let min_alleles_posterior = params.value_of("min_alleles_posterior").unwrap_or("5");
    let min_alleles_posterior = min_alleles_posterior.to_string().parse::<usize>().unwrap();

    let params = Params {
        ref_mtx: ref_mtx,
        alt_mtx: alt_mtx,
        barcodes: barcodes,
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


