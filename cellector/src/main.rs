#[macro_use]
extern crate clap;
extern crate hashbrown;
extern crate rand;
extern crate statrs;
extern crate itertools;
extern crate rayon;
extern crate vcf;
extern crate flate2;



extern crate csv;
use csv::Writer;
use std::error::Error;
use std::f64::consts::E;


use std::f32::consts::PI;

use std::fs::OpenOptions;
use std::io::prelude::*;


use flate2::read::GzDecoder;
use flate2::read::MultiGzDecoder;
use vcf::*;


use rayon::prelude::*;

use rand::Rng;
use rand::rngs::StdRng;
use rand::SeedableRng;

use clap::App;
use std::f32;

use std::ffi::OsStr;
use std::io::Read;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use hashbrown::{HashMap,HashSet};
use itertools::izip;

use statrs::distribution::Exp;
use rand::distributions::Distribution;





fn main() {

    let params = load_params();

    let (ground_truth_barcode_to_assignment, ground_truth_barcode_and_assignments) = load_ground_truth(&params).unwrap();

    let cell_barcodes = load_barcodes(&params);

    let (loci_used, cells, coefficients, locus_counts, i) = load_cell_data(&params, ground_truth_barcode_to_assignment.clone(), cell_barcodes.clone());
    
    cellector(loci_used, cells, coefficients, locus_counts);

}




fn cellector( loci_used: usize, cell_data: Vec<CellStruct>, coefficients: Vec<Vec<f64>>, locus_counts: Vec<Vec<[usize; 3]>>) {

    let mut alphas_betas_pairs = init_alphas_betas_pairs(loci_used,&locus_counts);

    let mut normalized_log_likelihoods: Vec<f64> = vec![0.0; cell_data.len()];


    let mut minority_cluster_cells: HashSet<usize> = HashSet::new();


    normalized_log_likelihoods = calculate_normalized_log_likelihoods(&locus_counts, &coefficients, alphas_betas_pairs.clone(), &cell_data);
    let (mut minority_cluster, mut majority_cluster, minority_counts) = iqr_detector(&normalized_log_likelihoods);

    for &cell_index in minority_cluster.iter() {
        minority_cluster_cells.insert(cell_index);
    }

    let mut writer = Writer::from_path("output_iteration_1.csv").unwrap();
    writer.write_record(&["barcode", "ground_truth", "likelihood", "loci"]).unwrap();

    for (i, cell) in cell_data.iter().enumerate() {

        let record = vec![
            cell.barcode.clone(),
            cell.ground_truth.clone(),
            normalized_log_likelihoods[i].to_string(),
            locus_counts[i].len().to_string(),
        ];
    
        writer.write_record(&record).unwrap();
    }
    println!("outlierdetection...");
    println!("first iteration done...");

    
    let mut iteration_count = 2;

    loop {
        println!("Iteration {} done...", iteration_count);

        let new_alphas_betas_pairs = reset_alpha_beta_pairs(loci_used, &locus_counts, majority_cluster.clone());

        
        let new_normalized_log_likelihoods = calculate_normalized_log_likelihoods_filtered(
            &locus_counts, 
            &coefficients, 
            new_alphas_betas_pairs, 
            &cell_data, 
            majority_cluster.clone(), 
            normalized_log_likelihoods.clone()
        );
        
        let (mut minority_cluster_filtered, mut majority_cluster_filtered, minority_counts) = iqr_detector_filterd(&new_normalized_log_likelihoods, majority_cluster.clone());
        for &cell_index in minority_cluster_filtered.iter() {
            minority_cluster_cells.insert(cell_index);
        }

        let mut writer = Writer::from_path(format!("output_iteration_{}.csv", iteration_count)).unwrap();
        writer.write_record(&["barcode", "ground_truth", "likelihood", "loci"]).unwrap();
        for (i, cell) in cell_data.iter().enumerate() {

            let record = vec![
                cell.barcode.clone(),
                cell.ground_truth.clone(),
                new_normalized_log_likelihoods[i].to_string(),
                locus_counts[i].len().to_string(),
            ];
        
            writer.write_record(&record).unwrap();
        }


        normalized_log_likelihoods = new_normalized_log_likelihoods;
        majority_cluster = majority_cluster_filtered;

        if minority_counts == 0 {
            break;
        }

        iteration_count += 1;
    }


    posterioir_calc(loci_used, cell_data, coefficients, locus_counts, minority_cluster_cells, majority_cluster);



}

fn posterioir_calc(loci_used: usize, cell_data: Vec<CellStruct>, coefficients: Vec<Vec<f64>>,
     locus_counts: Vec<Vec<[usize; 3]>>, minority_cluster: HashSet<usize>, majority_cluster: Vec<usize>){

    let (minority_alpha_beta_pairs, majority_alpha_beta_pairs) = alpha_beta_pars_for_majority_and_minority(loci_used, majority_cluster.clone(), minority_cluster.clone(), locus_counts.clone());
    let mut cells_log_likelihoods_for_posterior: Vec<(f64, f64)> = vec![(0.0, 0.0); cell_data.len()];

    let excluded_cells_fraction = minority_cluster.len() as f64 / cell_data.len() as f64;

    
    // loop over all the cells
    for (cell_index, cell) in cell_data.iter().enumerate(){
        
        let mut cell_log_likelihood = 0.0;
        let mut cell_log_likelihood_majority = 0.0;
        let mut cell_log_likelihood_minority = 0.0;
        let cell_loci_counts = locus_counts[cell_index].clone(); 
        
        for (inner_locus_index, locus_counts) in cell_loci_counts.iter().enumerate() {
            
            let locus_index = locus_counts[0];
            let ref_count = locus_counts[1];
            let alt_count = locus_counts[2];
            let coefficient = coefficients[cell_index][inner_locus_index];
            let alpha_minority = minority_alpha_beta_pairs[0][locus_index];
            let beta_minority = minority_alpha_beta_pairs[1][locus_index];
            let alpha_majority = majority_alpha_beta_pairs[0][locus_index];
            let beta_majority = majority_alpha_beta_pairs[1][locus_index];
            let cell_log_likelihood_locus_minority = log_beta_binomial_PMF_for_posterior(alt_count, ref_count, alpha_minority as f64 * excluded_cells_fraction, beta_minority as f64 * excluded_cells_fraction, coefficient);
            let cell_log_likelihood_locus_majority = log_beta_binomial_PMF(alt_count, ref_count, alpha_majority, beta_majority, coefficient);   
            cell_log_likelihood_minority += cell_log_likelihood_locus_minority;
            cell_log_likelihood_majority += cell_log_likelihood_locus_majority;

        }


        // let normalized_likelihood_minority = cell_log_likelihood_minority / cell_loci_counts.len() as f64;
        // let normalized_likelihood_majority = cell_log_likelihood_majority / cell_loci_counts.len() as f64;
        // cells_log_likelihoods_for_posterior[cell_index] = (normalized_likelihood_minority, normalized_likelihood_majority);#FLAG

        cells_log_likelihoods_for_posterior[cell_index] = (cell_log_likelihood_minority, cell_log_likelihood_majority);
        


    }


    let minority_prior = excluded_cells_fraction;
    let majority_prior = 1.0 - excluded_cells_fraction;
    let log_minority_prior = minority_prior.ln();
    let log_majority_prior = majority_prior.ln();


    let mut denoms: Vec<f64> = vec![0.0; cell_data.len()];

    for (cell_index, cell) in cell_data.iter().enumerate() {

        let (minority_likelihood, majority_likelihood) = cells_log_likelihoods_for_posterior[cell_index];
        let minority_posterior_log = minority_likelihood + log_minority_prior;
        let majority_posterior_log = majority_likelihood + log_majority_prior;
        denoms[cell_index] = logsumexp(minority_posterior_log, majority_posterior_log);
        let posterior_minority = minority_posterior_log - denoms[cell_index];
        let posterior_majority = majority_posterior_log - denoms[cell_index];
        cells_log_likelihoods_for_posterior[cell_index] = (posterior_minority, posterior_majority);


    }

    let treshhold = 0.999;


    let writer = OpenOptions::new().write(true).create(true).open("cellector_assignments.tsv").unwrap();
    let mut writer = BufWriter::new(writer);
    writer.write_all(b"barcode\tground_truth\tassignment\n").unwrap();

    let mut final_assignments: Vec<String> = vec![String::new(); cell_data.len()];

    let mut cells_posteriors = vec![(0.0, 0.0); cell_data.len()];
    for (cell_index, cell) in cell_data.iter().enumerate() {
        let (minority_posterior, majority_posterior) = cells_log_likelihoods_for_posterior[cell_index];
        let minority_posterior = minority_posterior.exp();
        let majority_posterior = majority_posterior.exp();
        cells_posteriors[cell_index] = (minority_posterior, majority_posterior);

        let assignment = if minority_posterior > treshhold {
            "0"
        } else if majority_posterior > treshhold {
            "1"
        } else {
            "Unassigned"
        };
        writer.write_all(format!("{}\t{}\t{}\n", cell.barcode, cell.ground_truth, assignment).as_bytes()).unwrap();
        final_assignments[cell_index] = assignment.to_string();

    }

    let mut barcode_ground_truth_assignment: Vec<Vec<String>> = Vec::new();

    for cell_index in 0..cell_data.len() {
        let assignment = final_assignments[cell_index].clone();
        let ground_truth = cell_data[cell_index].ground_truth.clone();

        barcode_ground_truth_assignment.push(vec![cell_data[cell_index].barcode.clone(), ground_truth.clone(), assignment.clone()]);
        println!("{}\t{}\t{}", cell_data[cell_index].barcode, ground_truth, assignment);
    }





    let mut locus_to_plot: Vec<Vec<f64>> = vec![vec![0.0; 9]; loci_used];
    for (cell_index, cell) in cell_data.iter().enumerate() {

        let cellector_assignments = final_assignments[cell_index].clone();
        let cell_loci_counts = locus_counts[cell_index].clone();



        for (inner_locus_index, locus_counts) in cell_loci_counts.iter().enumerate() {

            let locus_index = locus_counts[0];
            let ref_count = locus_counts[1];
            let alt_count = locus_counts[2];

            if cellector_assignments == "0" {
                locus_to_plot[locus_index][0] += 1 as f64;
                locus_to_plot[locus_index][1] += ref_count as f64;
                locus_to_plot[locus_index][2] += alt_count as f64;
            } else if cellector_assignments == "1" {
                locus_to_plot[locus_index][3] += 1 as f64;
                locus_to_plot[locus_index][4] += ref_count as f64;
                locus_to_plot[locus_index][5] += alt_count as f64;
            } else {
                locus_to_plot[locus_index][6] += 1 as f64;
                locus_to_plot[locus_index][7] += ref_count as f64;
                locus_to_plot[locus_index][8] += alt_count as f64;
            }

        }


    }


    

    
        

}


fn logsumexp(val_1: f64, val_2: f64) -> f64 {
    let max = val_1.max(val_2);
    let sum = (val_1 - max).exp() + (val_2 - max).exp();
    max + sum.ln()

}

// fn logsumexp(values: &[f64]) -> f64 {
//     let max = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
//     let sum = values.iter().map(|x| (x - max).exp()).sum::<f64>();
//     max + sum.ln()
// }







fn alpha_beta_pars_for_majority_and_minority( loci_used: usize, majority_cluster: Vec<usize>, minority_cluster: HashSet<usize>, locus_counts: Vec<Vec<[usize; 3]>>
) -> (Vec<Vec<usize>>, Vec<Vec<usize>>){

    let mut minority_alpha_beta_pairs: Vec<Vec<usize>> = vec![vec![1; loci_used], vec![1; loci_used]];
    let mut majority_alpha_beta_pairs: Vec<Vec<usize>> = vec![vec![1; loci_used], vec![1; loci_used]];


    // minority 
    for cell_index in minority_cluster.iter() {

        let cell_loci = locus_counts[*cell_index].clone();

        for locus in cell_loci.iter() {

            let locus_index = locus[0];
            let ref_count = locus[1];
            let alt_count = locus[2];
            minority_alpha_beta_pairs[0][locus_index] += alt_count;
            minority_alpha_beta_pairs[1][locus_index] += ref_count;

        }

    }

    // Majority
    for cell_index in majority_cluster.iter() {

        let cell_loci = locus_counts[*cell_index].clone();

        for locus in cell_loci.iter() {

            let locus_index = locus[0];
            let ref_count = locus[1];
            let alt_count = locus[2];
            majority_alpha_beta_pairs[0][locus_index] += alt_count;
            majority_alpha_beta_pairs[1][locus_index] += ref_count;


        }

    }

    (minority_alpha_beta_pairs, majority_alpha_beta_pairs)

}







fn calculate_normalized_log_likelihoods_filtered
    (locus_counts: &Vec<Vec<[usize; 3]>>, coefficients: &Vec<Vec<f64>>, 
    alphas_betas_pairs: Vec<Vec<usize>>, cell_data: &Vec<CellStruct>, majority_cluster: Vec<usize>, normalized_log_likelihoods: Vec<f64>
)-> Vec<f64>{

    let mut updated_normalized_log_likelihoods = normalized_log_likelihoods.clone();

    for i in majority_cluster.iter() {
        updated_normalized_log_likelihoods[*i] = 0.0;
    }

    for (cell_index, cells_locus_counts) in locus_counts.iter().enumerate() {

        if majority_cluster.contains(&cell_index) {

            let mut normalized_likelihood = 0.0;
            let mut cell_log_likelihood = 0.0;
            let cell_loci = coefficients[cell_index].len();

            for (cell_locus_index, locus_count) in cells_locus_counts.iter().enumerate() {

                let locus_index = locus_count[0];
                let ref_count = locus_count[1];
                let alt_count = locus_count[2];
                let coefficient = coefficients[cell_index][cell_locus_index];
                let alpha = alphas_betas_pairs[0][locus_index];
                let beta = alphas_betas_pairs[1][locus_index];
                let cell_log_likelihood_locus = log_beta_binomial_PMF(alt_count, ref_count, alpha, beta, coefficient);
                cell_log_likelihood += cell_log_likelihood_locus;
            }
            // println!("barcode{}\t iterations \t likelihood{}", cell_data[cell_index].barcode, cell_log_likelihood);

            
            normalized_likelihood = cell_log_likelihood / cell_loci as f64;
            updated_normalized_log_likelihoods[cell_index] = normalized_likelihood;

        }
    }
    

    updated_normalized_log_likelihoods

}




fn calculate_normalized_log_likelihoods
    (locus_counts: &Vec<Vec<[usize; 3]>>, coefficients: &Vec<Vec<f64>>, alphas_betas_pairs: Vec<Vec<usize>>, cell_data: &Vec<CellStruct>)
     -> Vec<f64>{

    let mut normalized_log_likelihoods: Vec<f64> = vec![0.0; cell_data.len()];
    
    for (cell_index, cells_locus_counts) in locus_counts.iter().enumerate() {

        let mut normalized_likelihood = 0.0;
        let mut cell_log_likelihood = 0.0;
        let cell_loci = coefficients[cell_index].len();

        for (cell_locus_index, locus_count) in cells_locus_counts.iter().enumerate() {

            let locus_index = locus_count[0];
            let ref_count = locus_count[1];
            let alt_count = locus_count[2];
            let coefficient = coefficients[cell_index][cell_locus_index];
            let alpha = alphas_betas_pairs[0][locus_index];
            let beta = alphas_betas_pairs[1][locus_index];
            let cell_log_likelihood_locus = log_beta_binomial_PMF(alt_count, ref_count, alpha , beta, coefficient);
            cell_log_likelihood += cell_log_likelihood_locus;
        }

        if cell_log_likelihood.is_nan() || cell_log_likelihood.is_infinite() {
            println!("cell_log_likelihood is NaN or infinite");
        }
        normalized_likelihood = cell_log_likelihood / cell_loci as f64;
        normalized_log_likelihoods[cell_index] = normalized_likelihood;


    }


    normalized_log_likelihoods
}


fn reset_alpha_beta_pairs(loci_used: usize,locus_counts: &Vec<Vec<[usize; 3]>>, majority: Vec<usize>) ->  Vec<Vec<usize>>{

    let mut alphas_betas_pairs: Vec<Vec<usize>> = vec![vec![1; loci_used], vec![1; loci_used]];

    for cell_index in majority.iter() {

        let cell_loci = locus_counts[*cell_index].clone();

        for locus in cell_loci.iter() {

            let locus_index = locus[0];
            let ref_count = locus[1];
            let alt_count = locus[2];
            alphas_betas_pairs[0][locus_index] += alt_count;
            alphas_betas_pairs[1][locus_index] += ref_count;

        }

    }

    alphas_betas_pairs
}








fn iqr_detector_filterd(data: &Vec<f64>, majority_cluster_filtered: Vec<usize>) -> (Vec<usize>, Vec<usize>, usize) {
    
    let mut sorted_data = data.clone();

    sorted_data.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Greater));

    let mut minority_counts = 0;

    let q1 = percentile(&sorted_data, 25.0);
    let q3 = percentile(&sorted_data, 75.0);
    let iqr = q3 - q1;

    let lower_bound = q1 - 4.0 * iqr;
    let meadian = sorted_data[sorted_data.len() / 2];
    println!("median: {},IQR: {},  q1-4*iqr: {} len {}, max{}, min{}", meadian,iqr, lower_bound, sorted_data.len(), sorted_data[sorted_data.len()-1], sorted_data[0]);


    let mut minority_cluster = Vec::new();
    let mut majority_cluster = Vec::new();

    for i in majority_cluster_filtered.iter() {
        if data[*i] < lower_bound {
            minority_cluster.push(*i); 
            minority_counts += 1;
            
        } else {
            majority_cluster.push(*i); 
        }
    }

    println!("minority counts:new {}", minority_counts);
    (minority_cluster, majority_cluster, minority_counts)
}

fn iqr_detector(data: &Vec<f64>) -> (Vec<usize>, Vec<usize>, usize) {
    
    let mut sorted_data = data.clone();
    // sorted_data.sort_by(|a, b| a.partial_cmp(b).unwrap());
    sorted_data.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Greater));

    let mut minority_counts = 0;

    let q1 = percentile(&sorted_data, 25.0);
    let q3 = percentile(&sorted_data, 75.0);
    let iqr = q3 - q1;

    let lower_bound = q1 - 4.0 * iqr;

    let meadian = sorted_data[(sorted_data.len() / 2)+2];
    println!("median: {},IQR: {},  q1-4*iqr: {} len {}, max{}, min{}", meadian,iqr, lower_bound, sorted_data.len(), sorted_data[sorted_data.len()-1], sorted_data[0]);


    let mut minority_cluster = Vec::new();
    let mut majority_cluster = Vec::new();


    data.iter().enumerate().for_each(|(index, &value)| {
        if value < lower_bound {
            minority_cluster.push(index); 
            minority_counts += 1;
            
        } else {
            majority_cluster.push(index); 
        }
    });

    println!("minority counts:* {}", minority_counts);
    (minority_cluster, majority_cluster, minority_counts)
}

 
fn percentile(data: &Vec<f64>, percentile: f64) -> f64 {
    let k = (percentile / 100.0) * ((data.len() - 1) as f64);
    let f = k.floor() as usize;
    let c = k.ceil() as usize;

    if f != c {
        data[f] + (k - f as f64) * (data[c] - data[f])
    } else {
        data[f]
    }
}




fn log_beta_binomial_PMF_for_posterior(alt_count: usize, ref_count: usize, alpha: f64, beta: f64, coefficient: f64) -> f64 {

    let num_params = alt_count as f64 + alpha;
    let denom_params = ref_count as f64 + beta;
    let num = log_beta_calc_for_posterior(num_params, denom_params);
    let denom = log_beta_calc_for_posterior(alpha, beta);
    let log_likelihood = num as f64 - denom as f64; 
    log_likelihood

}


fn log_beta_binomial_PMF(alt_count: usize, ref_count: usize, alpha: usize, beta: usize, coefficient: f64) -> f64 {

    let num_params = alt_count + alpha;
    let denom_params = ref_count + beta;
    let num = log_beta_calc(num_params, denom_params);
    let denom = log_beta_calc(alpha, beta);
    let log_likelihood = (coefficient as f64) + num as f64 - denom as f64; 
    log_likelihood

}


fn init_alphas_betas_pairs(loci_used: usize,locus_counts: &Vec<Vec<[usize; 3]>>) ->  Vec<Vec<usize>>{

    let mut alphas_betas_pairs: Vec<Vec<usize>> = vec![vec![1; loci_used], vec![1; loci_used]];

    for (cell_index, cell_loci) in locus_counts.iter().enumerate() {

        for locus in cell_loci.iter() {

            let locus_index = locus[0];
            let ref_count = locus[1];
            let alt_count = locus[2];
            alphas_betas_pairs[0][locus_index] += alt_count;
            alphas_betas_pairs[1][locus_index] += ref_count;

        }

    }

    alphas_betas_pairs
}


fn calculate_variance(values: &[f64], mean: f64) -> f64 {
    if values.len() <= 1 {
        return 0.0; 
    }

    let mut sum_of_squared_diffs = 0.0;
    for &value in values {
        let diff = value - mean;
        sum_of_squared_diffs += diff * diff;
    }

    let variance = sum_of_squared_diffs / (values.len() - 1) as f64;
    variance
}

fn mean_calc(data: Vec<f64>) -> f64 {
    let sum: f64 = data.iter().sum();
    let mean = sum / data.len() as f64;
    mean
}


fn log_beta_calc_for_posterior(alpha: f64, beta: f64) -> f64 {

    let log_gamma_alpha = statrs::function::gamma::ln_gamma(alpha);
    let log_gamma_beta = statrs::function::gamma::ln_gamma(beta);
    let log_gamma_alpha_beta = statrs::function::gamma::ln_gamma(alpha + beta);

    log_gamma_alpha + log_gamma_beta - log_gamma_alpha_beta
}

fn log_beta_calc(alpha: usize, beta: usize) -> f64 {

    let log_gamma_alpha = statrs::function::gamma::ln_gamma(alpha as f64);
    let log_gamma_beta = statrs::function::gamma::ln_gamma(beta as f64);
    let log_gamma_alpha_beta = statrs::function::gamma::ln_gamma((alpha + beta) as f64);

    log_gamma_alpha + log_gamma_beta - log_gamma_alpha_beta
}








fn load_params() -> Params{

    let yaml = load_yaml!("params.yml");
    let params = App::from_yaml(yaml).get_matches();
    let ref_mtx = params.value_of("ref").unwrap().to_string();
    let alt_mtx = params.value_of("alt").unwrap().to_string();
    let barcodes = params.value_of("barcodes").unwrap().to_string();
    let min_alt = params.value_of("min_alt").unwrap_or("4");
    let min_alt = min_alt.to_string().parse::<u32>().unwrap();
    let min_ref = params.value_of("min_ref").unwrap_or("4");
    let min_ref = min_ref.to_string().parse::<u32>().unwrap();
    let ground_truth = params.value_of("ground_truth").unwrap().to_string();
    let outlier_detector = params.value_of("outlier_detector").unwrap_or("iqr");


    Params {

        ref_mtx: ref_mtx,
        alt_mtx: alt_mtx,
        barcodes: barcodes,
        min_alt: min_alt,
        min_ref: min_ref,
        ground_truth: ground_truth, 
        outlier_detector: outlier_detector.to_string(),

    }
}



fn load_cell_data(params: &Params, ground_truth_barcode_to_assignment: HashMap<String, String>, cell_barcodes: Vec<String>) 
-> (usize, Vec<CellStruct>, Vec<Vec<f64>>, Vec<Vec<[usize; 3]>>, HashMap<usize, usize>){

    let alt_reader = File::open(params.alt_mtx.to_string()).expect("cannot open alt mtx file");
    let ref_reader = File::open(params.ref_mtx.to_string()).expect("cannot open ref mtx file");

    let alt_reader = BufReader::new(alt_reader);
    let ref_reader = BufReader::new(ref_reader);


    let mut total_loci = 0;
    let mut total_cells = 0;
    let mut line_number = 0;
    let mut locus_to_locus_index: HashMap<usize, usize> = HashMap::new();
    let mut locus_alt_ref_cell_counts: HashMap<usize, [u32; 2]> = HashMap::new();
    let mut unique_cell_ids: HashSet<usize> = HashSet::new();
    let mut cell_index = 0;
    let mut cell_id_to_barcode_index: HashMap<usize, usize> = HashMap::new();
    let mut cell_id_to_locus_counts: HashMap<usize, Vec<[usize; 3]>> = HashMap::new();


    

    //Loading loci
    eprintln!("loading loci...");
    for (alt_line, ref_line) in izip!(alt_reader.lines(), ref_reader.lines()) {
        let alt_line = alt_line.expect("cannot read alt mtx");
        let ref_line = ref_line.expect("cannot read ref mtx");
        if line_number > 2 {

            let alt_tokens: Vec<&str> = alt_line.split_whitespace().collect();
            let ref_tokens: Vec<&str> = ref_line.split_whitespace().collect();

            let locus = alt_tokens[0].to_string().parse::<usize>().unwrap() - 1; // 0-indexed Locus_ID
            let ref_count = ref_tokens[2].to_string().parse::<u32>().unwrap();
            let alt_count = alt_tokens[2].to_string().parse::<u32>().unwrap();

            assert!(locus < total_loci);

            let temp_count = locus_alt_ref_cell_counts.entry(locus).or_insert([0; 2]);
            if ref_count > 0 { temp_count[0] += 1; }
            if alt_count > 0 { temp_count[1] += 1; }


            let cell = alt_tokens[1].to_string().parse::<usize>().unwrap() - 1; // 0-indexed Cell_ID
            assert!(cell < total_cells);
            unique_cell_ids.insert(cell);

            if ref_count + alt_count > 0 { 
            cell_id_to_locus_counts.entry(cell).or_insert(Vec::new()).push([locus, ref_count as usize, alt_count as usize]);
            }
            
        } else if line_number == 2 {
            let tokens: Vec<&str> = alt_line.split_whitespace().collect();
            total_loci = tokens[0].to_string().parse::<usize>().unwrap();
            total_cells = tokens[1].to_string().parse::<usize>().unwrap();
        }
        line_number += 1;
    }

    eprintln!("filtering loci...");

    let mut locus_index = 0;
    let mut index_to_locus_map : HashMap<usize, usize> = HashMap::new();

    for (locus, counts) in &locus_alt_ref_cell_counts {
        if counts[0] >= params.min_ref && counts[1] >= params.min_alt {
            locus_to_locus_index.insert(*locus, locus_index);
            index_to_locus_map.insert(locus_index, *locus);
            locus_index += 1;
        }
    }

    eprintln!("loading cell data...");

    


    let mut sorted_cell_ids: Vec<usize> = Vec::new();
    for cell_id in unique_cell_ids.iter() {
        sorted_cell_ids.push(*cell_id);
    }
    sorted_cell_ids.sort();



    let mut cell_structs_vec: Vec<CellStruct> = Vec::new();
    let mut cell_index_to_coefficients: Vec<Vec<f64>> = Vec::new(); 
    let mut cell_index_to_locus_counts: Vec<Vec<[usize; 3]>> = Vec::new(); 


    let mut cell_id_to_cell_index : HashMap<usize, usize> = HashMap::new();



    for (cell_index, cell_id) in sorted_cell_ids.iter().enumerate() {

        cell_id_to_cell_index.insert(*cell_id, cell_index);
        let barcode = cell_barcodes[cell_index].clone();
        let ground_truth = ground_truth_barcode_to_assignment.get(&barcode).map_or_else(|| "N".to_string(), |s| s.clone());
        cell_structs_vec.push(CellStruct::new(*cell_id, barcode, ground_truth));


        let mut alt_ref_counts: Vec<[usize; 3]> = Vec::new();
        let mut coefficients: Vec<f64> = Vec::new();

        let locus_counts = cell_id_to_locus_counts.get(cell_id).unwrap();
        for (inner_cell_locus_index, locus_count) in locus_counts.iter().enumerate() {

            
            if locus_to_locus_index.contains_key(&locus_count[0]) {

                let locus = locus_to_locus_index.get(&locus_count[0]).unwrap();
                let ref_count = locus_count[1];
                let alt_count = locus_count[2];
                let coefficient = statrs::function::factorial::ln_binomial((ref_count + alt_count) as u64, alt_count as u64) as f64;
                coefficients.push(coefficient);
                alt_ref_counts.push([*locus, ref_count, alt_count]);

            }
        
        }         

        
        cell_index_to_locus_counts.push(alt_ref_counts);
        cell_index_to_coefficients.push(coefficients);
    }


    (locus_to_locus_index.len(), cell_structs_vec, cell_index_to_coefficients, cell_index_to_locus_counts, index_to_locus_map)


}










#[derive(Clone)]
struct CellStruct {
    barcode: String,
    ground_truth: String,
    cell_id: usize,
}

impl CellStruct {
    fn new(ID: usize, barcode: String, ground_truth: String) -> CellStruct {
        CellStruct{
            cell_id: ID,
            barcode: barcode,
            ground_truth: ground_truth,
        }
    }
}
         



#[derive(Clone)]
struct CellData {
    cell_id: usize,
    barcode: String,
    allele_fractions: Vec<f32>,
    log_binomial_coefficient: Vec<f64>,
    alt_counts: Vec<u32>,
    ref_counts: Vec<u32>, 
    loci: Vec<usize>,
    total_alleles: f32,

}

impl CellData {
    fn new(ID: usize) -> CellData {
        CellData{
            cell_id: ID,
            barcode: String::new(),
            allele_fractions: Vec::new(),
            log_binomial_coefficient: Vec::new(),
            alt_counts: Vec::new(),
            ref_counts: Vec::new(),
            loci: Vec::new(),
            total_alleles: 0.0,
        }
    }
}





fn load_barcodes(params: &Params) -> Vec<String> {

    let reader = reader(&params.barcodes);
    let mut cell_barcodes: Vec<String> = Vec::new();
    for line in reader.lines() {
        let line = line.expect("Unable to read line");
        cell_barcodes.push(line.to_string());
    }
    cell_barcodes
}



fn load_ground_truth(params: &Params) -> io::Result<(HashMap<String, String>, Vec<Vec<String>>)> {

    eprintln!("Loading ground truth");
    let reader = reader(&params.ground_truth);
    let mut ground_truth_map: HashMap<String, String> = HashMap::new();
    let mut ground_truth_vec: Vec<Vec<String>> = Vec::new();

    for line_result in reader.lines() {
        let line = line_result.expect("Unable to read a line in the ground truth file");
        let columns: Vec<&str> = line.split('\t').collect(); 

        if columns.len() == 2 {

            let barcode = columns[0].to_string();
            let assignment = columns[1].to_string();

            ground_truth_map.insert(barcode.clone(), assignment.clone());
            ground_truth_vec.push(vec![barcode, assignment]);

        } else {

            eprintln!("Invalid line format: {}", line);
            eprintln!("The correct format is: barcode\tassignment");
        }
    }

    Ok((ground_truth_map, ground_truth_vec))
}


pub fn reader(filename: &str) -> Box<dyn BufRead> {
    let path = Path::new(filename);
    let file = match File::open(&path) {
        Err(why) => panic!("couldn't open file {}", filename),
        Ok(file) => file,
    };
    if path.extension() == Some(OsStr::new("gz")) {
        Box::new(BufReader::with_capacity(128 * 1024, MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::with_capacity(128 * 1024, file))
    }
}




struct Params {
    ref_mtx: String,
    alt_mtx: String,
    barcodes: String,
    min_alt: u32,
    min_ref: u32,
    ground_truth: String,
    outlier_detector: String,
}


struct output{
    barcode: String,
    ground_truth: String,
    likelihood: f64,
    assignment: String,
    CellData: CellData,
}

impl output {
    fn new(barcode: String, ground_truth: String, likelihood: f64, assignment: String, CellData: CellData) -> output {
        output{
            barcode: barcode,
            ground_truth: ground_truth,
            likelihood: likelihood,
            assignment: assignment,
            CellData: CellData,
        }
    }
}


