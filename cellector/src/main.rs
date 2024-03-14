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

    let cell_barcodes = load_barcodes(&params);

    let (locus_to_cell_data ,loci_used, total_cells, cell_data, index_to_locus, locus_to_index, locus_to_cell_ids) = load_cell_data(&params,cell_barcodes.clone());
    
    let (ground_truth_barcode_to_assignment, ground_truth_barcode_and_assignments) = load_ground_truth(&params).unwrap();

    println!("loci used: {}", loci_used.len());
    println!("total cells: {}", total_cells);
    println!("cell barcodes: {}", cell_barcodes.len());
    println!("ground truth: {}", ground_truth_barcode_to_assignment.len());

    cellector(cell_data, loci_used, &locus_to_cell_ids, locus_to_cell_data, ground_truth_barcode_and_assignments, ground_truth_barcode_to_assignment, cell_barcodes, index_to_locus, locus_to_index);

}



fn cellector(cell_data: Vec<CellData>, loci_used: HashSet<usize>, locus_to_cell_ids: &HashMap<usize, Vec<usize>>, 
    locus_to_cell_data: HashMap<usize,
    HashMap<usize, [u32; 2]>>,
    ground_truth_vec: Vec<Vec<String>>, ground_truth_barcode_to_assignment: HashMap<String, String>,
    cell_barcodes: Vec<String>,
    index_to_locus: Vec<usize>,
    locus_to_index : HashMap<usize, usize>,
    ){

    // make sure the data is loaded correctly

    eprintln!("cellector...");

    let mut alphas_betas_pairs = init_alphas_betas_pairs(&cell_data, &loci_used, &locus_to_cell_data);

    let mut normalized_log_likelihoods: Vec<f64> = vec![0.0; cell_data.len()];


    let mut cell_index_to_outputs: HashMap<usize, output> = HashMap::new();



    let mut barcode_index = 0;


    for (i, cell) in cell_data.iter().enumerate() {

        let mut barcode = cell_barcodes[barcode_index].clone();
        let ground_truth = ground_truth_barcode_to_assignment.get(&barcode).unwrap().clone();
        
        let mut normalized_likelihood = 0.0;
        let mut cell_likelihood = 0.0;
        let assignment = String::from("unassigned");

        

        for (i, locus_index) in cell.loci.iter().enumerate() {

           
            let log_likelihood = log_beta_binomial_PMF(locus_index, alphas_betas_pairs.clone(), cell, i, index_to_locus.clone());


            if i == 0 {
                cell_likelihood = log_likelihood;
            }
            else {
            cell_likelihood = log_sum_exp(cell_likelihood, log_likelihood);
            }

            

        }

        normalized_likelihood = cell_likelihood / cell.total_alleles as f64;
        normalized_log_likelihoods[i] = normalized_likelihood;

        
        cell_index_to_outputs.insert(i, output::new(barcode.clone(), ground_truth.clone(), normalized_likelihood, assignment, cell.clone()));
        // println!("{}\tbarcode: {} ground_truth: {} likelihood: {}", i,barcode, ground_truth, normalized_likelihood);
        barcode_index += 1;
        // println!("{}\tbarcode: {} ground_truth: {} likelihood: {}", i,barcode, ground_truth, normalized_likelihood);
        let cell_hashing = ground_truth_barcode_to_assignment.get(&barcode).unwrap().clone();
        println!("({}){}\t{}\t{}", i, barcode, cell_hashing, normalized_likelihood);


    }

    println!("first iteration done\n");
    println!("removing sigma 4 outliers\n");

    let mean: f64 = normalized_log_likelihoods.iter().sum::<f64>() / normalized_log_likelihoods.len() as f64;
    let variance: f64 = normalized_log_likelihoods.iter()
        .map(|&value| {
            let diff = value - mean;
            diff * diff
        })
        .sum::<f64>() / (normalized_log_likelihoods.len() - 1) as f64;
    let std_dev = variance.sqrt();

     for i in 0..normalized_log_likelihoods.len() {
        if normalized_log_likelihoods[i] > mean + 4.0 * std_dev {
            println!("removing outlier: {}", i);
        }
     }





}

fn mean_calc(data: Vec<f64>) -> f64 {
    let sum: f64 = data.iter().sum();
    let mean = sum / data.len() as f64;
    mean
}

fn log_sum_exp(likelihood_1: f64, likelihood_2: f64) -> f64 {
    let max = likelihood_1.max(likelihood_2);
    let sum = (likelihood_1 - max).exp() + (likelihood_2 - max).exp();
    sum.ln() + max
}


fn log_beta_binomial_PMF(locus_index: &usize, alphas_betas_pairs: HashMap<usize, Vec<usize>>, cell: &CellData, i: usize, index_to_locus: Vec<usize>
) -> f64 {

    let mut log_likelihood = 0.0 as f64;

    let locus = index_to_locus[*locus_index];
    let alpha = alphas_betas_pairs.get(&locus).unwrap()[0];
    let beta = alphas_betas_pairs.get(&locus).unwrap()[1];

    let log_binomial_coefficient = cell.log_binomial_coefficient[i];
    let ref_count = cell.ref_counts[i];
    let alt_count = cell.alt_counts[i];
    
    let num_params = alt_count as usize + alpha ;
    let denom_params =  ref_count as usize + beta ;


    let num = log_beta_calc(num_params, denom_params);
    let denom = log_beta_calc(alpha, beta);

    log_likelihood = (log_binomial_coefficient as f64) + num - denom;


    log_likelihood
}



fn log_beta_calc(alpha: usize, beta: usize) -> f64 {

    let log_gamma_alpha = statrs::function::gamma::ln_gamma(alpha as f64);
    let log_gamma_beta = statrs::function::gamma::ln_gamma(beta as f64);
    let log_gamma_alpha_beta = statrs::function::gamma::ln_gamma((alpha + beta) as f64);

    log_gamma_alpha + log_gamma_beta - log_gamma_alpha_beta
}





fn init_alphas_betas_pairs(cell_data: &Vec<CellData>, loci_used: &HashSet<usize>, locus_to_cell_data: &HashMap<usize, HashMap<usize, [u32; 2]>>) -> HashMap<usize, Vec<usize>> {

    // let mut alphas_betas_pairs : Vec<Vec<usize>> = vec![vec![0; loci_used.len()]; 3];
    let mut alphas_betas_pairs : HashMap<usize, Vec<usize>> = HashMap::new();

    for locus in loci_used.iter() {
        

        let mut alpha = 1;
        let mut beta = 1;

        for (cell, counts) in locus_to_cell_data.get(locus).unwrap() {
            alpha += counts[1];
            beta += counts[0];
        }

        alphas_betas_pairs.insert(*locus, vec![alpha as usize, beta as usize]);

    }

    // for (i, locus) in loci_used.iter().enumerate() {

    //     alphas_betas_pairs[0][i] = locus.clone();

    //     let mut alpha = 1;
    //     let mut beta = 1;

    //     for (cell, counts) in locus_to_cell_data.get(locus).unwrap() {
    //         alpha += counts[1];
    //         beta += counts[0];
    //     }

    //     alphas_betas_pairs[1][i] = alpha as usize;
    //     alphas_betas_pairs[2][i] = beta as usize;
        
    // }

    alphas_betas_pairs
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


    Params {

        ref_mtx: ref_mtx,
        alt_mtx: alt_mtx,
        barcodes: barcodes,
        min_alt: min_alt,
        min_ref: min_ref,
        ground_truth: ground_truth, 

    }
}


fn load_cell_data(params: &Params, barcodes: Vec<String>) -> (HashMap<usize, HashMap<usize, [u32; 2]>>, HashSet<usize>, usize, Vec<CellData>, Vec<usize>, HashMap<usize, usize>, HashMap<usize, Vec<usize>>) {
    let alt_reader = File::open(params.alt_mtx.to_string()).expect("cannot open alt mtx file");

    let alt_reader = BufReader::new(alt_reader);
    let ref_reader = File::open(params.ref_mtx.to_string()).expect("cannot open ref mtx file");

    let mut locus_to_cell_ids : HashMap<usize, Vec<usize>> = HashMap::new();
    let mut locus_to_cell_index : HashMap<usize, usize> = HashMap::new();
    let mut cell_id_to_CellData : HashMap<usize, CellData> = HashMap::new();


    
    let ref_reader = BufReader::new(ref_reader);
    let mut used_loci: HashSet<usize> = HashSet::new();
    let mut line_number = 0;
    let mut total_loci = 0;
    let mut total_cells = 0;
    let mut all_loci: HashSet<usize> = HashSet::new();
    let mut locus_cell_counts: HashMap<usize, [u32; 2]> = HashMap::new();
    let mut locus_umi_counts: HashMap<usize, [u32; 2]> = HashMap::new();
    let mut locus_counts: HashMap<usize, HashMap<usize, [u32; 2]>> = HashMap::new();
    for (alt_line, ref_line) in izip!(alt_reader.lines(), ref_reader.lines()) {
        let alt_line = alt_line.expect("cannot read alt mtx");
        let ref_line = ref_line.expect("cannot read ref mtx");
        if line_number > 2 {
            let alt_tokens: Vec<&str> = alt_line.split_whitespace().collect();
            let ref_tokens: Vec<&str> = ref_line.split_whitespace().collect();
            let locus = alt_tokens[0].to_string().parse::<usize>().unwrap() - 1; // Locus_ID
            all_loci.insert(locus);
            let cell = alt_tokens[1].to_string().parse::<usize>().unwrap() - 1; // Cell_ID

            locus_to_cell_ids.entry(locus).or_insert(Vec::new()).push(cell);  


            let ref_count = ref_tokens[2].to_string().parse::<u32>().unwrap();
            let alt_count = alt_tokens[2].to_string().parse::<u32>().unwrap();
            assert!(locus < total_loci);
            assert!(cell < total_cells);
            let cell_counts = locus_cell_counts.entry(locus).or_insert([0; 2]);
            let umi_counts = locus_umi_counts.entry(locus).or_insert([0; 2]);
            if ref_count > 0 { cell_counts[0] += 1; umi_counts[0] += ref_count; }
            if alt_count > 0 { cell_counts[1] += 1; umi_counts[1] += alt_count; }
            let cell_counts = locus_counts.entry(locus).or_insert(HashMap::new());
            cell_counts.insert(cell, [ref_count, alt_count]);
        } else if line_number == 2 {
            let tokens: Vec<&str> = alt_line.split_whitespace().collect();
            total_loci = tokens[0].to_string().parse::<usize>().unwrap();
            total_cells = tokens[1].to_string().parse::<usize>().unwrap();
        }
        line_number += 1;
    }

    let mut all_loci2: Vec<usize> = Vec::new();
    for loci in all_loci {
        all_loci2.push(loci);
    }
    let mut all_loci = all_loci2;

    all_loci.sort();
    let mut index_to_locus: Vec<usize> = Vec::new();
    let mut locus_to_index: HashMap<usize, usize> = HashMap::new();
    let mut cell_data: Vec<CellData> = Vec::new();

    for _cell in 0..total_cells {
        cell_data.push(CellData::new(0));
    }

    let mut locus_index = 0;
    let mut barcode_index = 0;
    for locus in all_loci {
        let cell_counts = locus_cell_counts.get(&locus).unwrap();
        let umi_counts = locus_umi_counts.get(&locus).unwrap();
        if cell_counts[0] >= params.min_ref && cell_counts[1] >= params.min_alt && umi_counts[0] >= 0 && umi_counts[1] >= 0 {
            used_loci.insert(locus);
            index_to_locus.push(locus);
            locus_to_index.insert(locus, locus_index);
            for (cell, counts) in locus_counts.get(&locus).unwrap() {
                
                if counts[0]+counts[1] == 0 { continue; }
                
                cell_data[*cell].barcode = barcodes[barcode_index].clone();
                cell_data[*cell].cell_id = *cell;
                // barcode_index += 1;
                cell_data[*cell].alt_counts.push(counts[1]);
                cell_data[*cell].ref_counts.push(counts[0]);
                cell_data[*cell].loci.push(locus_index);
                cell_data[*cell].allele_fractions.push((counts[1] as f32)/((counts[0] + counts[1]) as f32));
                cell_data[*cell].log_binomial_coefficient.push(
                     statrs::function::factorial::ln_binomial((counts[1]+counts[0]) as u64, counts[1] as u64) as f32);
                cell_data[*cell].total_alleles += (counts[0] + counts[1]) as f32;


                

            }
            locus_index += 1;
        }
    }
    eprintln!("total loci used {}",used_loci.len());
    
    (locus_counts, used_loci, total_cells, cell_data, index_to_locus, locus_to_index, locus_to_cell_ids)
}




         



#[derive(Clone)]
struct CellData {
    cell_id: usize,
    barcode: String,
    allele_fractions: Vec<f32>,
    log_binomial_coefficient: Vec<f32>,
    alt_counts: Vec<u32>,
    ref_counts: Vec<u32>, // 
    loci: Vec<usize>,// ID
    total_alleles: f32,
}

impl CellData {
    fn new(ID: usize) -> CellData {
        CellData{
            cell_id: 0,
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
        let columns: Vec<&str> = line.split('\t').collect(); // Use tab as the separator

        if columns.len() == 2 {
            let barcode = columns[0].to_string();
            let assignment = columns[1].to_string();

            ground_truth_map.insert(barcode.clone(), assignment.clone());

            ground_truth_vec.push(vec![barcode, assignment]);
        } else {

            eprintln!("Invalid line format: {}", line);
            //add the correct forpat
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



