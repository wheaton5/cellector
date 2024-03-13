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
    let (locus_to_cell_data ,loci_used, total_cells, cell_data, index_to_locus, locus_to_index, locus_to_cell_ids) = load_cell_data(&params);
    let cell_barcodes = load_barcodes(&params);
    let (ground_truth_barcode_to_assignment, ground_truth_vec) = load_ground_truth(&params).unwrap();

    println!("loci used: {}", loci_used.len());
    println!("total cells: {}", total_cells);
    println!("cell barcodes: {}", cell_barcodes.len());
    println!("ground truth: {}", ground_truth_barcode_to_assignment.len());

    cellector(cell_data, loci_used, &locus_to_cell_ids, locus_to_cell_data);

}



fn cellector(cell_data: Vec<CellData>, loci_used: HashSet<usize>, locus_to_cell_ids: &HashMap<usize, Vec<usize>>, locus_to_cell_data: HashMap<usize, HashMap<usize, [u32; 2]>>){

    // make sure the data is loaded correctly

    eprintln!("cellector...");

    let mut alphas_betas_pairs = init_alphas_betas_pairs(&cell_data, &loci_used, &locus_to_cell_data);

    let mut log_likelihoods: Vec<f64> = vec![0.0; cell_data.len()]; 


    for (i, cell) in cell_data.iter().enumerate() {

        let likelihood = log_beta_binomial_PMF(alphas_betas_pairs[1][0], alphas_betas_pairs[2][0], cell.alt_counts[0] as usize, cell.alt_counts[0] as usize + cell.ref_counts[0] as usize);
        let normilized_liklihood = likelihood / ( cell.loci.len() as f64);  // normalize by the number of loci
        log_likelihoods[i] = normilized_liklihood;
    }





    let mean_likelihood = log_likelihoods.iter().sum::<f64>() / log_likelihoods.len() as f64;
    let sigma = log_likelihoods.iter().map(|x| (x - mean_likelihood).powi(2)).sum::<f64>() / log_likelihoods.len() as f64;
    let sigma = sigma.sqrt();




}




fn log_beta_binomial_PMF(alpha: usize, beta: usize, k: usize, n: usize) -> f64 {

    let mut log_likelihood = 0.0;

    for i in 0..k {

        log_likelihood += statrs::function::factorial::ln_binomial(n as u64, i as u64) as f64;
        log_likelihood += statrs::function::factorial::ln_binomial(alpha as u64, i as u64) as f64;
        log_likelihood += statrs::function::factorial::ln_binomial(beta as u64, (n-i) as u64) as f64;

    }

    log_likelihood

}



fn init_alphas_betas_pairs(cell_data: &Vec<CellData>, loci_used: &HashSet<usize>, locus_to_cell_data: &HashMap<usize, HashMap<usize, [u32; 2]>>) -> Vec<Vec<usize>>{

    let mut alphas_betas_pairs : Vec<Vec<usize>> = vec![vec![0; loci_used.len()]; 3];

    for (i, locus) in loci_used.iter().enumerate() {

        alphas_betas_pairs[0][i] = locus.clone();

        let mut alpha = 1;
        let mut beta = 1;

        for (cell, counts) in locus_to_cell_data.get(locus).unwrap() {

            alpha += counts[1];
            beta += counts[0];
            

        }

        alphas_betas_pairs[1][i] = alpha as usize;
        alphas_betas_pairs[2][i] = beta as usize;
        
    }

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


fn load_cell_data(params: &Params) -> (HashMap<usize, HashMap<usize, [u32; 2]>>, HashSet<usize>, usize, Vec<CellData>, Vec<usize>, HashMap<usize, usize>, HashMap<usize, Vec<usize>>) {
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
    for locus in all_loci {
        let cell_counts = locus_cell_counts.get(&locus).unwrap();
        let umi_counts = locus_umi_counts.get(&locus).unwrap();
        if cell_counts[0] >= params.min_ref && cell_counts[1] >= params.min_alt && umi_counts[0] >= 0 && umi_counts[1] >= 0 {
            used_loci.insert(locus);
            index_to_locus.push(locus);
            locus_to_index.insert(locus, locus_index);
            for (cell, counts) in locus_counts.get(&locus).unwrap() {
                
                if counts[0]+counts[1] == 0 { continue; }
                
                cell_data[*cell].alt_counts.push(counts[1]);
                cell_data[*cell].ref_counts.push(counts[0]);
                cell_data[*cell].loci.push(locus_index);
                cell_data[*cell].allele_fractions.push((counts[1] as f32)/((counts[0] + counts[1]) as f32));
                cell_data[*cell].log_binomial_coefficient.push(
                     statrs::function::factorial::ln_binomial((counts[1]+counts[0]) as u64, counts[1] as u64) as f32);
                cell_data[*cell].total_alleles += (counts[0] + counts[1]) as f32;

                // cell is cell index in the celldata, what I want is to be able to acces this by cell_id
                

            }
            locus_index += 1;
        }
    }
    eprintln!("total loci used {}",used_loci.len());
    
    (locus_counts, used_loci, total_cells, cell_data, index_to_locus, locus_to_index, locus_to_cell_ids)
}




         




struct CellData {
    cell_id: usize,
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





