#[macro_use]

extern crate hashbrown;
extern crate rand;
extern crate statrs;
extern crate itertools;
extern crate rayon;
extern crate vcf;
extern crate flate2;



//Reza
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


use std::f32;

use std::ffi::OsStr;
use std::io::Read;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use hashbrown::{HashMap,HashSet};
use itertools::izip;




extern crate clap;
use clap::App;





fn main() {

    let params = load_params();


}





fn load_params() -> Params{

    let yaml = load_yaml!("params.yml");
    let params = App::from_yaml(yaml).get_matches();

    let ref_mtx = params.value_of("ref_mtx").unwrap().to_string();
    let alt_mtx = params.value_of("alt_mtx").unwrap().to_string();
    let barcodes = params.value_of("barcodes").unwrap().to_string();
    let min_alt = params.value_of("min_alt").unwrap_or("4").parse::<u32>().unwrap();
    let min_ref = params.value_of("min_ref").unwrap_or("4").parse::<u32>().unwrap();
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


fn load_cell_data(params: &Params) {
  
    let alt_reader = File::open(params.alt_mtx.to_string()).expect("cannot open alt mtx file");
    let alt_reader = BufReader::new(alt_reader);
    let ref_reader = File::open(params.ref_mtx.to_string()).expect("cannot open ref mtx file");
    let ref_reader = BufReader::new(ref_reader);

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

        for (alt_line, ref_line) in izip!(alt_reader.lines(), ref_reader.lines()) {
            let alt_line = alt_line.expect("cannot read alt mtx");
            let ref_line = ref_line.expect("cannot read ref mtx");
            if line_number > 2 {
                let alt_tokens: Vec<&str> = alt_line.split_whitespace().collect();
                let ref_tokens: Vec<&str> = ref_line.split_whitespace().collect();
                let locus = alt_tokens[0].to_string().parse::<usize>().unwrap() - 1;
                all_loci.insert(locus);
                let cell = alt_tokens[1].to_string().parse::<usize>().unwrap() - 1;
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

struct CellData {
    allele_fractions: Vec<f32>,
    log_binomial_coefficient: Vec<f32>,
    alt_counts: Vec<u32>,
    ref_counts: Vec<u32>, 
    loci: Vec<usize>,
    total_alleles: f32,
}
