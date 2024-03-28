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

    let (ground_truth_barcode_to_assignment, ground_truth_barcode_and_assignments) = load_ground_truth(&params).unwrap();

    let cell_barcodes = load_barcodes(&params);

    let (loci_used, cells, coefficients, locus_counts, i) = load_cell_data(&params, ground_truth_barcode_to_assignment.clone(), cell_barcodes.clone());

    // let x = i.get(&10707).unwrap();
    // println!("{:?}", x);
    // let y = i.get(&6831).unwrap();
    // println!("{:?}", y);
    
    cellector(loci_used, cells, coefficients, locus_counts);

}




fn cellector( loci_used: usize, cell_data: Vec<CellStruct>, coefficients: Vec<Vec<f32>>, locus_counts: Vec<Vec<[usize; 3]>>) {

    let mut alphas_betas_pairs = init_alphas_betas_pairs(loci_used,&locus_counts);
    let mut normalized_log_likelihoods: Vec<f32> = vec![0.0; cell_data.len()];

    // wanna make an stack that willhold new minorities in each iteration which is just cell indices
    let mut all_minorities: Vec<usize> = Vec::new();

    // println!("loci used: {}", loci_used);
    // println!("cell data: {:?}", cell_data.len());
    // println!("coefficients: {:?}", coefficients.len());
    // println!("locus counts: {:?}", locus_counts.len());
    println!("\n\ncellector started...");

    normalized_log_likelihoods = calculate_normalized_log_likelihoods(&locus_counts, &coefficients, alphas_betas_pairs.clone(), &cell_data);
    let (mut minority_cluster, mut majority_cluster) = iqr_detector(&normalized_log_likelihoods);
    all_minorities.append(&mut minority_cluster);
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

    println!("first iteration done...\n");

    
    let mut iteration_count = 2;

    loop {
        println!("Iteration {} done...\n", iteration_count);

        let new_alphas_betas_pairs = reset_alpha_beta_pairs(loci_used, &locus_counts, majority_cluster.clone());
        
        let new_normalized_log_likelihoods = calculate_normalized_log_likelihoods_filtered(
            &locus_counts, 
            &coefficients, 
            new_alphas_betas_pairs, 
            &cell_data, 
            majority_cluster.clone(), 
            normalized_log_likelihoods.clone()
        );
        
        let (mut minority_cluster_filtered, mut majority_cluster_filtered) = iqr_detector(&new_normalized_log_likelihoods);
        all_minorities.append(&mut minority_cluster_filtered);

        // just wanna check if there is a new minority in each iteration
        

        //wanna wtite to file for each iteration similar to above
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

        // Prepare for next iteration or break loop
        normalized_log_likelihoods = new_normalized_log_likelihoods;
        majority_cluster = majority_cluster_filtered;



        // if minority is less than 3 break
        if iteration_count > 8 {
            break;
        }

        iteration_count += 1;
    }


}



fn calculate_normalized_log_likelihoods_filtered
    (locus_counts: &Vec<Vec<[usize; 3]>>, coefficients: &Vec<Vec<f32>>, 
    alphas_betas_pairs: Vec<Vec<usize>>, cell_data: &Vec<CellStruct>, majority_cluster: Vec<usize>, normalized_log_likelihoods: Vec<f32>
)-> Vec<f32>{

    let mut updated_normalized_log_likelihoods = normalized_log_likelihoods.clone();
    //set the majority indeces to 0 in likelihoods
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

            
            normalized_likelihood = cell_log_likelihood / cell_loci as f32;
            updated_normalized_log_likelihoods[cell_index] = normalized_likelihood;

        }
    }
    

    updated_normalized_log_likelihoods

}


fn calculate_normalized_log_likelihoods
    (locus_counts: &Vec<Vec<[usize; 3]>>, coefficients: &Vec<Vec<f32>>, alphas_betas_pairs: Vec<Vec<usize>>, cell_data: &Vec<CellStruct>)
     -> Vec<f32>{

    let mut normalized_log_likelihoods: Vec<f32> = vec![0.0; cell_data.len()];
    
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
            let cell_log_likelihood_locus = log_beta_binomial_PMF(alt_count, ref_count, alpha, beta, coefficient);
            cell_log_likelihood += cell_log_likelihood_locus;
        }

        if cell_log_likelihood.is_nan() || cell_log_likelihood.is_infinite() {
            println!("cell_log_likelihood is NaN or infinite");
        }
        normalized_likelihood = cell_log_likelihood / cell_loci as f32;
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



fn outlier_detector(data: &Vec<f32>, params: &Params) -> (Vec<usize>, Vec<usize>) {
    let (minority, majority) = if params.outlier_detector == "iqr" {
        iqr_detector(data)
    } else {
        panic!("Invalid outlier detector method specified.")
    };

    (minority, majority)
}

fn print_special_floats(data: &Vec<f32>) {
    for (index, &value) in data.iter().enumerate() {
        if value.is_nan() {
            println!("NaN found at index: {}", index);
        } else if value.is_infinite() {
            if value.is_sign_positive() {
                println!("Positive infinity found at index: {}", index);
            } else {
                println!("Negative infinity found at index: {}", index);
            }
        }
    }
}




fn iqr_detector(data: &Vec<f32>) -> (Vec<usize>, Vec<usize>) {
    
    let mut sorted_data = data.clone();
    // sorted_data.sort_by(|a, b| a.partial_cmp(b).unwrap());
    sorted_data.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Greater));


    let q1 = percentile(&sorted_data, 25.0);
    let q3 = percentile(&sorted_data, 75.0);
    let iqr = q3 - q1;

    let lower_bound = q1 - 4.0 * iqr;


    let mut minority_cluster = Vec::new();
    let mut majority_cluster = Vec::new();


    data.iter().enumerate().for_each(|(index, &value)| {
        if value < lower_bound {
            minority_cluster.push(index); 
        } else {
            majority_cluster.push(index); 
        }
    });

    (minority_cluster, majority_cluster)
}

 
fn percentile(data: &Vec<f32>, percentile: f32) -> f32 {
    let k = (percentile / 100.0) * ((data.len() - 1) as f32);
    let f = k.floor() as usize;
    let c = k.ceil() as usize;

    if f != c {
        data[f] + (k - f as f32) * (data[c] - data[f])
    } else {
        data[f]
    }
}







fn log_beta_binomial_PMF(alt_count: usize, ref_count: usize, alpha: usize, beta: usize, coefficient: f32) -> f32 {

    let num_params = alt_count + alpha;
    let denom_params = ref_count + beta;
    let num = log_beta_calc(num_params, denom_params);
    let denom = log_beta_calc(alpha, beta);
    let log_likelihood = (coefficient as f32) + num as f32 - denom as f32;
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
-> (usize, Vec<CellStruct>, Vec<Vec<f32>>, Vec<Vec<[usize; 3]>>, HashMap<usize, usize>){

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

            // obtaining cell_id to 1-generate index 2-sort in order to acces barcodes.
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
    let mut cell_index_to_coefficients: Vec<Vec<f32>> = Vec::new(); // cell_id to a vec of log binomial coefficients indexed based ob locus_index of the cell data
    let mut cell_index_to_locus_counts: Vec<Vec<[usize; 3]>> = Vec::new(); // cell_id to a vec of locus counts indexed based ob locus_index of the cell data


    let mut cell_id_to_cell_index : HashMap<usize, usize> = HashMap::new();



    for (cell_index, cell_id) in sorted_cell_ids.iter().enumerate() {

        cell_id_to_cell_index.insert(*cell_id, cell_index);
        let barcode = cell_barcodes[cell_index].clone();
        let ground_truth = ground_truth_barcode_to_assignment.get(&barcode).map_or_else(|| "N".to_string(), |s| s.clone());
        cell_structs_vec.push(CellStruct::new(*cell_id, barcode, ground_truth));


        let mut alt_ref_counts: Vec<[usize; 3]> = Vec::new();
        let mut coefficients: Vec<f32> = Vec::new();

        let locus_counts = cell_id_to_locus_counts.get(cell_id).unwrap();
        for (inner_cell_locus_index, locus_count) in locus_counts.iter().enumerate() {

            
            if locus_to_locus_index.contains_key(&locus_count[0]) {

                let locus = locus_to_locus_index.get(&locus_count[0]).unwrap();
                let ref_count = locus_count[1];
                let alt_count = locus_count[2];
                let coefficient = statrs::function::factorial::ln_binomial((ref_count + alt_count) as u64, alt_count as u64) as f32;
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
    log_binomial_coefficient: Vec<f32>,
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



//Previous verision:

// fn load_cell_data(params: &Params, barcodes: Vec<String>) -> (HashMap<usize, HashMap<usize, [u32; 2]>>, HashSet<usize>, usize, Vec<CellData>, Vec<usize>, HashMap<usize, usize>, HashMap<usize, Vec<usize>>) {
//     let alt_reader = File::open(params.alt_mtx.to_string()).expect("cannot open alt mtx file");

//     let alt_reader = BufReader::new(alt_reader);
//     let ref_reader = File::open(params.ref_mtx.to_string()).expect("cannot open ref mtx file");

//     let mut locus_to_cell_ids : HashMap<usize, Vec<usize>> = HashMap::new();
//     let mut locus_to_cell_index : HashMap<usize, usize> = HashMap::new();
//     let mut cell_id_to_CellData : HashMap<usize, CellData> = HashMap::new();


    
//     let ref_reader = BufReader::new(ref_reader);
//     let mut used_loci: HashSet<usize> = HashSet::new();
//     let mut line_number = 0;
//     let mut total_loci = 0;
//     let mut total_cells = 0;
//     let mut all_loci: HashSet<usize> = HashSet::new();
//     let mut locus_cell_counts: HashMap<usize, [u32; 2]> = HashMap::new();
//     let mut locus_umi_counts: HashMap<usize, [u32; 2]> = HashMap::new();
//     let mut locus_counts: HashMap<usize, HashMap<usize, [u32; 2]>> = HashMap::new();
//     for (alt_line, ref_line) in izip!(alt_reader.lines(), ref_reader.lines()) {
//         let alt_line = alt_line.expect("cannot read alt mtx");
//         let ref_line = ref_line.expect("cannot read ref mtx");
//         if line_number > 2 {
//             let alt_tokens: Vec<&str> = alt_line.split_whitespace().collect();
//             let ref_tokens: Vec<&str> = ref_line.split_whitespace().collect();
//             let locus = alt_tokens[0].to_string().parse::<usize>().unwrap() - 1; // Locus_ID
//             all_loci.insert(locus);
//             let cell = alt_tokens[1].to_string().parse::<usize>().unwrap() - 1; // Cell_ID

//             locus_to_cell_ids.entry(locus).or_insert(Vec::new()).push(cell);  


//             let ref_count = ref_tokens[2].to_string().parse::<u32>().unwrap();
//             let alt_count = alt_tokens[2].to_string().parse::<u32>().unwrap();
//             assert!(locus < total_loci);
//             assert!(cell < total_cells);
//             let cell_counts = locus_cell_counts.entry(locus).or_insert([0; 2]);
//             let umi_counts = locus_umi_counts.entry(locus).or_insert([0; 2]);
//             if ref_count > 0 { cell_counts[0] += 1; umi_counts[0] += ref_count; }
//             if alt_count > 0 { cell_counts[1] += 1; umi_counts[1] += alt_count; }
//             let cell_counts = locus_counts.entry(locus).or_insert(HashMap::new());
//             cell_counts.insert(cell, [ref_count, alt_count]);
//         } else if line_number == 2 {
//             let tokens: Vec<&str> = alt_line.split_whitespace().collect();
//             total_loci = tokens[0].to_string().parse::<usize>().unwrap();
//             total_cells = tokens[1].to_string().parse::<usize>().unwrap();
//         }
//         line_number += 1;
//     }

//     let mut all_loci2: Vec<usize> = Vec::new();
//     for loci in all_loci {
//         all_loci2.push(loci);
//     }
//     let mut all_loci = all_loci2;

//     all_loci.sort();
//     let mut index_to_locus: Vec<usize> = Vec::new();
//     let mut locus_to_index: HashMap<usize, usize> = HashMap::new();
//     let mut cell_data: Vec<CellData> = Vec::new();

//     for _cell in 0..total_cells {
//         cell_data.push(CellData::new(0));
//     }

//     let mut locus_index = 0;
//     let mut barcode_index = 0;
//     for locus in all_loci {
//         let cell_counts = locus_cell_counts.get(&locus).unwrap();
//         let umi_counts = locus_umi_counts.get(&locus).unwrap();
//         if cell_counts[0] >= params.min_ref && cell_counts[1] >= params.min_alt && umi_counts[0] >= 0 && umi_counts[1] >= 0 {
//             used_loci.insert(locus);
//             index_to_locus.push(locus);
//             locus_to_index.insert(locus, locus_index);
//             for (cell, counts) in locus_counts.get(&locus).unwrap() {
                
//                 if counts[0]+counts[1] == 0 { continue; }
                
//                 cell_data[*cell].barcode = barcodes[barcode_index].clone();
//                 cell_data[*cell].cell_id = *cell;
//                 // barcode_index += 1;
//                 cell_data[*cell].alt_counts.push(counts[1]);
//                 cell_data[*cell].ref_counts.push(counts[0]);
//                 cell_data[*cell].loci.push(locus_index);
//                 cell_data[*cell].allele_fractions.push((counts[1] as f32)/((counts[0] + counts[1]) as f32));
//                 cell_data[*cell].log_binomial_coefficient.push(
//                      statrs::function::factorial::ln_binomial((counts[1]+counts[0]) as u64, counts[1] as u64) as f32);
//                 cell_data[*cell].total_alleles += (counts[0] + counts[1]) as f32;


                

//             }
//             locus_index += 1;
//         }
//     }
//     eprintln!("total loci used {}",used_loci.len());

//     //wanna write down all the used loci in a file
//     let mut writer = Writer::from_path("loci.csv").unwrap();
//     // writer.write_record(&["locus"]).unwrap();
//     for locus in used_loci.iter() {
//         writer.write_record(&[locus.to_string().as_str()]).unwrap();
//     }

    
//     (locus_counts, used_loci, total_cells, cell_data, index_to_locus, locus_to_index, locus_to_cell_ids)
// }



// fn cellector_2(cell_data: Vec<CellData>, loci_used: HashSet<usize>, locus_to_cell_ids: &HashMap<usize, Vec<usize>>, 
//     locus_to_cell_data: HashMap<usize,
//     HashMap<usize, [u32; 2]>>,
//     ground_truth_vec: Vec<Vec<String>>, ground_truth_barcode_to_assignment: HashMap<String, String>,
//     cell_barcodes: Vec<String>,
//     index_to_locus: Vec<usize>,
//     locus_to_index : HashMap<usize, usize>,
//     ){

//     // make sure the data is loaded correctly

//     eprintln!("cellector...");

//     let mut alphas_betas_pairs = init_alphas_betas_pairs(&cell_data, &loci_used, &locus_to_cell_data);
//     let mut normalized_log_likelihoods: Vec<f64> = vec![0.0; cell_data.len()];
//     let mut cell_index_to_outputs: HashMap<usize, output> = HashMap::new();
//     let mut barcode_index = 0;


//     for (i, cell) in cell_data.iter().enumerate() {

//         let mut barcode = cell_barcodes[barcode_index].clone();
//         let ground_truth = ground_truth_barcode_to_assignment.get(&barcode).unwrap().clone();
        
//         let mut normalized_likelihood = 0.0;
//         let mut cell_log_likelihood = 0.0;
//         let assignment = String::from("unassigned");

        

//         for (j, locus_index) in cell.loci.iter().enumerate() {

//             // let locus = index_to_locus[*locus_index];
//             // if !loci_used.contains(&locus) {
//             //     continue;
//             // }

           
//             let log_likelihood = log_beta_binomial_PMF(locus_index, alphas_betas_pairs.clone(), cell, j, index_to_locus.clone());

//             assert!(log_likelihood <= 0.0);
            
//             cell_log_likelihood += log_likelihood;
           

            

//         }

//         normalized_likelihood = cell_log_likelihood / cell.alt_counts.len() as f64;
//         normalized_log_likelihoods[i] = normalized_likelihood;

        
//         cell_index_to_outputs.insert(i, output::new(barcode.clone(), ground_truth.clone(), normalized_likelihood, assignment, cell.clone()));
//         println!("{}\tbarcode: {} ground_truth: {} likelihood: {}", i,barcode, ground_truth, normalized_likelihood);
//         barcode_index += 1;
      


//     }


//     let mut writer = Writer::from_path("output.csv").unwrap();
//     writer.write_record(&["barcode", "ground_truth", "likelihood", "loci"]).unwrap();

//     for (i, output) in cell_index_to_outputs.iter() {

//         let loci = output.CellData.ref_counts.len().to_string();
//         writer.write_record(&[&output.barcode, &output.ground_truth, &output.likelihood.to_string(), &loci]).unwrap();
//     }

//     println!("first iteration done\n");
//     println!("removing sigma 4 outliers\n");

//     let mean: f64 = normalized_log_likelihoods.iter().sum::<f64>() / normalized_log_likelihoods.len() as f64;
//     let variance: f64 = calculate_variance(&normalized_log_likelihoods, mean);
//     let std_dev = variance.sqrt();

    
//     let mut num_removed = 0;
//     let mut removed_outliers: Vec<usize> = Vec::new();
//     let mut removing_step = 0;

//     loop{

//         if removing_step == 0{
//             for i in 0..normalized_log_likelihoods.len() {
//                 if normalized_log_likelihoods[i] > mean + 4.0 * std_dev {
//                     removed_outliers.push(i);
//                     num_removed += 1;
//                 }
//             }
//             removing_step = 1;
//         }
//         else{
//             // // recalculate the mean without removed cells
//             // let mut new_mean = 0.0;
//             // let mut new_variance = 0.0;
//             // let mut new_std_dev = 0.0;
//             // let mut new_normalized_log_likelihoods: Vec<f64> = Vec::new();
//             // for i in 0..normalized_log_likelihoods.len() {
//             //     if !removed_outliers.contains(&i) {
//             //         new_normalized_log_likelihoods.push(normalized_log_likelihoods[i]);
//             //     }
//             // }
//             // new_mean = new_normalized_log_likelihoods.iter().sum::<f64>() / new_normalized_log_likelihoods.len() as f64;
//             // new_variance = calculate_variance(&new_normalized_log_likelihoods, new_mean);
//             // new_std_dev = new_variance.sqrt();

//             // //remove new 4-sigma outliers
//             // for i in 0..new_normalized_log_likelihoods.len() {
//             //     if new_normalized_log_likelihoods[i] > new_mean + 4.0 * new_std_dev {
//             //         removed_outliers.push(i);
//             //         num_removed += 1;
//             //     }
//             // }
//         }
    
//         if num_removed == 0 {
//             break; 
//         }
//         println!("{} outliers removed", num_removed);
//         num_removed = 0;
//     }

//     for i in removed_outliers.iter() {
//         let mut output = cell_index_to_outputs.get_mut(i).unwrap();
//         println!("{}\t{}\t{}", output.barcode, output.likelihood,output.ground_truth);
//     }
    


// }




// fn init_alphas_betas_pairs(cell_data: &Vec<CellData>, loci_used: &HashSet<usize>, locus_to_cell_data: &HashMap<usize, HashMap<usize, [u32; 2]>>) -> HashMap<usize, Vec<usize>> {

//     // let mut alphas_betas_pairs : Vec<Vec<usize>> = vec![vec![0; loci_used.len()]; 3];
//     let mut alphas_betas_pairs : HashMap<usize, Vec<usize>> = HashMap::new();

//     for locus in loci_used.iter() {
        

//         let mut alpha = 1;
//         let mut beta = 1;

//         for (cell, counts) in locus_to_cell_data.get(locus).unwrap() {
//             alpha += counts[1];
//             beta += counts[0];
//         }

//         alphas_betas_pairs.insert(*locus, vec![alpha as usize, beta as usize]); // alt = alpha, ref = beta



//     }

//     alphas_betas_pairs
// }




// fn log_beta_binomial_PMF(locus_index: &usize, alphas_betas_pairs: HashMap<usize, Vec<usize>>, cell: &CellData, i: usize, index_to_locus: Vec<usize>) -> f64 {

//     let mut log_likelihood = 0.0 as f64;

//     let locus = index_to_locus[*locus_index];
//     let alpha = alphas_betas_pairs.get(&locus).unwrap()[0];
//     let beta = alphas_betas_pairs.get(&locus).unwrap()[1];

//     let log_binomial_coefficient = cell.log_binomial_coefficient[i];
//     let ref_count = cell.ref_counts[i];
//     let alt_count = cell.alt_counts[i];
    
//     let num_params = alt_count as usize + alpha ;
//     let denom_params =  ref_count as usize + beta ;


//     let num = log_beta_calc(num_params, denom_params);
//     let denom = log_beta_calc(alpha, beta);

//     log_likelihood = (log_binomial_coefficient as f64) + num - denom;


//     log_likelihood
// }
