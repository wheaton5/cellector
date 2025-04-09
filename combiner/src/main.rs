#[macro_use]
extern crate clap;
extern crate hashbrown;
extern crate rand;
extern crate flate2;
extern crate itertools;

use clap::App;
use rand::Rng;
use rand::seq::IteratorRandom;
use rand::rngs::StdRng;
use rand::SeedableRng;
use flate2::read::MultiGzDecoder;
use std::io::{BufWriter, Write, BufRead, BufReader};
use std::path::Path;
use std::ffi::OsStr;
use std::fs::File;
use std::process::Command;
use itertools::izip;

use hashbrown::{HashMap,HashSet};

fn main() {
    let params = load_params();
    create_output_dir(&params);
    let (locus_data2_to_data1, total_loci_out) = get_locus_mapping(&params);
    
    let (mut alt_reader, mut ref_reader) = (reader(&params.alt1), reader(&params.ref1));
    let (_total_loci, total_cells) = consume_mtx_header(&mut alt_reader, &mut ref_reader);
    let cells1 = select_cells(&params, params.num_cells_1, total_cells);
    let (mut alt_reader, mut ref_reader) = (reader(&params.alt2), reader(&params.ref2));
    let (_total_loci, total_cells) = consume_mtx_header(&mut alt_reader, &mut ref_reader);
    let cells2: Vec<usize>;
    let num_cells_2: usize;
    if params.dataset2_mask.is_some() {
        cells2 = select_cells_by_barcode(&params, params.dataset2_mask.clone().unwrap());
        num_cells_2 = cells2.len();
    } else {
        assert!(params.num_cells_2.is_some(), "missing argument num_cells_2 or dataset2_mask");
        cells2 = select_cells(&params, params.num_cells_2.unwrap(), total_cells);
        num_cells_2 = params.num_cells_2.unwrap();
    }
    // decide which cells are doublets and that mapping
    let (cell_ids_data1, cell_ids_data2) = get_cell_ids_and_output_barcodes(&params, &cells1, &cells2);
    output_new_mtxs(&params, &locus_data2_to_data1, total_loci_out, &cell_ids_data1, &cell_ids_data2, num_cells_2);
    // so now we need to make a new barcodes file
    println!("{},{}", params.num_cells_1, num_cells_2);
    // and now a new 
    //println!("{:?}", cells1);  
}

fn output_new_mtxs(params: &Params, locus_data2_to_data1: &HashMap<usize, usize>, total_loci_out: usize, cell_ids_data1: &HashMap<usize, usize>, cell_ids_data2: &HashMap<usize, usize>, num_cells_2: usize) {
    let filename = format!("{}/alt.mtx", params.output_directory);
    let filehandle = File::create(&filename).expect(&format!("Unable to create file {}", &filename));
    let mut alt_writer = BufWriter::new(filehandle);
    let filename = format!("{}/ref.mtx", params.output_directory);
    let filehandle = File::create(&filename).expect(&format!("Unable to create file {}", &filename));
    let mut ref_writer = BufWriter::new(filehandle);
    let tmpseed = params.seed.to_be_bytes();//[params.seed; 32]; // 4 guaranteed random number by fair dice roll https://xkcd.com/221/
    let mut seed = [0u8;32];
    for i in 0..8 { seed[i] = tmpseed[i]; }
    let mut rng: StdRng = SeedableRng::from_seed(seed);

    let total_cells = params.num_cells_1 + num_cells_2;
    let total_loci = total_loci_out;
    let total_entries = 0;
    alt_writer.write_all(b"%%MatrixMarket matrix coordinate real general\n% written by sprs\n");
    ref_writer.write_all(b"%%MatrixMarket matrix coordinate real general\n% written by sprs\n");
    alt_writer.write_all(format!("{}\t{}\t{}\n", total_loci, total_cells, total_entries).as_bytes());
    ref_writer.write_all(format!("{}\t{}\t{}\n", total_loci, total_cells, total_entries).as_bytes());

    let mut lines: Vec<(usize, usize, usize, usize)> = Vec::new(); // locus_id, cell_id, refcount, altcount

    let (mut alt_reader, mut ref_reader) = (reader(&params.alt1), reader(&params.ref1));
    let (total_loci, _total_cells) = consume_mtx_header(&mut alt_reader, &mut ref_reader);
    for (alt_line, ref_line) in izip!(alt_reader.lines(), ref_reader.lines()) {
        let (alt_line, ref_line) = (alt_line.expect("cannot read alt mtx"), ref_line.expect("cannot read ref mtx"));
        let mut data: VartrixDatum = read_mtx_lines(alt_line, ref_line);
        if let Some(&cell_id) = cell_ids_data1.get(&data.cell_id) {
            if data.locus < 1 {println!("what1? {}",data.locus); }
            let refcount = data.ref_count.clone();
            let altcount = data.alt_count.clone();
            for _ in 0..refcount {
                if rng.gen::<f64>() < params.downsample_rate { data.ref_count -= 1; }
            }
            for _ in 0..altcount {
                if rng.gen::<f64>() < params.downsample_rate { data.alt_count -= 1; }
            }
            lines.push((data.locus, cell_id, data.ref_count, data.alt_count));
        }
    }

    let (mut alt_reader, mut ref_reader) = (reader(&params.alt2), reader(&params.ref2));
    consume_mtx_header(&mut alt_reader, &mut ref_reader);
    for (alt_line, ref_line) in izip!(alt_reader.lines(), ref_reader.lines()) {
        let (alt_line, ref_line) = (alt_line.expect("cannot read alt mtx"), ref_line.expect("cannot read ref mtx"));
        let mut data: VartrixDatum = read_mtx_lines(alt_line, ref_line);
        if let Some(&cell_id) = cell_ids_data2.get(&data.cell_id) {
            if *locus_data2_to_data1.get(&data.locus).unwrap() < 1 { println!("what? {}", locus_data2_to_data1.get(&data.locus).unwrap());}
            let refcount = data.ref_count.clone();
            let altcount = data.alt_count.clone();
            for _ in 0..refcount {
                if rng.gen::<f64>() < params.downsample_rate { data.ref_count -= 1; }
            }
            for _ in 0..altcount {
                if rng.gen::<f64>() < params.downsample_rate { data.alt_count -= 1; }
            }
            lines.push((*locus_data2_to_data1.get(&data.locus).unwrap(), cell_id, data.ref_count, data.alt_count));
        }
    }
    lines.sort();
    for (locus_id, cell_id, ref_count, alt_count) in lines {
        alt_writer.write_all(format!("{}\t{}\t{}\n", locus_id, cell_id, alt_count).as_bytes());
        ref_writer.write_all(format!("{}\t{}\t{}\n", locus_id, cell_id, ref_count).as_bytes());
    }
}

struct VartrixDatum {
    locus: usize,
    cell_id: usize,
    alt_count: usize,
    ref_count: usize,
}

fn read_mtx_lines(alt_line: String, ref_line: String) ->
    VartrixDatum {
    let alt_tokens: Vec<&str> = alt_line.split_whitespace().collect();
    let ref_tokens: Vec<&str> = ref_line.split_whitespace().collect();
    let locus = alt_tokens[0].to_string().parse::<usize>().unwrap(); // 0-indexed Locus_ID
    let ref_count = ref_tokens[2].to_string().parse::<usize>().unwrap();
    let alt_count = alt_tokens[2].to_string().parse::<usize>().unwrap();
    let cell_id = alt_tokens[1].to_string().parse::<usize>().unwrap();
    return VartrixDatum {
        locus: locus,
        cell_id: cell_id,
        ref_count: ref_count,
        alt_count: alt_count,
    };
}

fn get_cell_ids_and_output_barcodes(params: &Params, cells1: &Vec<usize>, cells2: &Vec<usize>) -> (HashMap<usize, usize>, HashMap<usize, usize>) {
    let mut barcodes1: Vec<String> = Vec::new();
    let mut barcodes2: Vec<String> = Vec::new();
    let reader1 = reader(&params.barcodes1);
    for line in reader1.lines() {
        let line = line.expect("cannot read barcodes1");
        barcodes1.push(line.to_string());
    }
    let reader2 = reader(&params.barcodes2);
    for line in reader2.lines() {
        let line = line.expect("cannot read barcodes2");
        barcodes2.push(line.to_string());
    }
    let filename = format!("{}/barcodes.tsv",params.output_directory);
    let filehandle = File::create(&filename).expect(&format!("Unable to create file {}", &filename));
    let mut writer = BufWriter::new(filehandle);
    let filename = format!("{}/gt.tsv",params.output_directory);
    let filehandle = File::create(&filename).expect(&format!("Unable to create file {}", &filename));
    let mut gt_writer = BufWriter::new(filehandle);
    
    let mut cell_mapping_data1: HashMap<usize, usize> = HashMap::new();
    let mut cell_mapping_data2: HashMap<usize, usize> = HashMap::new();
    let mut cell_id_out = 1;
    for cell_id_in in cells1 {
        cell_mapping_data1.insert(*cell_id_in, cell_id_out);
        let mut bc = barcodes1[cell_id_in-1].to_string();
        bc.push_str("\n");
        writer.write_all(bc.as_bytes());
        cell_id_out += 1;
        let mut gt_line = barcodes1[cell_id_in-1].to_string();
        gt_line.push_str(&format!("\tmajority\n"));
        gt_writer.write_all(gt_line.as_bytes());
    }
    for cell_id_in in cells2 {
        cell_mapping_data2.insert(*cell_id_in, cell_id_out);
        let mut bc = barcodes2[cell_id_in-1].to_string();
        bc.pop();
        bc.push_str("2\n");
        let mut gt_line = bc.clone();
        gt_line.pop();
        writer.write_all(bc.as_bytes());
        cell_id_out += 1;
        gt_line.push_str(&format!("\tminority\n"));
        gt_writer.write_all(gt_line.as_bytes());
    }
    
    return (cell_mapping_data1, cell_mapping_data2);
}

pub fn create_output_dir(params: &Params) {
    Command::new("mkdir")
        .arg(&params.output_directory)
        .output()
        .expect("failed to create output directory");
}

fn get_locus_mapping(params: &Params) -> (HashMap<usize, usize>, usize) {
    let mut locus_data2_to_data1: HashMap<usize, usize> = HashMap::new(); // mapping from locus_id in dataset2 to locus_id in dataset1
    let mut chr_pos_to_locus_id: HashMap<(String, usize), usize> = HashMap::new();    
    let reader1 = reader(&params.vcf1);
    let mut record_number = 1; // because stupid matrix market format is 1 indexed, we start with 1 here
    for line in reader1.lines() {
        let line = line.expect("Unable to read a line in vcf1");
        if line.starts_with("#") { continue; }   
        let toks: Vec<&str> = line.split('\t').collect();
        let chrom = toks[0].to_string();
        let pos = toks[1].to_string();
        let pos = pos.parse::<usize>().unwrap();
        chr_pos_to_locus_id.insert((chrom, pos), record_number);
        record_number += 1;
    }
    let mut record_number2 = 1; // reset counter for new vcf
    let reader2 = reader(&params.vcf2);
    for line in reader2.lines() {
        let line = line.expect("Unable to read a line in vcf2");
        if line.starts_with("#") { continue; }
        let toks: Vec<&str> = line.split('\t').collect();
        let chrom = toks[0].to_string();
        let pos = toks[1].to_string();
        let pos = pos.parse::<usize>().unwrap();
        
        if let Some(locus_id) = chr_pos_to_locus_id.get(&(chrom, pos)) {
            locus_data2_to_data1.insert(record_number2, *locus_id);
        } else {
            locus_data2_to_data1.insert(record_number2, record_number);
            record_number += 1;
        }
        record_number2 += 1;
    }
    return (locus_data2_to_data1, record_number - 1);
}

pub fn reader(filename: &str) -> Box<dyn BufRead> {
    let path = Path::new(filename);
    let file = match File::open(&path) {
        Err(_why) => panic!("couldn't open file {}", filename),
        Ok(file) => file,
    };
    if path.extension() == Some(OsStr::new("gz")) {
        Box::new(BufReader::with_capacity(128 * 1024, MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::with_capacity(128 * 1024, file))
    }
}

fn select_cells(params: &Params, num_cells_to_use: usize, total_cells: usize) -> Vec<usize> {
    assert!(num_cells_to_use <= total_cells, "cant ask for more cells than exist in dataset");
    let tmpseed = params.seed.to_be_bytes();//[params.seed; 32]; // 4 guaranteed random number by fair dice roll https://xkcd.com/221/
    let mut seed = [0u8;32];
    for i in 0..8 { seed[i] = tmpseed[i]; }
    let mut rng: StdRng = SeedableRng::from_seed(seed);
    // do something
    let cell_ids = (1..(total_cells+1)).choose_multiple(&mut rng, num_cells_to_use);
    return cell_ids;
}

fn select_cells_by_barcode(params: &Params, dataset2_mask: String) -> Vec<usize> {
    let mut barcodes2: Vec<String> = Vec::new();
    let mut barcodes_mask: Vec<String> = Vec::new();
    let reader2 = reader(&params.barcodes2);
    for line in reader2.lines() {
        let line = line.expect("cannot read barcodes2");
        barcodes2.push(line.to_string());
    }
    let reader_mask = reader(&dataset2_mask);
    for line in reader_mask.lines() {
        let line = line.expect("cannot read dataset2_mask");
        barcodes_mask.push(line.to_string());
    }
    let mut cell_ids: Vec<usize> = Vec::new();
    for (id, barcode) in barcodes2.iter().enumerate() {
        if barcodes_mask.contains(barcode) {
            cell_ids.push(id+1);
        }
    }
    return cell_ids;
}

fn consume_mtx_header(alt_reader: &mut Box<dyn BufRead>, ref_reader: &mut Box<dyn BufRead>) ->
    (usize, usize) {
    // total_loci, total_cells
    let mut line = String::new();
    let mut total_loci = 0;
    let mut total_cells = 0;
    for x in 0..3 {
        alt_reader.read_line(&mut line).expect("cannot read line from alt matrix market file");
        line.clear();
        ref_reader.read_line(&mut line).expect("cannot read line from ref matrix market file");
        if x == 2 {
            let toks: Vec<&str> = line.split_whitespace().collect();
            total_loci  = toks[0].to_string().parse::<usize>().unwrap();
            total_cells = toks[1].to_string().parse::<usize>().unwrap();
        }
    }
    return (total_loci, total_cells);
}

pub struct Params {
    vcf1: String,
    vcf2: String,
    alt1: String,
    ref1: String,
    alt2: String,
    ref2: String,
    barcodes1: String,
    barcodes2: String,
    num_cells_1: usize,
    num_cells_2: Option<usize>,
    dataset2_mask: Option<String>,
    output_directory: String,
    seed: usize,
    downsample_rate: f64,
}

fn load_params() -> Params{
    let yaml = load_yaml!("params.yml");
    let params = App::from_yaml(yaml).get_matches();
    let vcf1 = params.value_of("vcf1").unwrap().to_string();
    let vcf2 = params.value_of("vcf2").unwrap().to_string();
    let alt1 = params.value_of("alt1").unwrap().to_string();
    let ref1 = params.value_of("ref1").unwrap().to_string();
    let alt2 = params.value_of("alt2").unwrap().to_string();
    let ref2 = params.value_of("ref2").unwrap().to_string();
    let barcodes1 = params.value_of("barcodes1").unwrap().to_string();
    let barcodes2 = params.value_of("barcodes2").unwrap().to_string();
    let num_cells_1 = params.value_of("num_cells_1").unwrap().to_string();
    let num_cells_1 = num_cells_1.parse::<usize>().unwrap();
    let num_cells_2 = match params.value_of("num_cells_2") {
        Some(x) => Some(x.to_string().parse::<usize>().unwrap()),
        None => None,
    };
    let dataset2_mask = match params.value_of("dataset2_mask") {
        Some(x) => Some(x.to_string()),
        None => None,
    };
    let output_directory = params.value_of("output_directory").unwrap().to_string();
    let seed = params.value_of("seed").unwrap_or("4").to_string();
    let seed = seed.parse::<usize>().unwrap();
    let downsample_rate = params.value_of("downsample_rate").unwrap_or("0.0").to_string();
    let downsample_rate = downsample_rate.parse::<f64>().unwrap();

    return Params {
        vcf1: vcf1,
        vcf2: vcf2,
        alt1: alt1,
        ref1: ref1,
        alt2: alt2,
        ref2: ref2,
        barcodes1: barcodes1,
        barcodes2: barcodes2,
        num_cells_1: num_cells_1,
        num_cells_2: num_cells_2,
        dataset2_mask: dataset2_mask,
        output_directory: output_directory,
        seed: seed,
        downsample_rate: downsample_rate,
    }
}
