use flate2::read::MultiGzDecoder;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::ffi::OsStr;
use File;
use izip;
use stats;
use std::process::Command;
use hashbrown::{HashMap,HashSet};
use Params;
use AlleleCount;

pub struct CellLocusData {
    pub locus_index: usize,
    pub locus_id: usize,
    pub log_binomial_coefficient: f64,
    pub alt_count: f64,
    pub ref_count: f64,
    pub total: usize,
}

pub struct CellData {
    pub cell_id: usize,
    pub barcode: String,
    pub assignment: String,
    pub cell_loci_data: Vec<CellLocusData>,
}

pub struct VcfLocusData {
    pub locus_index: usize,
    pub chrom: String,
    pub pos: String,
    pub ref_allele: String,
    pub alt_allele: String,
}

pub fn load_vcf_data(params: &Params) -> Option<Vec<VcfLocusData>> {
    if params.vcf.is_none() {
        return None;
    }
    let vcf = &params.vcf.as_ref().unwrap().to_string();
    let mut to_return: Vec<VcfLocusData> = Vec::new();
    let mut record_index: usize = 0;
    let reader = reader(&vcf);
    for line in reader.lines() {
        let line = line.expect("Unable to read a line in the ground truth file");
        if line.starts_with("#") { continue; }
        let toks: Vec<&str> = line.split('\t').collect();
        let chrom = toks[0].to_string();
        let pos = toks[1].to_string();
        let ref_allele = toks[3].to_string();
        let alt_allele = toks[4].to_string();
        to_return.push(VcfLocusData{
            locus_index: record_index,
            chrom: chrom,
            pos: pos,
            ref_allele: ref_allele,
            alt_allele: alt_allele,
        });
        record_index += 1;
    }
    return Some(to_return);
}


pub fn create_output_dir(params: &Params) {
    Command::new("mkdir")
        .arg(&params.output_directory)
        .output()
        .expect("failed to create output directory");
}

pub fn load_barcodes(params: &Params) -> (Vec<String>, HashMap<String, usize>) {
    let reader = reader(&params.barcodes);
    let mut cell_barcodes: Vec<String> = Vec::new();
    let mut barcode_cell_id: HashMap<String, usize> = HashMap::new();
    for (index, line) in reader.lines().enumerate() {
        let line = line.expect("Unable to read line");
        cell_barcodes.push(line.to_string());
        barcode_cell_id.insert(line.to_string(), index);
    }
    return (cell_barcodes, barcode_cell_id);
}

pub fn load_ground_truth(params: &Params, barcode_to_cell_id: &HashMap<String, usize>) ->
    Vec<String> {
    // vector of assignment (indexed by cell_id)

    let mut cell_id_to_ground_truth: Vec<String> = Vec::new();
    for _i in 0..barcode_to_cell_id.len() {
        cell_id_to_ground_truth.push("na".to_string());
    }
    if let Some(ground_truth) = &params.ground_truth {
        let reader = reader(&ground_truth);
        for line_result in reader.lines() {
            let line = line_result.expect("Unable to read a line in the ground truth file");
            let columns: Vec<&str> = line.split('\t').collect();
            assert!(columns.len() == 2, "Invalid line format: {}\nThe correct format is: barcode\tassignment", line);
            let barcode = columns[0].to_string();
            let assignment = columns[1].to_string();
            if let Some(cell_id) = barcode_to_cell_id.get(&barcode) {
                cell_id_to_ground_truth[*cell_id] = assignment.clone();
            }
        }
    }
    return cell_id_to_ground_truth;
}

pub fn load_mtx_final(params: &Params, excluded_cells: &HashSet<usize>) -> (Vec<AlleleCount>, Vec<AlleleCount>) {
    let mut locus_alleles_minority: Vec<AlleleCount> = Vec::new();
    let mut locus_alleles_majority: Vec<AlleleCount> = Vec::new();
    let (mut alt_reader, mut ref_reader) = (reader(&params.alt_mtx), reader(&params.ref_mtx));

    let (total_loci, _total_cells) = consume_mtx_header(&mut alt_reader, &mut ref_reader);
    for _ in 0..total_loci {
        locus_alleles_minority.push(AlleleCount{ alt_count: 0, ref_count: 0 });
        locus_alleles_majority.push(AlleleCount{ alt_count: 0, ref_count: 0 });
    }
    for (alt_line, ref_line) in izip!(alt_reader.lines(), ref_reader.lines()) {
        let (alt_line, ref_line) = (alt_line.expect("cannot read alt mtx"), ref_line.expect("cannot read ref mtx"));
        let data: VartrixDatum = read_mtx_lines(alt_line, ref_line);
        if excluded_cells.contains(&data.cell_id) {
            locus_alleles_minority[data.locus].alt_count += data.alt_count;
            locus_alleles_minority[data.locus].ref_count += data.ref_count;
        } else {
            locus_alleles_majority[data.locus].alt_count += data.alt_count;
            locus_alleles_majority[data.locus].ref_count += data.ref_count;
        }
    }

    return (locus_alleles_minority, locus_alleles_majority);
}

pub fn load_cell_data(params: &Params, cell_id_to_barcode: &Vec<String>, cell_id_to_assignment: &Vec<String>) ->
    (Vec<bool>, Vec<usize>, Vec<CellData>, Vec<[f64; 2]>, Vec<Vec<f64>>) { 
    // loci_used, vec of celldata, locus_counts (vec indexed by locus of [refcount, altcount])
    // 2 pass on mtx file. First to get loci_used then to get the cell data only for loci_used
    let (num_loci_used, loci_used, locus_to_used_index) = get_loci_used(params);
    let (mut alt_reader, mut ref_reader) = (reader(&params.alt_mtx), reader(&params.ref_mtx));
    let (_total_loci, total_cells) = consume_mtx_header(&mut alt_reader, &mut ref_reader);
    let mut cell_data = init_cell_data(total_cells, cell_id_to_barcode, cell_id_to_assignment);
    let mut locus_counts: Vec<[f64; 2]> = Vec::new();
    for _i in 0..num_loci_used { locus_counts.push([0.0;2]); }
    let mut locus_ids: Vec<usize> = Vec::new();
    for locus_id in 0..loci_used.len() {
        if loci_used[locus_id] { locus_ids.push(locus_id); }
    }
    // precompute some log_binomial_coefficients
    let max_n = 100;
    let precomputed_log_binomial_coefficients: Vec<Vec<f64>> = stats::precompute_log_binomial_coefficients(max_n);
    for (alt_line, ref_line) in izip!(alt_reader.lines(), ref_reader.lines()) {
        let (alt_line, ref_line) = (alt_line.expect("cannot read alt mtx"), ref_line.expect("cannot read ref mtx"));
        let data: VartrixDatum = read_mtx_lines(alt_line, ref_line);
        if !loci_used[data.locus] { continue; }
        let total = data.ref_count + data.alt_count;
        let used_locus_index = locus_to_used_index[data.locus];
        locus_counts[used_locus_index][0] += data.ref_count as f64;
        locus_counts[used_locus_index][1] += data.alt_count as f64;
        let log_coefficient: f64;
        if data.ref_count + data.alt_count <= max_n {
            log_coefficient = precomputed_log_binomial_coefficients[data.ref_count + data.alt_count][data.alt_count];
        } else {
            log_coefficient = statrs::function::factorial::ln_binomial((data.alt_count + data.ref_count) as u64, data.alt_count as u64) as f64;
        }
        cell_data[data.cell_id].cell_loci_data.push(
            CellLocusData{
                locus_index: used_locus_index,
                locus_id: data.locus,
                alt_count: data.alt_count as f64,
                ref_count: data.ref_count as f64,
                total: total,
                log_binomial_coefficient: log_coefficient,
            });
    }
    // now that we are only using used_loci, all loci are used (until later when we may filter loci bc they are problematic)
    let mut loci_used: Vec<bool> = Vec::new();
    for _locus in 0..num_loci_used {
        loci_used.push(true);
    }
    return (loci_used, locus_ids, cell_data, locus_counts, precomputed_log_binomial_coefficients);
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
    let locus = alt_tokens[0].to_string().parse::<usize>().unwrap() - 1; // 0-indexed Locus_ID
    let ref_count = ref_tokens[2].to_string().parse::<usize>().unwrap();
    let alt_count = alt_tokens[2].to_string().parse::<usize>().unwrap();
    let cell_id = alt_tokens[1].to_string().parse::<usize>().unwrap() - 1;
    return VartrixDatum {
        locus: locus,
        cell_id: cell_id,
        ref_count: ref_count,
        alt_count: alt_count,
    };
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

fn init_cell_data(total_cells: usize, cell_id_to_barcode: &Vec<String>, cell_id_to_assignment: &Vec<String>) ->
    Vec<CellData> {
    let mut cell_data: Vec<CellData> = Vec::new();
    for cell_id in 0..total_cells {
        cell_data.push(CellData{
            cell_id: cell_id,
            barcode: cell_id_to_barcode[cell_id].clone(),
            assignment: cell_id_to_assignment[cell_id].clone(),
            cell_loci_data: Vec::new(),
        });
    }
    return cell_data;
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


fn get_loci_used(params: &Params) -> (usize, Vec<bool>, Vec<usize>) {
    let (mut alt_reader, mut ref_reader) = (reader(&params.alt_mtx), reader(&params.ref_mtx));
    let (total_loci, _total_cells) = consume_mtx_header(&mut alt_reader, &mut ref_reader);
    let mut locus_counts: Vec<[usize; 2]> = Vec::new();
    let mut loci_used: Vec<bool> = Vec::new();
    let mut locus_to_used_index: Vec<usize> = Vec::new();
    for _i in 0..total_loci {
        locus_counts.push([0;2]);
        loci_used.push(false);
        locus_to_used_index.push(usize::MAX);
    }
    for (alt_line, ref_line) in izip!(alt_reader.lines(), ref_reader.lines()) {
        let (alt_line, ref_line) = (alt_line.expect("cannot read alt mtx"), ref_line.expect("cannot read ref mtx"));
        let data: VartrixDatum = read_mtx_lines(alt_line, ref_line);
        if data.ref_count > 0 { locus_counts[data.locus][0] += 1; }
        if data.alt_count > 0 { locus_counts[data.locus][1] += 1; }
    }
    let mut num_loci_used = 0;
    for locus in 0..total_loci {
        if locus_counts[locus][0] >= params.min_ref && locus_counts[locus][1] >= params.min_alt {
            loci_used[locus] = true;
            locus_to_used_index[locus] = num_loci_used;
            num_loci_used += 1;
        }
    }
    return (num_loci_used, loci_used, locus_to_used_index);
}
