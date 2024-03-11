#[macro_use]


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








struct Params {
    ref_mtx: String,
    alt_mtx: String,
    barcodes: String,
    min_alt: u32,
    min_ref: u32,
    ground_truth: String,
}