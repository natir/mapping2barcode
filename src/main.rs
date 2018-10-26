extern crate bio;
extern crate clap;
extern crate rust_htslib;

use clap::{App, Arg};
use rust_htslib::bam;
use rust_htslib::bam::Read;
use bio::io::fastq;

use std::slice;
use std::ffi::CStr;
use std::collections::HashMap;

fn main() {
    let matches = App::new("mapping2barcodegraph")
        .version("0.1")
        .author("Pierre Marijon <pierre.marijon@inria.fr>")
        .about("Use mapping of barcode 10x read to assembly to build a barcode graph")
        .arg(Arg::with_name("reads")
             .short("r")
             .long("reads")
             .display_order(10)
             .takes_value(true)
             .help("Reads with this header format [>|@]read_id barcode"))
        .arg(Arg::with_name("mapping")
             .short("m")
             .long("mapping")
             .display_order(20)
             .takes_value(true)
             .help("Mapping of reads against assembly only bam file")
        )
        .arg(Arg::with_name("output")
             .short("o")
             .long("output")
             .display_order(30)
             .takes_value(true)
             .help("Where graph is write")
        )
        .get_matches();

    let mut tig2reads: HashMap<String, Vec<String>> = HashMap::new();

    let mut readsAtig = bam::Reader::from_path(matches.value_of("mapping").unwrap()).unwrap();
    let header = readsAtig.header().clone();
    
    let mut nb_discard = 0;
    let mut nb_bad_mapq = 0;
    for r in readsAtig.records() {
        let record = r.unwrap();
        if record.is_secondary() || record.is_unmapped() {
            nb_discard += 1;
            continue
        }
        
        if record.mapq() < 60 {
            nb_discard += 1;
            nb_bad_mapq += 1;
            continue
        }

        let ref_name = String::from_utf8_lossy(header.target_names()[record.tid() as usize]);
        let que_name = String::from_utf8_lossy(record.qname());
        
        println!("{:?} {:?}", ref_name, que_name);
        tig2reads.entry(ref_name.to_string()).or_insert(Vec::new()).push(que_name.to_string());
    }

    let mut read2barcode: HashMap<String, String> = HashMap::new();
    let reader = fastq::Reader::from_file(matches.value_of("reads").unwrap()).unwrap();
    for r in reader.records() {
        let record = r.unwrap();

        reader.insert(record.id(), record.desc().unwrap_or("NA"));
    }

    
    
    println!("Hello, world!");
}
