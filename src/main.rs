extern crate bio;
extern crate xz2;
extern crate clap;
extern crate bzip2;
extern crate flate2;
extern crate itertools;
extern crate rust_htslib;

#[macro_use]
extern crate enum_primitive;

mod file;

use clap::{App, Arg};
use rust_htslib::bam;
use rust_htslib::bam::Read;
use bio::io::fastq;
use itertools::Itertools;

use std::ffi::CStr;
use std::io::Write;
use std::collections::{HashMap, HashSet};

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
        .arg(Arg::with_name("threshold")
             .short("t")
             .long("threshold")
             .display_order(40)
             .takes_value(true)
             .default_value("20")
             .help("Number of read map against contig to add barcode in clique")
        )
        .get_matches();

    let mut tig2reads: HashMap<String, HashSet<String>> = HashMap::new();

    let mut read_agains_tig = bam::Reader::from_path(matches.value_of("mapping").expect("Error durring mapping file access")).expec("Error durring mapping file reading");
    let header = read_agains_tig.header().clone();
    
    let mut nb_discard = 0;
    let mut nb_bad_mapq = 0;
    for r in read_agains_tig.records() {
        let record = r.expect("Trouble durring bam parsing");
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
        
        tig2reads.entry(ref_name.to_string()).or_insert(HashSet::new()).insert(que_name.to_string());
    }

    let mut read2barcode: HashMap<String, String> = HashMap::new();
    let (reader, _) = file::get_readable_file(matches.value_of("reads").expect("We have problem to determinate compression type"));
    let parser = fastq::Reader::new(reader);
    for r in parser.records() {
        let record = r.expect("Error durring fastq sequence parsing");

        read2barcode.insert(record.id().to_string(), record.desc().unwrap_or("NA").to_string());
    }

    let threshold = matches.value_of("threshold").expect("Error in threshold access").parse::<u32>().expect("Error in threshold parsing");
    let mut writer = std::fs::File::create(matches.value_of("output").expect("Error in output path access")).expect("Error durring output file creation");
    for (tig, reads) in tig2reads.iter() {
        let mut barcodes: HashMap<String, u32> = HashMap::new();

        for read in reads {
            if read2barcode.contains_key(read) {
                *barcodes.entry(read2barcode.get(read).expect("read isn't in barcode dict").to_string()).or_insert(0) += 1;
            }
        }

        let valid_barcodes = barcodes.into_iter().filter(|x| x.1 > threshold).map(|x| x.0).collect::<Vec<String>>();
        for (a, b) in valid_barcodes.iter().cartesian_product(valid_barcodes.iter()) {
            if a == b {
                continue;
            }
            writer.write_fmt(format_args!("{},{}\n", a, b));
        }
    }
    

}
