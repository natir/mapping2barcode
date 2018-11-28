extern crate bio;
extern crate xz2;
extern crate csv;
extern crate clap;
extern crate bzip2;
extern crate flate2;
extern crate petgraph;
extern crate itertools;

#[macro_use]
extern crate enum_primitive;

/* project mod */
mod file;
mod work;

/* crates use */
use clap::{App, Arg};

/* std use */
use std::collections::{HashMap, HashSet};

fn main() {
    let matches = App::new("mapping2barcodegraph")
        .version("0.2")
        .author("Pierre Marijon <pierre.marijon@inria.fr>")
        .about("Use mapping of barcode 10x read to assembly to build a barcode graph")
        .arg(Arg::with_name("reads")
             .short("r")
             .long("reads")
             .required(true)
             .display_order(10)
             .takes_value(true)
             .help("Reads with this header format [>|@]read_id barcode"))
        .arg(Arg::with_name("mapping")
             .short("m")
             .required(true)
             .long("mapping")
             .display_order(20)
             .takes_value(true)
             .help("Mapping of reads against assembly only bam file")
        )
        .arg(Arg::with_name("graph")
             .short("g")
             .required(true)
             .long("graph")
             .display_order(25)
             .takes_value(true)
             .help("Contig graph of reads in fasta format")
        )
        .arg(Arg::with_name("output")
             .short("o")
             .required(true)
             .long("output")
             .display_order(30)
             .takes_value(true)
             .help("Output prefix")
        )
        .arg(Arg::with_name("threshold")
             .short("t")
             .long("threshold")
             .display_order(40)
             .takes_value(true)
             .default_value("100000")
             .help("Number of read map against contig to add barcode in clique")
        )
        .arg(Arg::with_name("min_mapq")
             .short("M")
             .long("minimal-mapq")
             .display_order(50)
             .takes_value(true)
             .default_value("50")
             .help("If mapping quality is less than this threshold the mapping is discard")
        )
        .get_matches();

    let reads_path = matches.value_of("reads").expect("Error durring reads path access").to_string();
    let mapping_path = matches.value_of("mapping").expect("Error durring map path access").to_string();
    let graph_path = matches.value_of("graph").expect("Error durring graph path access").to_string();
    
    let threshold = matches.value_of("threshold").expect("Error durring threshold access").parse::<u64>().expect("Error durring threshold parsing");
    let min_mapq = matches.value_of("min_mapq").expect("Error durring minimal mapq access").parse::<u8>().expect("Error durring minimal mapq parsing");

    eprintln!("read ema info\n\tbegin");
    let (premolecule2tig_pos, barcode2premolecule, premolecule2reads, reads2barcode) = work::get_ema_info(reads_path);
    eprintln!("\tend");
    
    eprintln!("read assembly graph\n\tbegin");
    let (tig_graph, tig2len, tig2index) = work::parse_graph(graph_path);
    eprintln!("\tend");
    
    eprintln!("build pre molecule graph\n\tbegin");
    let (premolecule_graph, premolecule2index) = work::build_premolecule_graph(tig_graph, tig2len, premolecule2tig_pos, barcode2premolecule, tig2index, threshold);
    eprintln!("\tend");
    
    eprintln!("label read with molecule\n\tbegin");
    for (id, cc) in petgraph::algo::kosaraju_scc(&premolecule_graph).iter().enumerate() {
        for node in cc {
            let premolecule = premolecule_graph.node_weight(*node).unwrap();
            for read in premolecule2reads.get(premolecule).unwrap() {
                let barcode = reads2barcode.get(read).unwrap();

                println!("{}\t{}\t{}", barcode, read, id);
            }
        }
    }
    eprintln!("\tend");
}

