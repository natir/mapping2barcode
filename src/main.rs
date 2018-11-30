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
mod parse_info;
mod premolecule;

/* crates use */
use clap::{App, Arg};

/* std use */
use std::io::Write;

fn main() {
    let matches = App::new("mapping2barcodegraph")
        .version("0.2")
        .author("Pierre Marijon <pierre.marijon@inria.fr>")
        .about("Use mapping of barcode 10x read to assembly to build a barcode graph")
        .arg(Arg::with_name("ema")
             .short("e")
             .long("ema_info")
             .required(true)
             .display_order(10)
             .takes_value(true)
             .help("Summary of ema mapping result in tsv: read_id  contig  mapping_position  barcode_id  premolecule_id"))
        .arg(Arg::with_name("graph")
             .short("g")
             .long("graph")
             .required(true)
             .display_order(25)
             .takes_value(true)
             .help("Contig graph in gfa1 format")
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
        .get_matches();

    let ema_path = matches.value_of("ema").expect("Error durring ema info path access").to_string();
    let graph_path = matches.value_of("graph").expect("Error durring graph path access").to_string();
    let output_prefix = matches.value_of("output").expect("Error durring output prefix access").to_string();
    
    let threshold = matches.value_of("threshold").expect("Error durring threshold access").parse::<u64>().expect("Error durring threshold parsing");

    eprintln!("read ema info\n\tbegin");
    let (premolecule2tig_pos, barcode2premolecule, premolecule2reads, reads2barcode) = parse_info::ema(ema_path);
    eprintln!("\tend");

    for (p, t) in premolecule2tig_pos.iter() {
        println!("{:?} {:?}", p, t);
    }
    
    eprintln!("read assembly graph\n\tbegin");
    let (tig_graph, tig2len, tig2index) = parse_info::graph(graph_path);
    eprintln!("\tend");
    
    eprintln!("build premolecule graph\n\tbegin");
    let (premolecule_graph, _) = premolecule::build_graph(&tig_graph, &tig2len, &premolecule2tig_pos, &barcode2premolecule, &tig2index, threshold);
    eprintln!("\tend");
    
    eprintln!("write premolecule graph\n\tbegin");
    let mut graph_writer = std::io::BufWriter::new(std::fs::File::create(format!("{}_premolecule_graph.edges", output_prefix)).expect("Can't create graph file"));
       graph_writer.write(b"Source,Target,Weight\n").expect("Error durring premolecule graph header write");

    for e in premolecule_graph.raw_edges() {
        graph_writer.write_fmt(format_args!("{},{},{}\n", premolecule_graph[e.source()], premolecule_graph[e.target()], e.weight)).expect("Error durring premolecule graph write");
    }
    eprintln!("\tend");
    
    eprintln!("label reads with molecule\n\tbegin");
    let mut nb_molecule = 0;
    let mut assignation_writer = std::io::BufWriter::new(std::fs::File::create(format!("{}_read2molecule.tsv", output_prefix)).expect("Can't create result file"));
    for (id, cc) in petgraph::algo::kosaraju_scc(&premolecule_graph).iter().enumerate() {
        for node in cc {
            let premolecule = premolecule_graph.node_weight(*node).unwrap();
            for read in premolecule2reads.get(premolecule).unwrap() {
                let barcode = reads2barcode.get(read).unwrap();

                assignation_writer.write_fmt(format_args!("{}\t{}\t{}\n", barcode, read, id)).expect("Error durring read to molecule write");
            }
        }

        nb_molecule = id;
    }
    eprintln!("\tend");

    eprintln!("statistique:");
    eprintln!("\tnumber of premolecule\t{}", premolecule2reads.len());
    eprintln!("\tnumber of barcode\t{}", barcode2premolecule.len());
    eprintln!("\tnumber of molecule\t{}", nb_molecule + 1);
}

