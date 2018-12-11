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
use std::time;
use std::collections::HashMap;

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
             .default_value("10000")
             .help("Number of read map against contig to add barcode in clique")
        )
        .get_matches();

    let ema_path = matches.value_of("ema").expect("Error durring ema info path access").to_string();
    let graph_path = matches.value_of("graph").expect("Error durring graph path access").to_string();
    let output_prefix = matches.value_of("output").expect("Error durring output prefix access").to_string();
    
    let threshold = matches.value_of("threshold").expect("Error durring threshold access").parse::<u64>().expect("Error durring threshold parsing");

    /* Read ema information */
    eprintln!("read ema info\n\tbegin");
    let mut begin = time::Instant::now();

    let (premolecule2tig_pos, barcode2premolecule, premolecule2reads, reads2barcode) = parse_info::ema(ema_path);

    let mut duration = time::Instant::now() - begin;
    eprintln!("\tend {}s{}", duration.as_secs(), duration.subsec_millis());

    
    /* Read contig graph information */
    eprintln!("read assembly graph\n\tbegin");
    begin = time::Instant::now();

    let (tig_graph, tig2len, tig2index) = parse_info::graph(graph_path);

    duration = time::Instant::now() - begin;
    eprintln!("\tend {}s{}", duration.as_secs(), duration.subsec_millis());

    
    /* Build premolecule graph */
    eprintln!("build premolecule graph\n\tbegin");
    begin = time::Instant::now();

    let (premolecule_graph, premol_node2index) = premolecule::build_graph(&tig_graph, &tig2len, &premolecule2tig_pos, &barcode2premolecule, &tig2index, threshold);

    duration = time::Instant::now() - begin;
    eprintln!("\tend {}s{}", duration.as_secs(), duration.subsec_millis());

    
    /* Clean premolecule graph */
    eprintln!("clean premolecule graph\n\tbegin");
    begin = time::Instant::now();

    let clean_graph = premolecule::clean_graph(&premolecule_graph, &premol_node2index);
    premolecule::write_graph(&clean_graph, format!("{}_clean_graph.edges", output_prefix));
    
    duration = time::Instant::now() - begin;
    eprintln!("\tend {}s{}", duration.as_secs(), duration.subsec_millis());
    
    
    /* Write graph */
    eprintln!("write premolecule graph\n\tbegin");
    begin = time::Instant::now();

    premolecule::write_graph(&premolecule_graph, format!("{}_premolecule_graph.edges", output_prefix));

    duration = time::Instant::now() - begin;
    eprintln!("\tend {}s{}", duration.as_secs(), duration.subsec_millis());


    /* Labels reads */
    eprintln!("label reads with molecule\n\tbegin");
    begin = time::Instant::now();

    let mut nb_molecule = 0;
    let mut nb_pair_multi_assign = 0;
    let mut read_assign: HashMap<String, (String, String, usize)> = HashMap::new();

    let mut multi_assign_writer = std::io::BufWriter::new(std::fs::File::create(format!("{}_multi_assign.lst", output_prefix)).expect("Can't create multi assign read pairs file"));
    let mut read2premole_writer = std::io::BufWriter::new(std::fs::File::create(format!("{}_read2premolecule.lst", output_prefix)).expect("Can't create multi assign read pairs file"));    

    for (id, cc) in petgraph::algo::kosaraju_scc(&premolecule_graph).iter().enumerate() {
        for node in cc {
            let premolecule = premolecule_graph.node_weight(*node).unwrap();
            for read in premolecule2reads.get(premolecule).unwrap() {
                let barcode = reads2barcode.get(read).unwrap();
                let basic_read = &read[0..read.rfind("_").unwrap()];
                
                read2premole_writer.write_fmt(format_args!("{},{}\n", basic_read, premolecule)).unwrap();

                if read_assign.contains_key(basic_read) {
                    if id != read_assign.get(basic_read).expect("You can't be her").2 {
                        nb_pair_multi_assign += 1;
                        read_assign.remove(basic_read);
                        multi_assign_writer.write_fmt(format_args!("{}\n", basic_read)).expect("Error durring write in multi assignation file");
                    }
                    continue;
                } else {
                    read_assign.insert(basic_read.to_string(), (barcode.to_string(), basic_read.to_string(), id));
                }
            }
        }
        nb_molecule = id;
    }

    duration = time::Instant::now() - begin;
    eprintln!("\tend {}s{}", duration.as_secs(), duration.subsec_millis());


    /* Write result */
    eprintln!("write result file\n\tbegin");
    begin = time::Instant::now();
    
    let mut assignation_writer = std::io::BufWriter::new(std::fs::File::create(format!("{}_read2molecule.tsv", output_prefix)).expect("Can't create result file"));

    for (_, value) in read_assign.iter() {
        assignation_writer.write_fmt(format_args!("{}\t{}\t{}\n", value.0, value.1, value.2)).expect("Error durring read to molecule write");
    }
    
    duration = time::Instant::now() - begin;
    eprintln!("\tend {}s{}", duration.as_secs(), duration.subsec_millis());


    /* Write statistic */
    eprintln!("statistics:");
    eprintln!("\tnumber of reads mapped\t\t{}", reads2barcode.len());
    eprintln!("\tnumber of premolecule\t\t{}", premolecule2reads.len());
    eprintln!("\tnumber of barcode\t\t{}", barcode2premolecule.len());
    eprintln!("\tnumber of molecule\t\t{}", nb_molecule + 1); // molecule nb begin at 0
    eprintln!("\tnumber of reads pairs assigned\t{}", read_assign.len());
    eprintln!("\tnumber of pair multiple assign\t{}", nb_pair_multi_assign);
}

