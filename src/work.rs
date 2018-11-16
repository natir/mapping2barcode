/* project use */
use file;

/* crates use */
use bio::io::{fastq, fasta};
use std::collections::{HashMap, HashSet};
use rust_htslib::bam;
use rust_htslib::bam::Read;
use itertools::Itertools;


pub fn get_read_barcode_info(reads_path: String) -> (HashMap<String, String>, HashMap<String, HashSet<String>>) {
    
    let mut reads2barcode: HashMap<String, String> = HashMap::new();
    let mut barcode2reads: HashMap<String, HashSet<String>> = HashMap::new();
    
    let (reader, _) = file::get_readable_file(&reads_path);
    
    let parser = fastq::Reader::new(reader);
    for r in parser.records() {
        let record = r.expect("Error durring fastq sequence parsing");

        let read_id = record.id().to_string().split("/").next().expect("Error durring fastq header parssing").to_string();

        let raw_barcode = record.desc().unwrap_or("NA").to_string();
        let mut barcode = "NA".to_string();

        for b in raw_barcode.split(" ") {
            if b.starts_with("BX") {
                barcode = b.to_string();
                break;
            }
        }
        
        reads2barcode.insert(read_id.clone(), barcode.clone());
        barcode2reads.entry(barcode).or_insert(HashSet::new()).insert(read_id);
    }

    return (reads2barcode, barcode2reads);
}

pub fn get_read_mapping_info(mapping_path: String, min_mapq: u8) -> (HashMap<String, (String, u64, u64)>, HashMap<String, HashSet<(String, u64, u64)>>) {

    let mut reads2contig_pos : HashMap<String, (String, u64, u64)> = HashMap::new();
    let mut contigs2read_pos : HashMap<String, HashSet<(String, u64, u64)>> = HashMap::new();
    
    let mut read_agains_tig = bam::Reader::from_path(mapping_path).expect("Error durring mapping file reading");

    let header = read_agains_tig.header().clone();

    for r in read_agains_tig.records() {
        let record = r.expect("Trouble durring bam parsing");
        if record.is_secondary() || record.is_unmapped() || record.mapq() < min_mapq {
            continue
        }
        
        let read_name = String::from_utf8_lossy(record.qname()).to_string();
        let ctig_name = String::from_utf8_lossy(header.target_names()[record.tid() as usize]).to_string();
        let begin_pos = record.pos() as u64;
        let endin_pos = begin_pos + record.seq().len() as u64;

        reads2contig_pos.insert(read_name.clone(), (ctig_name.clone(), begin_pos, endin_pos));
        contigs2read_pos.entry(ctig_name).or_insert(HashSet::new()).insert((read_name, begin_pos, endin_pos));
    }

    return (reads2contig_pos, contigs2read_pos)
} 

pub fn build_premolecule(barcode2reads: HashMap<String, HashSet<String>>, reads2barcode: &HashMap<String, String>, reads2contig_pos: HashMap<String, (String, u64, u64)>, contigs2read_pos: HashMap<String, HashSet<(String, u64, u64)>>, max_mol_len: u64) -> (HashMap<String, HashSet<String>>, HashMap<String, HashSet<String>>, HashMap<String, Vec<String>>, HashMap<String, String>, HashMap<String, HashSet<(String, u64, u64)>>, HashMap<String, (String, u64, u64)>) {

    let mut premolecule2barcode : HashMap<String, HashSet<String>> = HashMap::new();
    let mut barcode2premolecule : HashMap<String, HashSet<String>> = HashMap::new();

    let mut premolecule2reads : HashMap<String, Vec<String>> = HashMap::new();
    let mut reads2premolecule : HashMap<String, String> = HashMap::new();

    let mut tig2premolecule_pos : HashMap<String, HashSet<(String, u64, u64)>> = HashMap::new();
    let mut premolecule2tig_pos : HashMap<String, (String, u64, u64)> = HashMap::new();

    for (tig, reads) in contigs2read_pos.iter() {
        let mut local_barcodes2reads: HashMap<String, Vec<(String, u64, u64)>> = HashMap::new();
        for read in reads {
              local_barcodes2reads.entry(reads2barcode.get(&read.0).unwrap_or(&"None".to_string()).to_string()).or_insert(Vec::new()).push(read.clone())
        }

        
        for (barcode, reads_pos) in local_barcodes2reads.iter_mut() {
            for (reads, begin, end) in split_read(reads_pos, max_mol_len) {
                let premolecule_id = vec![tig.clone(), barcode.clone(), begin.to_string(), end.to_string()].join("_");

                for read in reads.iter() {
                    reads2premolecule.insert(read.to_string(), premolecule_id.clone());
                }
                premolecule2tig_pos.insert(premolecule_id.clone(), (tig.clone(), begin, end));
                
                barcode2premolecule.entry(barcode.clone()).or_insert(HashSet::new()).insert(premolecule_id.clone());
                premolecule2barcode.entry(premolecule_id.clone()).or_insert(HashSet::new()).insert(barcode.clone());
                premolecule2reads.entry(premolecule_id.clone()).or_insert(Vec::new()).append(&mut reads.clone());
                tig2premolecule_pos.entry(tig.clone()).or_insert(HashSet::new()).insert((premolecule_id.clone(), begin, end));
            }
        }
        
    }

    return (premolecule2barcode, barcode2premolecule, premolecule2reads, reads2premolecule, tig2premolecule_pos, premolecule2tig_pos);
}

fn split_read(reads_pos: &mut Vec<(String, u64, u64)>, threshold: u64) -> HashSet<(Vec<String>, u64, u64)> {
    let mut cluster_of_reads: HashSet<(Vec<String>, u64, u64)> = HashSet::new();

    reads_pos.sort_by(|a, b| a.1.cmp(&b.1));

    let mut clusterd_read: Vec<String> = Vec::new();
    let mut first_read_beg_pos = reads_pos[0].1;
    let mut last_read_beg_pos = first_read_beg_pos;
    for (read, begin_pos, _) in reads_pos {
        if *begin_pos - last_read_beg_pos > threshold {
            cluster_of_reads.insert((clusterd_read, first_read_beg_pos, last_read_beg_pos));
            clusterd_read = Vec::new();
            first_read_beg_pos = *begin_pos;
        }
        clusterd_read.push(read.to_string());
        
        last_read_beg_pos = *begin_pos;
    }
    cluster_of_reads.insert((clusterd_read, first_read_beg_pos, last_read_beg_pos));

    return cluster_of_reads;
}


pub fn parse_graph(graph_path: String) -> (petgraph::Graph<String, String>, HashMap<String, u64>, HashMap<String, petgraph::graph::NodeIndex>) {
    let mut contig_len: HashMap<String, u64> = HashMap::new();
    let mut node2index: HashMap<String, petgraph::graph::NodeIndex> = HashMap::new();
    let mut contig_graph: petgraph::Graph<String, String> = petgraph::Graph::new();

    let (reader, _) = file::get_readable_file(&graph_path);
    
    let parser = fasta::Reader::new(reader);
    for r in parser.records() {
        let record = r.expect("Error in contig fasta parssing");
        let start = record.id();
        for desc in record.desc().unwrap_or("").split(" ") {
            if desc.starts_with("LN:i:") {
                contig_len.insert(start.to_string(), desc[5..].to_string().parse::<u64>().unwrap());
            }
            
            if desc.starts_with("L:") {
                let mut tmp = desc.split(":");
                tmp.next();
                let mut edge = String::new();
                edge.push(tmp.next().unwrap().chars().next().unwrap());
                let dest = tmp.next().unwrap();
                edge.push(tmp.next().unwrap().chars().next().unwrap());
                
                add_edge(&mut contig_graph, &mut node2index, start.to_string(), dest.to_string(), edge);
            }
        }
    }

    return (contig_graph, contig_len, node2index);
}

fn add_node(graph: &mut petgraph::Graph<String, String>, node: String, node2index: &mut HashMap<String, petgraph::graph::NodeIndex>) -> petgraph::graph::NodeIndex {
    return if node2index.contains_key(&node) {
        *node2index.get(&node).expect("You can't be her")
    } else {
        let index = graph.add_node(node.clone());
        node2index.insert(node, index);
        index
    };
}

fn add_edge(graph: &mut petgraph::Graph<String, String>, node2index: &mut HashMap<String, petgraph::graph::NodeIndex>, node_a: String, node_b: String, new_edge: String) {
  
    let n_a = add_node(graph, node_a, node2index);
    let n_b = add_node(graph, node_b, node2index);

    graph.add_edge(n_a, n_b, new_edge);
}

pub fn build_premolecule_graph(tig_graph: petgraph::Graph<String, String>, tig2len: HashMap<String, u64>, premolecule2tig: HashMap<String, (String, u64, u64)>, barcode2premolecule: HashMap<String, HashSet<String>>, tig2index: HashMap<String, petgraph::graph::NodeIndex>, threshold: u64) -> (petgraph::Graph<String, u64>, HashMap<String, petgraph::graph::NodeIndex>) {
    
    let mut premolecule_graph: petgraph::graph::Graph<String, u64> = petgraph::graph::Graph::new();
    let mut node2index: HashMap<String, petgraph::graph::NodeIndex> = HashMap::new();
    
    for (_, premolecule) in barcode2premolecule.iter() {
        for (p1, p2) in premolecule.iter().cartesian_product(premolecule.iter()) {
            if p1 == p2 {
                continue
            }
            
            let tig1 = *tig2index.get(&premolecule2tig.get(p1).unwrap().0).unwrap();
            let tig2 = *tig2index.get(&premolecule2tig.get(p2).unwrap().0).unwrap();
            let weight = if tig1 == tig2 {
                same_tig_dist(premolecule2tig.get(p1).unwrap(), premolecule2tig.get(p2).unwrap())
            } else {
                other_tig_dist(tig1, premolecule2tig.get(p1).unwrap(), tig2, premolecule2tig.get(p2).unwrap(), &tig_graph, &tig2len, threshold)
            };

            if weight == std::u64::MAX {
                continue;
            }
            add_edge_u64(&mut premolecule_graph, &mut node2index, p1.to_string(), p2.to_string(), weight);
        }
    }

    return (premolecule_graph, node2index);
}
fn add_node_u64(graph: &mut petgraph::Graph<String, u64>, node: String, node2index: &mut HashMap<String, petgraph::graph::NodeIndex>) -> petgraph::graph::NodeIndex {
    return if node2index.contains_key(&node) {
        *node2index.get(&node).expect("You can't be her")
    } else {
        let index = graph.add_node(node.clone());
        node2index.insert(node, index);
        index
    };
}

fn add_edge_u64(graph: &mut petgraph::Graph<String, u64>, node2index: &mut HashMap<String, petgraph::graph::NodeIndex>, node_a: String, node_b: String, new_edge: u64) {
  
    let n_a = add_node_u64(graph, node_a, node2index);
    let n_b = add_node_u64(graph, node_b, node2index);

    graph.add_edge(n_a, n_b, new_edge);
}

fn same_tig_dist(p1: &(String, u64, u64), p2: &(String, u64, u64)) -> u64 {
    return if p1.1 > p2.1 {
        p1.1 - p2.1
    } else {
        p2.1 - p1.1
    };
}

fn other_tig_dist(node1: petgraph::graph::NodeIndex, p1: &(String, u64, u64), node2: petgraph::graph::NodeIndex, p2: &(String, u64, u64), graph: &petgraph::graph::Graph<String, String>, tig2len: &HashMap<String, u64>, threshold: u64) -> u64 {

    let mut cumulative_len = 0;
    let mut relative_ori = "+";
    
    let (_, mut path) = petgraph::algo::astar(graph, node1, |finish| finish == node2, |_| 1, |_| 0).unwrap();

    let node_name1 = graph.node_weight(path.remove(0)).unwrap();
    cumulative_len += tig2len.get(node_name1).unwrap() - p1.2;

    let node_name2 = graph.node_weight(path.pop().unwrap()).unwrap();

    let mut previous_node = node1;
    for node in path {
        let tig = graph.node_weight(node).unwrap();
        cumulative_len += tig2len.get(tig).unwrap();

        if cumulative_len > threshold {
            return std::u64::MAX;
        }

        let edge = graph.edge_weight(graph.find_edge(previous_node, node).unwrap()).unwrap();
        
        if edge.as_str().get(0..0) != edge.as_str().get(1..1) {
           relative_ori = flip(relative_ori);
        }
        
        previous_node = node;
    }

    match relative_ori {
        "-" => cumulative_len += tig2len.get(node_name2).unwrap() - p2.1,
        "+" | _ => cumulative_len += p2.1,
    }
    
    return cumulative_len;
}

fn flip(a: &str) -> &str {
    match a {
        "+" => "-",
        "-" => "+",
        _ => "+",
    }
}
