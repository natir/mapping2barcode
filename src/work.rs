/* project use */
use file;

/* crates use */
use csv;
use bio::io::{fastq, fasta};
use std::collections::{HashMap, HashSet};
use itertools::Itertools;

pub fn parse_graph(graph_path: String) -> (petgraph::Graph<String, String>, HashMap<String, u64>, HashMap<String, petgraph::graph::NodeIndex>) {
    let mut contig_len: HashMap<String, u64> = HashMap::new();
    let mut node2index: HashMap<String, petgraph::graph::NodeIndex> = HashMap::new();
    let mut contig_graph: petgraph::Graph<String, String> = petgraph::Graph::new();

    let (reader, _) = file::get_readable_file(&graph_path);
    
    let mut parser = csv::ReaderBuilder::new().delimiter(b'\t').flexible(true).from_reader(reader);
    for r in parser.records() {
        let record = r.expect("Error in contig gfa parssing");

        if &record[0] == "S" {
            contig_len.insert(record[1].to_string(), record[2].len() as u64);
        }
        
        if &record[0] == "L" {
            add_edge(&mut contig_graph, &mut node2index, record[1].to_string(), record[3].to_string(), format!("{}{}", record[2].to_string(), record[4].to_string()));
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

pub fn build_premolecule_graph(tig_graph: petgraph::Graph<String, String>, tig2len: HashMap<String, u64>, premolecule2tig: HashMap<String, (String, Vec<u64>)>, barcode2premolecule: HashMap<String, HashSet<String>>, tig2index: HashMap<String, petgraph::graph::NodeIndex>, threshold: u64) -> (petgraph::Graph<String, u64>, HashMap<String, petgraph::graph::NodeIndex>) {
    
    let mut premolecule_graph: petgraph::graph::Graph<String, u64> = petgraph::graph::Graph::new();
    let mut node2index: HashMap<String, petgraph::graph::NodeIndex> = HashMap::new();
    
    for (_, premolecule) in barcode2premolecule.iter() {
        for (p1, p2) in premolecule.iter().cartesian_product(premolecule.iter()) {
            if p1 == p2 {
                continue
            }

            let tig_ = &premolecule2tig.get(p1).expect("tig1_").0;
            let tig1 = *tig2index.get(tig_).expect("tig1");
            let tig2 = *tig2index.get(&premolecule2tig.get(p2).expect("tig2_").0).expect("tig2");
            let weight = if tig1 == tig2 {
                same_tig_dist(premolecule2tig.get(p1).expect("same tig p1"), premolecule2tig.get(p2).expect("same tig p2"))
            } else {
                other_tig_dist(tig1, premolecule2tig.get(p1).expect("other tig p1"), tig2, premolecule2tig.get(p2).expect("other tig p2"), &tig_graph, &tig2len, threshold)
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

fn same_tig_dist(p1: &(String, Vec<u64>), p2: &(String, Vec<u64>)) -> u64 {    
    let p1_begin = p1.1.iter().min().unwrap();
    let p2_begin = p2.1.iter().min().unwrap();


    return if p1_begin > p2_begin {
        p1_begin - p2_begin
    } else {
        p2_begin - p1_begin
    };
}

fn other_tig_dist(node1: petgraph::graph::NodeIndex, p1: &(String, Vec<u64>), node2: petgraph::graph::NodeIndex, p2: &(String, Vec<u64>), graph: &petgraph::graph::Graph<String, String>, tig2len: &HashMap<String, u64>, threshold: u64) -> u64 {

    let mut cumulative_len = 0;
    let mut relative_ori = "+";

    let p1_end = p1.1.iter().max().unwrap();
    let p2_begin = p2.1.iter().min().unwrap();    

    let (_, mut path) = petgraph::algo::astar(graph, node1, |finish| finish == node2, |_| 1, |_| 0).unwrap();

    let node_name1 = graph.node_weight(path.remove(0)).unwrap();
    cumulative_len += tig2len.get(node_name1).unwrap() - p1_end;

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
        "-" => cumulative_len += tig2len.get(node_name2).unwrap() - p2_begin,
        "+" | _ => cumulative_len += p2_begin,
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


pub fn get_ema_info(tsv_path: String) -> (HashMap<String, (String, Vec<u64>)>, HashMap<String, HashSet<String>>, HashMap<String, HashSet<String>>, HashMap<String, String>) {

    let mut premolecule2tig_pos : HashMap<String, (String, Vec<u64>)> = HashMap::new();
    let mut barcode2premolecule: HashMap<String, HashSet<String>> = HashMap::new();
    let mut premolecule2reads: HashMap<String, HashSet<String>> = HashMap::new();
    let mut reads2barcode: HashMap<String, String> = HashMap::new();
    
    let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').flexible(true).from_path(tsv_path).unwrap();
    for result in reader.records() {
        let record = result.expect("Error during ema parsing");

        if record.len() != 5 {
            continue;
        } 
        
        let read_id = record[0].to_string();
        let tig_id = record[1].to_string();
        let pos = record[2].parse::<u64>().unwrap();
        let barcode_id = record[3].to_string();
        let premolecule_id = format!("{}_{}", &record[3], &record[4]);
        
        premolecule2tig_pos.entry(premolecule_id.clone()).or_insert((tig_id.clone(), Vec::new())).1.push(pos);
        barcode2premolecule.entry(barcode_id.clone()).or_insert(HashSet::new()).insert(premolecule_id.clone());
        premolecule2reads.entry(premolecule_id.clone()).or_insert(HashSet::new()).insert(read_id.clone());
        reads2barcode.insert(read_id, barcode_id);
    }

    return (premolecule2tig_pos, barcode2premolecule, premolecule2reads, reads2barcode);
}
