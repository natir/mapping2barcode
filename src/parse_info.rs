/* project use */
use file;

/* std use */
use std::collections::{HashMap, HashSet};

/* crates use */
use csv;

pub fn graph(graph_path: String) -> (petgraph::Graph<String, String>, HashMap<String, u64>, HashMap<String, petgraph::graph::NodeIndex>) {
    let mut contig_len: HashMap<String, u64> = HashMap::new();
    let mut node2index: HashMap<String, petgraph::graph::NodeIndex> = HashMap::new();
    let mut contig_graph: petgraph::Graph<String, String> = petgraph::Graph::new();

    let (reader, _) = file::get_readable_file(&graph_path);
    
    let mut parser = csv::ReaderBuilder::new().delimiter(b'\t').flexible(true).from_reader(reader);
    for r in parser.records() {
        let record = r.expect("Error in contig gfa parssing");

        if &record[0] == "S" {
            contig_len.insert(record[1].to_string(), record[2].len() as u64);
            add_node(&mut contig_graph, record[1].to_string(), &mut node2index);
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

pub fn ema(tsv_path: String) -> (HashMap<String, (String, Vec<u64>)>, HashMap<String, HashSet<String>>, HashMap<String, HashSet<String>>, HashMap<String, String>) {

    let mut premolecule2tig_pos : HashMap<String, (String, Vec<u64>)> = HashMap::new();
    let mut barcode2premolecule: HashMap<String, HashSet<String>> = HashMap::new();
    let mut premolecule2reads: HashMap<String, HashSet<String>> = HashMap::new();
    let mut reads2barcode: HashMap<String, String> = HashMap::new();

    let (reader, _) = file::get_readable_file(&tsv_path);

    let mut parser = csv::ReaderBuilder::new().delimiter(b'\t').flexible(true).has_headers(false).from_reader(reader);
    for result in parser.records() {
        let record = result.expect("Error during ema parsing");

        if record.len() != 5 {
            continue;
        } 
        
        let read_id = record[0].to_string();
        let tig_id = record[1].to_string();
        let pos = record[2].parse::<u64>().unwrap();
        let barcode_id = record[3].to_string();
        let premolecule_id = format!("{}_{}", &record[3], &record[4]);

        let new_read_id = found_read_id(&reads2barcode, &read_id);
        
        premolecule2tig_pos.entry(premolecule_id.clone()).or_insert((tig_id.clone(), Vec::new())).1.push(pos);
        barcode2premolecule.entry(barcode_id.clone()).or_insert(HashSet::new()).insert(premolecule_id.clone());
        premolecule2reads.entry(premolecule_id.clone()).or_insert(HashSet::new()).insert(new_read_id.clone());
        reads2barcode.insert(new_read_id, barcode_id);
    }

    return (premolecule2tig_pos, barcode2premolecule, premolecule2reads, reads2barcode);
}

fn found_read_id(reads2barcode: &HashMap<String, String>, read: &String) -> String {
    let mut i: u32 = 1;
    let mut new_id = format!("{}_{}", read.to_string(), i);
    
    while reads2barcode.contains_key(&new_id) {
        i += 1;
        new_id = format!("{}_{}", read, i);
    }

    return new_id;
}
