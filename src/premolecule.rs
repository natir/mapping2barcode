/* std use */
use std::collections::{HashMap, HashSet};

/* crates use */
use itertools::Itertools;

/* std use */

pub fn build_graph(tig_graph: &petgraph::Graph<String, String>, tig2len: &HashMap<String, u64>, premolecule2tig: &HashMap<String, (String, Vec<u64>)>, barcode2premolecule: &HashMap<String, HashSet<String>>, tig2index: &HashMap<String, petgraph::graph::NodeIndex>, threshold: u64) -> (petgraph::Graph<String, u64>, HashMap<String, petgraph::graph::NodeIndex>) {
    
    let mut premolecule_graph: petgraph::graph::Graph<String, u64> = petgraph::graph::Graph::new();
    let mut node2index: HashMap<String, petgraph::graph::NodeIndex> = HashMap::new();

    let mut edges: HashSet<(String, String, u64)> = HashSet::new(); 
    for premolecules in barcode2premolecule.values() {
        edges.extend(compute_edge_weight(premolecules, premolecule2tig, tig2index, tig_graph, tig2len, threshold));
    }    

    for (p1, p2, weight) in edges.iter() {
        if weight == &std::u64::MAX {
            continue;
        }

        if *weight == 0 {
            add_edge_u64(&mut premolecule_graph, &mut node2index, p1.to_string(), p2.to_string(), 1);
        } else {
            add_edge_u64(&mut premolecule_graph, &mut node2index, p1.to_string(), p2.to_string(), *weight);
        }
    }

    return (premolecule_graph, node2index);
}

fn compute_edge_weight(premolecules: &HashSet<String>, premolecule2tig: &HashMap<String, (String, Vec<u64>)>, tig2index: &HashMap<String, petgraph::graph::NodeIndex>, tig_graph: &petgraph::Graph<String, String>, tig2len: &HashMap<String, u64>, threshold: u64) -> HashSet<(String, String, u64)> {
    let mut result: HashSet<(String, String, u64)> = HashSet::new();
    let mut path_buffer: HashMap<(petgraph::graph::NodeIndex, petgraph::graph::NodeIndex), Vec<petgraph::graph::NodeIndex>> = HashMap::new();
    
    for (p1, p2) in premolecules.iter().cartesian_product(premolecules.iter()) {
        if p1 == p2 {
            continue
        }

        let (tig1, pos1) = premolecule2tig.get(p1).expect("tig1");
        let tig_index1 = *tig2index.get(tig1).expect("tig_index1");
        let (tig2, pos2) = premolecule2tig.get(p2).expect("tig1");
        let tig_index2 = *tig2index.get(tig2).expect("tig_index2");

        let weight = get_edge_weight(&tig1, &tig2, pos1, pos2, &tig_index1, &tig_index2, &tig_graph, &tig2len, threshold, &mut path_buffer);

        result.insert((p1.to_string(), p2.to_string(), weight));
    }

    return result;
}

fn get_edge_weight(tig1: &String, tig2: &String, pos1: &Vec<u64>, pos2: &Vec<u64>, node1: &petgraph::graph::NodeIndex, node2: &petgraph::graph::NodeIndex, tig_graph: &petgraph::graph::Graph<String, String>, tig2len: &HashMap<String, u64>, threshold: u64, buffer: &mut HashMap<(petgraph::graph::NodeIndex, petgraph::graph::NodeIndex), Vec<petgraph::graph::NodeIndex>>) -> u64
{
    return if tig1 == tig2 {
        same_tig_dist(pos1, pos2)
    } else {
        other_tig_dist(node1, pos1, node2, pos2, &tig_graph, &tig2len, threshold, buffer)
    };
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

fn same_tig_dist(p1: &Vec<u64>, p2: &Vec<u64>) -> u64 {    
    let p1_begin = p1.iter().min().unwrap();
    let p2_begin = p2.iter().min().unwrap();

    return if p1_begin > p2_begin {
        p1_begin - p2_begin
    } else {
        p2_begin - p1_begin
    };
}

fn other_tig_dist(node1: &petgraph::graph::NodeIndex, p1: &Vec<u64>, node2: &petgraph::graph::NodeIndex, p2: &Vec<u64>, graph: &petgraph::graph::Graph<String, String>, tig2len: &HashMap<String, u64>, threshold: u64, buffer: &mut HashMap<(petgraph::graph::NodeIndex, petgraph::graph::NodeIndex), Vec<petgraph::graph::NodeIndex>>) -> u64 {
    
    let path;
    if buffer.contains_key(&(*node1, *node2)) {
        path = buffer.get(&(*node1, *node2)).unwrap().clone();
    } else {
        let p = petgraph::algo::astar(graph, *node1, |finish| finish == *node2, |_| 1, |_| 0);

        match p.is_none() {
            true => return std::u64::MAX,
            false => path = p.unwrap().1,
        };
        buffer.insert((*node1, *node2), path.clone());
    }
    
    if path.len() == 2 {
        return other_tig_dist_one_edge(p1, p2, graph, tig2len, threshold, &path);
    } else {
        return other_tig_dist_x_edge(p1, p2, graph, tig2len, threshold, &path)
    }        
}

fn other_tig_dist_one_edge(p1: &Vec<u64>, p2: &Vec<u64>, graph: &petgraph::graph::Graph<String, String>, tig2len: &HashMap<String, u64>, threshold: u64, path: &Vec<petgraph::graph::NodeIndex>) -> u64 {
    let mut cumulative_len = 0;

    let mut edges = path.windows(2);
    let node_pair = edges.next().unwrap();
    let (edge, weight) = get_edge(node_pair, graph);

    let node1 = graph.node_weight(graph.raw_edges()[edge.index()].source()).unwrap();
    let length1 = tig2len.get(node1).unwrap();
    
    let node2 = graph.node_weight(graph.raw_edges()[edge.index()].target()).unwrap();
    let length2 = tig2len.get(node2).unwrap();

    if weight.as_bytes()[0] == b'+' {
        cumulative_len += length1 - p1.iter().max().unwrap();
    } else {
        cumulative_len += p1.iter().min().unwrap();
    };

    if weight.as_bytes()[1] == b'+' {
        cumulative_len += p2.iter().min().unwrap();
    } else {
        cumulative_len += length2 - p2.iter().max().unwrap();
    };

    if cumulative_len > threshold {
        return std::u64::MAX;
    } else {
        return cumulative_len;
    }
}

fn other_tig_dist_x_edge(p1: &Vec<u64>, p2: &Vec<u64>, graph: &petgraph::graph::Graph<String, String>, tig2len: &HashMap<String, u64>, threshold: u64, path: &Vec<petgraph::graph::NodeIndex>) -> u64 {
    let mut cumulative_len = 0;
    
    let mut edges = path.windows(2);
    let mut node_pair = edges.next().unwrap();
    let (edge, weight) = get_edge(node_pair, graph);

    let node1 = graph.node_weight(graph.raw_edges()[edge.index()].source()).unwrap();
    let length1 = tig2len.get(node1).unwrap();
    
    if weight.as_bytes()[0] == b'+' {
        cumulative_len += length1 - p1.iter().max().unwrap();
    } else {
        cumulative_len += p1.iter().min().unwrap();
    };
    
    for e in edges {
        node_pair = e;
        let (edge, _) = get_edge(node_pair, graph);

        let node = graph.node_weight(graph.raw_edges()[edge.index()].source()).unwrap();
        let length = tig2len.get(node).unwrap();

        cumulative_len += length;

        if cumulative_len > threshold {
            return std::u64::MAX;
        }
    }

    let (edge, weight) = get_edge(node_pair, graph);
    
    let node2 = graph.node_weight(graph.raw_edges()[edge.index()].target()).unwrap();
    let length2 = tig2len.get(node2).unwrap();
    
    if weight.as_bytes()[1] == b'+' {
        cumulative_len += p2.iter().min().unwrap();
    } else {
        cumulative_len += length2 - p2.iter().max().unwrap();
    };

    if cumulative_len > threshold {
        return std::u64::MAX;
    } else {
        return cumulative_len;
    }
}

fn get_edge(node_weight_couple: &[petgraph::prelude::NodeIndex], graph: &petgraph::graph::Graph<String, String>) -> (petgraph::graph::EdgeIndex, String) {
    let edge = graph.find_edge(node_weight_couple[0], node_weight_couple[1]).unwrap();
    let weight = graph.edge_weight(edge).unwrap();

    return (edge, weight.to_string());
}
