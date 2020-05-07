
/* project mod */
mod parse_info;
//mod premolecule;

/* crates use */
use structopt::StructOpt;

/* std use */
use std::io::Write;
use std::time;
use std::collections::HashSet;
use std::collections::HashMap;

#[derive(Debug, StructOpt)]
#[structopt(name = "mapping2barcode", about = "Use mapping of barcode 10x read to assembly to build a barcode graph", author = "Pierre Marijon <pmarijon@mpi-inf.mpg.de>")]
struct Command {
    #[structopt(short = "e", long = "ema_info", help = "Summary of ema mapping result in tsv: read_id  contig  mapping_position  barcode_id  premolecule_id")]
    ema: String,

    #[structopt(short = "a", long = "asm", help = "assembly in fasta format")]
    asm: String,

    #[structopt(short = "o", long = "output", help = "path where barcode graph is write")]
    output: String,

    #[structopt(short = "l", long = "overlap-length", help = "minimum overlap length", default_value = "7000")]
    threshold: u64,

    #[structopt(short = "p", long = "premolecule-threshold", help = "maximal distance between two read ", default_value = "3500")]
    premolecule: u64,
}


fn main() {
    let params = Command::from_args();
    
    /* Read ema information */
    eprintln!("read ema info\n\tbegin");
    let mut begin = time::Instant::now();

    let tig2barcode2premol2pos = parse_info::ema(params.ema, params.premolecule, params.threshold);

    let mut duration = time::Instant::now() - begin;
    eprintln!("\tend {}s{}", duration.as_secs(), duration.subsec_millis());

    
    /* Read contig graph information */
    eprintln!("read assembly\n\tbegin");
    begin = time::Instant::now();

    let tig2len = parse_info::assembly(params.asm);

    duration = time::Instant::now() - begin;
    eprintln!("\tend {}s{}", duration.as_secs(), duration.subsec_millis());


    eprintln!("found edge of barcode graph\n\tbegin");
    begin = time::Instant::now();
    
    let mut nodes = HashSet::new();
    let mut edges = HashMap::new();
    for (tig, value) in tig2barcode2premol2pos {
	if let Some(len) = tig2len.get(&tig) {
	    if len < &(params.threshold as usize) {
		continue;
	    }
	}

	for (barcode1, poss1) in value.iter() {
	    for (barcode2, poss2) in value.iter() {
	    
		if barcode1 == barcode2 {
		    continue;
		}

		for pos1 in poss1 {
		    for pos2 in poss2 {
			if let Some(ovl_len) = get_ovl(*pos1, *pos2) {
			    if ovl_len > params.threshold {
				nodes.insert(barcode1.clone());
				nodes.insert(barcode2.clone());
				if barcode1 < barcode2 {
				    let mut val = edges.entry((barcode1.clone(), barcode2.clone())).or_insert(ovl_len);
				    if *val < ovl_len {
					*val = ovl_len;
				    }
				} else {
				    let mut val = edges.entry((barcode2.clone(), barcode1.clone())).or_insert(ovl_len);

				    if *val < ovl_len {
					*val = ovl_len;
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
    duration = time::Instant::now() - begin;
    eprintln!("\tend {}s{}", duration.as_secs(), duration.subsec_millis());

    eprintln!("write barcode graph\n\tbegin");
    begin = time::Instant::now();

    let mut writer = std::io::BufWriter::new(
	std::fs::File::create(params.output).expect("error opening output file")
    );
    
    writeln!(writer, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>").expect("error durring gexf write");
    writeln!(writer, "<gexf xmlns=\"http://www.gexf.net/1.2draft\" version=\"1.2\">").expect("error durring gexf write");
    writeln!(writer, "<graph mode=\"static\" defaultedgetype=\"undirected\">").expect("error durring gexf write");

    writeln!(writer, "<nodes>").expect("error durring gexf write");
    for node in nodes {
	writeln!(writer, "<node id=\"{}\" label=\"{}\" />", node, node).expect("error durring gexf write");
    }
    writeln!(writer, "</nodes>").expect("error durring gexf write");

    writeln!(writer, "<edges>").expect("error durring gexf write");
    let mut id = 0;
    for (edge, val) in edges {
	writeln!(writer, "<edge id=\"{}\" source=\"{}\" target=\"{}\" weight=\"{}\" />", id, edge.0, edge.1, val).expect("error durring gexf write");
	id += 1;
    }
    
    writeln!(writer, "</edges>").expect("error durring gexf write");
    
    writeln!(writer, "</graph>").expect("error durring gexf write");
    writeln!(writer, "</gexf>").expect("error durring gexf write");
    
    duration = time::Instant::now() - begin;
    eprintln!("\tend {}s{}", duration.as_secs(), duration.subsec_millis());
}

pub fn get_ovl(pos1: (u64, u64), pos2: (u64, u64)) -> Option<u64> {
    // pos1 contains pos2
    if pos2.0 > pos1.0 && pos2.1 < pos1.1 {
	return Some(pos2.1 - pos2.0)
    }

    // pos2 contains pos1
    if pos1.0 > pos2.0 && pos1.1 < pos2.1 {
	return Some(pos1.1 - pos1.0)
    }

    if pos2.1 > pos1.0 && pos2.1 < pos1.1 {
	return Some(pos2.1 - pos1.0)
    }

    if pos1.1 > pos2.0 && pos1.1 < pos2.1 {
	return Some(pos1.1 - pos2.0)
    }

    return None
}

