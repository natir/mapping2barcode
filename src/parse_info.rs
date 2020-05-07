/* project use */
use file;

/* std use */
use std::collections::{HashMap, HashSet};

/* crates use */
use csv;

pub fn ema(tsv_path: String, premolecule_threshold: u64, ovl_threshold: u64) -> HashMap<String, HashMap<String, Vec<(u64, u64)>>> {

    let mut tig2barcode2poss: HashMap<String, HashMap<String, Vec<u64>>> = HashMap::new();

    let (reader, _) = niffler::from_path(&tsv_path).expect("ema file opening");

    let mut parser = csv::ReaderBuilder::new().delimiter(b'\t').flexible(true).has_headers(false).from_reader(reader);
    for result in parser.records() {
        let record = result.expect("Error during ema parsing");

        if record.len() != 5 {
            continue;
        } 
        
        let tig_id = record[1].to_string();
        let pos = record[2].parse::<u64>().unwrap();
        let barcode_id = record[3].to_string().split('-').next().unwrap().to_string();

        tig2barcode2poss.entry(tig_id.clone()).or_insert(HashMap::new()).entry(barcode_id.clone()).or_insert(Vec::new()).push(pos);
    }

    let mut tig2barcode2premol2pos: HashMap<String, HashMap<String, Vec<(u64, u64)>>> = HashMap::new();

    for (tig, value) in tig2barcode2poss {
	for (barcode, mut poss) in value {
	    let mut intervals = Vec::new();
	    
	    poss.sort();

	    if poss.len() < 2 {
		continue
	    }
	    
	    let mut iter = poss.iter();
	    let mut min = iter.next().unwrap();
	    let mut prev = min;
	    let mut last = min;
	
	    while let Some(next) = iter.next() {
		if next - prev > premolecule_threshold {
		    if *prev - *min > ovl_threshold {
			intervals.push((*min, *prev));
		    }
		    
		    min = next;
		}

		prev = next;
		last = next;
		
	    }
	    if *prev - *min > ovl_threshold {	
		intervals.push((*min, *prev));
	    }
	    
	    tig2barcode2premol2pos.entry(tig.clone()).or_insert(HashMap::new()).insert(barcode.clone(), intervals);
	}
    }
    
    return tig2barcode2premol2pos;
}


pub fn assembly(asm_path: String) -> HashMap<String, usize> {
    let mut tig2len: HashMap<String, usize> = HashMap::new();

    let (reader, _) = niffler::from_path(asm_path).expect("assembly file opening");
    let mut records = bio::io::fasta::Reader::new(reader).records();

    while let Some(Ok(record)) = records.next() {
	tig2len.insert(record.id().to_string(), record.seq().len());
    }

    tig2len
}
