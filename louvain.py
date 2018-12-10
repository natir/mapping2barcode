#!/usr/bin/env python3

import csv
import community
import networkx as nx

from collections import defaultdict

def main(graph_path, read2premol_path):
    graph = nx.Graph()
    read2barcode = dict()
    premolecule2reads = defaultdict(set)
    
    with open(graph_path, "r") as g:
        reader = csv.DictReader(g, delimiter=',')
        for r in reader:
            graph.add_edge(r["Source"], r["Target"], weight=1/float(r["Weight"]))

    with open(read2premol_path) as premole:
        reader = csv.reader(premole, delimiter=',')
        for r in reader:
            read_id = r[0]
            barcode_id = r[1].split("_")[0]

            read2barcode[read_id] = barcode_id
            premolecule2reads[r[1]].add(read_id)
            
    partition = community.best_partition(graph)
    for premol, parti in partition.items():
        for read in premolecule2reads[premol]:
            print("\t".join((read2barcode[read], read, str(parti))))
        
if __name__ == "__main__":
    import sys
    main(sys.argv[1], sys.argv[2])
