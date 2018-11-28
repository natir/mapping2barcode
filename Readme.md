
Generate read to position, barcode and premolecule with ema:
```bash
zcat 1MB_1M.fq.gz | paste - - - - - - - - | grep "BX:Z:" | tr '\t' '\n' | ema align -t8 -r barcoded_minia.contigs.fa -1 /dev/stdin | grep -v "^@" | cut -d$'\t' -f 1,3,4,13,15 > read_barcode_premolecule.tsv
```
