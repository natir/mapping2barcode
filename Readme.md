
# Requirements

- [rust environement setup](https://rustup.rs/)
- [ema](http://cb.csail.mit.edu/cb/ema/)

# Install

```
cargo install --git https://github.com/natir/mapping2barcode.git
```

Or use source

```
git clone https://github.com/natir/mapping2barcode.git
cd mapping2barcode
cargo install
```

# Usage

```bash
samtools faidx {reference}
bwa index {reference}
zcat {reads}.gz | paste - - - - - - - - | grep "BX:Z:" | tr '\t' '\n' | ema align -t8 -r {reference} -1 /dev/stdin | grep -v "^@" | cut -d$'\t' -f 1,3,4,13,15 > {output}.tsv
mapping2barcodegraph -a {reference} -e {output}.tsv -o {output}.gexf -l 9000 -p 5000
```

