# primula
Primer database and design

## Description
This program integrates a simple SQLite primer database and design tool.

## Usage
To design primers and query existing
> `primula.py <VCF> --design`

add primers to database
> `primula.py -p <FASTA>`

retrieve (and design) primer for a single location
> `primula.py -q chr1:202030 --design`

## Reference Genome
Make sure to use a Repeatmasked reference genome. Even better if common SNPs are masked as well (eg. >1%).
BEDTOOLS offers `maskfasta` for this purpose. Masking data in BED format can be obtained from UCSC (http://genome.ucsc.edu/cgi-bin/hgTables?command=start).
