# aa-table-with-af
Produce an 'amino acid table' file from a vcf, including depth and alt frequency info.

## Dependencies
This script uses [`pysam`](https://github.com/pysam-developers/pysam) to parse `.vcf` files.
pysam is available on both [pypi](https://pypi.org/project/pysam/) and [bioconda](https://anaconda.org/bioconda/pysam).

## Usage
```
usage: aa_table_with_af.py [-h] [-s SAMPLE_ID] vcf

positional arguments:
  vcf

optional arguments:
  -h, --help            show this help message and exit
  -s SAMPLE_ID, --sample_id SAMPLE_ID
```