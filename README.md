# LAPA (Alpha release)

![tests](https://github.com/mortazavilab/lapa/actions/workflows/python-app.yml/badge.svg)
[![codecov](https://codecov.io/gh/mortazavilab/lapa/branch/master/graph/badge.svg?token=MJQ88T8JWK)](https://codecov.io/gh/mortazavilab/lapa)
[![Documentation Status](https://readthedocs.org/projects/lapa/badge/?version=latest)](https://lapa.readthedocs.io/en/latest/?badge=latest)

Alternative polyadenylation detection from diverse data sources such as 3'-seq, long-read and short-reads.

![method](docs/method.png)

## Installation

```
pip install lapa
```

## Run

```
lapa --alignment {bam} --fasta {fasta} --annotation {gtf} --chrom_sizes {chrom_sizes} --output_dir {output}
```

Argument details:

```
lapa --help

Usage: lapa [OPTIONS]

Options:
  --alignment TEXT    Bam or Sam files
  --annotation TEXT   Standart transcriptome Annotation (Encode or Ensembl
                      gtf)
  --fasta TEXT        Genome reference (Encode or Ensembl fasta)
  --chrom_sizes TEXT  Chrom sizes files (can be generated with `faidx fasta -i
                      chromsizes > chrom_sizes`)
  --output_dir TEXT   Output directory
  --help              Show this message and exit.

```

Chromsize file can be generated with following comment:

```
faidx {fasta} -i chromsizes > {chrom_sizes}
```

For multiple sample run, csv file annotating samples replaces bam file.

`samples.csv`:

```
sample,path
sample1_rep1,bam_sample1_rep1_path.bam
sample1_rep2,bam_sample1_rep2_path.bam
sample2_rep1,bam_sample2_rep1_path.bam
sample2_rep2,bam_sample2_rep2_path.bam
```

Then:

```
lapa --alignment samples.csv --fasta {fasta} --annotation {gtf} --chrom_sizes {chrom_sizes} --output_dir {output}
```
