import pytest
import pandas as pd
from lapa import lapa, lapa_tss
from lapa.link import link_tss_to_tes
from lapa.utils.io import read_talon_read_annot


gtf_brca1 = 'tests/data/brca1.gtf'

nbr1_hepg2_read_annot = 'tests/data/nbr1_hepg2_read_annot.tsv'
read_annot = 'tests/data/chr17_read_annot.tsv'
fasta = 'tests/data/chr17.fa'
chr17_chrom_sizes = 'tests/data/chr17.chrom_sizes'
gtf = 'tests/data/chr17.gtf'
chrom_sizes = 'tests/data/chrom_sizes'

quantseq_gm12_bam = 'tests/data/quantseq3_gm12878_chr17_rep1.bam'
quantseq_both_gm12_bam = 'tests/data/quantseq3_gm12878_chr17_rep1.bam,tests/data/quantseq3_gm12878_chr17_rep2.bam'

sample_csv = 'tests/data/alignment_sample.csv'

agg_annotation_gene = [
    'tests/data/agg_gene_0.csv',
    'tests/data/agg_gene_1.csv',
    'tests/data/agg_gene_brca1.csv'
]

short_bam = 'tests/data/short_chr17.bam'
gtf_gm12_pb = 'tests/data/chr17_gm12_pb_talon.gtf'
read_annot_gm12_pb = 'tests/data/chr17_gm12_read_annot.tsv'

pb_brca1_bam = 'tests/data/brc1_pb.bam'


@pytest.fixture
def lapa_read_annot(tmp_path):
    output_dir = tmp_path / 'lapa'
    lapa(read_annot, fasta, gtf, chrom_sizes, output_dir)
    return output_dir


@pytest.fixture
def lapa_tss_read_annot(tmp_path):
    output_dir = tmp_path / 'lapa_tss'
    lapa_tss(read_annot, fasta, gtf, chrom_sizes, output_dir)
    return output_dir


@pytest.fixture
def lapa_pb_brca1(tmp_path):
    output_dir = tmp_path / 'lapa'
    lapa(pb_brca1_bam, fasta, gtf, chrom_sizes, output_dir)
    return output_dir


@pytest.fixture
def lapa_tss_pb_brca1(tmp_path):
    output_dir = tmp_path / 'lapa_tss'
    lapa_tss(pb_brca1_bam, fasta, gtf, chrom_sizes, output_dir)
    return output_dir


@pytest.fixture
def lapa_pb_chr17(tmp_path):
    output_dir = tmp_path / 'lapa'
    lapa(read_annot_gm12_pb, fasta, gtf, chrom_sizes, output_dir)
    return output_dir


@pytest.fixture
def lapa_tss_pb_chr17(tmp_path):
    output_dir = tmp_path / 'lapa_tss'
    lapa_tss(read_annot_gm12_pb, fasta, gtf, chrom_sizes, output_dir)
    return output_dir


@pytest.fixture
def lapa_links_chr17(tmp_path, lapa_pb_chr17, lapa_tss_pb_chr17):
    output = tmp_path / 'lapa_links.csv'
    link_tss_to_tes(read_annot_gm12_pb, lapa_pb_chr17,
                    lapa_tss_pb_chr17).to_csv(output, index=False)
    return output
