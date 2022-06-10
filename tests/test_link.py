import pytest
import numpy as np
from lapa.link import link_tss_to_tes
import pandas as pd
from conftest import read_annot, pb_brca1_bam


@pytest.fixture
def read_annot_link_path(tmp_path):
    read_annot_path = tmp_path / 'test_read_annot.tsv'

    pd.DataFrame({
        'read_name': ['r1', 'r2', 'r3'],
        'chrom': ['chr1', 'chr1', 'chr1'],
        'read_start': [5000, 5050, 90000],
        'read_end': [10001, 9970, 100000],
        'strand': ['+', '+', '+'],
        'annot_transcript_id': ['t1', 't1', 't2'],
        'annot_gene_id': ['g1', 'g1', 'g2'],
        'dataset': ['test', 'test', 'test']
    }).to_csv(read_annot_path, sep='\t',
              index=False)

    return read_annot_path


@pytest.fixture
def lapa_dir_link(tmp_path):
    lapa_dir = tmp_path / 'lapa'
    lapa_dir.mkdir()

    pd.DataFrame({
        'Chromosome': ['chr1', 'chr1'],
        'Start': [10000, 10000000],
        'End': [10010, 10000001],
        'polyA_site': [10005, 10000000],
        'count': [100, 100],
        'Strand': ['+', '+'],
        'Feature': ['three_prime_utr', 'three_prime_utr'],
        'gene_id': ['gene_a', 'gene_a'],
        'tpm': [1, 1],
        'gene_count': [100, 100],
        'usage': [0.5, 0.5],
        'fracA': [0, 0],
        'signal': ['AATAAA', 'AATAAA'],
        'annotated_site': [10005, 10000000]
    }).to_csv(lapa_dir / 'polyA_clusters.bed',
              index=False, sep='\t', header=False)

    return lapa_dir


@pytest.fixture
def lapa_tss_dir_link(tmp_path):
    lapa_tss_dir = tmp_path / 'lapa_tss'
    lapa_tss_dir.mkdir()

    pd.DataFrame({
        'Chromosome': ['chr1', 'chr1'],
        'Start': [4990, 100],
        'End': [5010, 110],
        'tss_site': [5000, 105],
        'count': [100, 100],
        'Strand': ['+', '+'],
        'Feature': ['five_prime_utr', 'five_prime_utr'],
        'gene_id': ['gene_a', 'gene_a'],
        'tpm': [1, 1],
        'gene_count': [2, 2],
        'usage': [0.5, 0.5],
        'annotated_site': [10005, 10000000]
    }).to_csv(lapa_tss_dir / 'tss_clusters.bed', index=False,
              sep='\t', header=False)

    return lapa_tss_dir


def test_link_reads_to_tes(read_annot_link_path,
                           lapa_dir_link,
                           lapa_tss_dir_link):

    df = link_tss_to_tes(read_annot_link_path,
                         lapa_dir_link,
                         lapa_tss_dir_link)

    pd.testing.assert_frame_equal(df, pd.DataFrame({
        'read_name': ['r1', 'r2', 'r3'],
        'Chromosome': pd.Categorical(['chr1', 'chr1', 'chr1']),
        'read_Start': [5000, 5050, 90000],
        'read_End': [10001, 9970, 100000],
        'Strand': pd.Categorical(['+', '+', '+']),
        'tss_site': [5000, 5000, -1],
        'polyA_site': [10005, 10005, -1]
    }))


def test_link_reads_to_tes_read_annot(lapa_read_annot,
                                      lapa_tss_read_annot):
    df = link_tss_to_tes(read_annot,
                         lapa_read_annot,
                         lapa_tss_read_annot,
                         distance=1000)

    _df = df.loc[(df['tss_site'] != -1) & (df['polyA_site'] != -1)]

    assert all(np.where(_df['Strand'] == '+',
                        _df['tss_site'] < _df['polyA_site'],
                        _df['polyA_site'] < _df['tss_site']))


def test_link_reads_to_tes_bam(lapa_pb_brca1,
                               lapa_tss_pb_brca1):
    df = link_tss_to_tes(pb_brca1_bam,
                         lapa_pb_brca1,
                         lapa_tss_pb_brca1,
                         distance=1000)

    _df = df.loc[(df['tss_site'] != -1) & (df['polyA_site'] != -1)]
    assert all(np.where(_df['Strand'] == '+',
                        _df['tss_site'] < _df['polyA_site'],
                        _df['polyA_site'] < _df['tss_site']))

    df2 = link_tss_to_tes(f'{pb_brca1_bam},{pb_brca1_bam}',
                          lapa_pb_brca1,
                          lapa_tss_pb_brca1,
                          distance=1000)

    _df2 = df2.loc[(df2['tss_site'] != -1) & (df2['polyA_site'] != -1)]
    assert all(np.where(_df2['Strand'] == '+',
                        _df2['tss_site'] < _df2['polyA_site'],
                        _df2['polyA_site'] < _df2['tss_site']))

    assert df.shape[0] * 2 == df2.shape[0]
