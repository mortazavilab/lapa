from collections import Counter
import pytest
import pyBigWig
import numpy as np
import pandas as pd
from lapa import lapa, read_polyA_cluster, read_apa_sample
from lapa.lapa import count_tes, TesClustering, tes_cluster_annotate, \
    tes_sample
from lapa.utils.io import read_talon_read_annot, cluster_col_order, \
    sample_col_order
from conftest import fasta, gtf, read_annot, chrom_sizes, \
    nbr1_hepg2_read_annot, quantseq_gm12_bam, sample_csv


def test_count_tes(tmp_path):
    df = read_talon_read_annot(nbr1_hepg2_read_annot)
    df_tes, tes = count_tes(df, chrom_sizes, tmp_path)

    assert df.shape[0] == df_tes['count'].sum()
    assert len(tes.keys()) == 1
    assert df.shape[0] == tes['hepg2']['count'].sum()

    bw = pyBigWig.open(str(tmp_path / 'all_tes_counts_pos.bw'))
    count = np.nansum(
        bw.values('chr17',
                  df_tes['End'].min() - 1,
                  df_tes['End'].max())
    )
    assert df.shape[0] == count


@pytest.fixture
def df_tes():
    return pd.DataFrame({
        'Chromosome': ['chr17', 'chr17', 'chrM', 'chr17', 'chr17'],
        'Start': [100099, 100100, 1100, 100150, 100199],
        'End': [100100, 100101, 1101, 100151, 100200],
        'Strand': ['+', '+', '+', '-', '-'],
        'count': [10, 5, 11, 3, 7]
    })


def test_TesClustering_cluster(df_tes):
    clusters = list(TesClustering(fasta).cluster(df_tes))

    assert len(clusters) == 3

    cluster = clusters[0]
    assert cluster.Chromosome == 'chr17'
    assert cluster.Start == 100099
    assert cluster.End == 100101
    assert cluster.Strand == '+'
    assert cluster.counts == [10, 5]

    cluster = clusters[1]
    assert str(cluster) == 'chr17:100199-100200:-'
    assert cluster.counts == [7]

    cluster = clusters[2]
    assert cluster.Chromosome == 'chrM'
    assert cluster.Start == 1100
    assert cluster.End == 1101
    assert cluster.Strand == '+'
    assert cluster.counts == [11]


def test_TesCluster_to_df(df_tes):
    df_clusters = TesClustering(fasta).to_df(df_tes)

    df_expected = pd.DataFrame({
        'Chromosome': ['chr17', 'chr17', 'chrM'],
        'Start': [100099, 100199, 1100],
        'End': [100101, 100200, 1101],
        'polyA_site': [100100, 100200, 1101],
        'count': [15, 7, 11],
        'Strand': ['+', '-', '+'],
        'fracA': [6, 2, -1],
        'signal': ['100157@GATAAA', 'None@None', 'None@None']
    })
    pd.testing.assert_frame_equal(df_clusters, df_expected)


def test_tes_cluster_annotation():

    df_cluster = pd.DataFrame({
        'Chromosome': ['chr17', 'chrM', 'chr17', 'chr17'],
        'Start': [100099, 1100, 43144800, 43145100],
        'End': [100101, 1101, 43144850, 43145130],
        'polyA_site': [100100, 1101, 43144825, 43145114],
        'count': [15, 11, 5, 10],
        'Strand': ['+', '+', '+', '+'],
        'fracA': [6, 6, 1, 1],
        'signal': ['100157@GATAAA', '100157@GATAAA', 'None@None', 'None@None']
    })

    df_annotate = tes_cluster_annotate(df_cluster, gtf)

    df_expected = pd.DataFrame({
        'Chromosome': ['chr17', 'chr17', 'chr17', 'chr17', 'chrM'],
        'Start': [100099, 43144800, 43145100, 43145100, 1100],
        'End': [100101, 43144850, 43145130, 43145130, 1101],
        'polyA_site': [100100, 43144825, 43145114, 43145114, 1101],
        'count': [15, 5, 10, 10, 11],
        'Strand': ['+', '+', '+', '+', '+'],
        'fracA': [6, 1, 1, 1, 6],
        'signal': ['100157@GATAAA', 'None@None', 'None@None', 'None@None', '100157@GATAAA'],
        'Feature': ['intergenic', 'exon', 'exon', 'exon', 'intergenic'],
        'gene_id': ['', 'ENSG00000198496.12', 'ENSG00000198496.12', 'ENSG00000267681.1', ''],
        'gene_name': ['', 'NBR2', 'NBR2', 'CTD-3199J23.6', ''],
        'canonical_site': [-1, -1, -1, -1, -1],
        'tpm': [
            15 * 1000000 / (15 + 11 + 10 + 5),
            5 * 1000000 / (15 + 11 + 10 + 5),
            10 * 1000000 / (15 + 11 + 10 + 5),
            10 * 1000000 / (15 + 11 + 10 + 5),
            11 * 1000000 / (15 + 11 + 10 + 5)
        ]
    })

    df_annotate = df_annotate.reset_index()
    del df_annotate['index']

    pd.testing.assert_frame_equal(df_annotate, df_expected,
                                  check_dtype=False, check_categorical=False)


def test_tes_sample():

    df_tes = pd.DataFrame({
        'Chromosome': ['chr17', 'chr17'],
        'Start': [43144811, 43145111],
        'End': [43144812, 43145112],
        'Strand': ['+', '+'],
        'count': [10, 6]
    })

    df_cluster = pd.DataFrame({
        'Chromosome': ['chr17', 'chr17', 'chr17', 'chrM'],
        'Start': [100099, 43144800, 43145100, 1100],
        'End': [100101, 43144850, 43145130, 1101],
        'Strand': ['+', '+', '+', '+'],
        'polyA_site': [100100, 43144825, 43145114, 1101],
        'count': [15, 6, 10, 11],
        'fracA': [6, 1, 10, 1],
        'signal': ['100157@GATAAA', 'None@None', 'None@None', '100157@GATAAA'],
        'Feature': ['intergenic', 'exon', 'exon', 'intergenic'],
        'gene_id': ['', 'ENSG00000198496.12', 'ENSG00000198496.12', ''],
        'gene_name': ['', 'NBR2', 'NBR2', ''],
        'canonical_site': [-1, -1, 43145114, -1],
        'canonical': [False, False, True, False],
        'tpm': [
            15 * 1000000 / (15 + 11 + 10 + 6),
            10 * 1000000 / (15 + 11 + 10 + 6),
            11 * 1000000 / (15 + 11 + 10 + 6),
            6 * 1000000 / (15 + 11 + 10 + 6)
        ]
    })

    df_apa = tes_sample(df_cluster, df_tes, filter_internal_priming=False)

    assert df_apa.shape[0] == 2
    assert all(df_apa['count'] == [10, 6])
    assert all(df_apa['usage'] == [10 / 16, 6 / 16])

    df_apa = tes_sample(df_cluster, df_tes)

    assert df_apa.shape[0] == 1

    row = df_apa.iloc[0]
    assert row['Chromosome'] == 'chr17'
    assert row['Start'] == 43144800
    assert row['End'] == 43144850
    assert row['Strand'] == '+'


def test_lapa(tmp_path):
    output_dir = tmp_path / 'lapa'

    df_read_annot = read_talon_read_annot(read_annot)

    lapa(read_annot, fasta, gtf, chrom_sizes, output_dir)

    df_cluster = read_polyA_cluster(str(output_dir / 'polyA_clusters.bed'))

    assert all(df_cluster.columns == cluster_col_order)

    assert all(df_cluster['polyA_site'] <= df_cluster['End'])
    assert all(df_cluster['polyA_site'] >= df_cluster['Start'])

    assert (df_cluster['End'] - df_cluster['Start']).max() < 150

    assert set(df_cluster['Chromosome']) == {'chr17', 'ERCC-00060'}

    counts = df_cluster.drop_duplicates([
        'Chromosome', 'Start', 'End', 'Strand'])['count'].sum()
    assert df_read_annot.shape[0] > counts

    df_apa = read_apa_sample(str(output_dir / 'gm12878_apa.bed'))
    assert all(df_apa == sample_col_order)

    counts = pd.Series(Counter(df_apa['Feature']))
    assert counts.idxmax() == 'three_prime_utr'

    counts = df_apa.drop_duplicates([
        'Chromosome', 'Start', 'End', 'Strand'])['count'].sum()
    _df = df_read_annot[df_read_annot['sample'] == 'gm12878']
    assert _df.shape[0] > counts


def test_lapa_bam(tmp_path):
    output_dir = tmp_path / 'lapa'

    lapa(quantseq_gm12_bam, fasta, gtf, chrom_sizes, output_dir)

    df_cluster = read_polyA_cluster(str(output_dir / 'polyA_clusters.bed'))

    assert all(df_cluster.columns == cluster_col_order)

    assert set(df_cluster['Chromosome']) == {'chr17'}
    assert (df_cluster['End'] - df_cluster['Start']).max() < 150

    assert all(df_cluster['fracA'] <= 10)
    assert all(df_cluster['fracA'] >= 0)

    df_apa = read_apa_sample(
        str(output_dir / 'quantseq3_gm12878_chr17_rep1_apa.bed'))

    assert all(df_apa.columns == sample_col_order)

    counts = pd.Series(Counter(df_apa['Feature']))
    assert counts.idxmax() == 'three_prime_utr'


def test_lapa_read_csv(tmp_path):
    output_dir = tmp_path / 'lapa'

    lapa(sample_csv, fasta, gtf, chrom_sizes, output_dir)

    df_cluster = read_polyA_cluster(str(output_dir / 'polyA_clusters.bed'))

    assert set(df_cluster['Chromosome']) == {'chr17'}

    assert (df_cluster['End'] - df_cluster['Start']).max() < 150

    assert all(df_cluster['fracA'] <= 10)
    assert all(df_cluster['fracA'] >= 0)

    df_apa = read_apa_sample(str(output_dir / 'gm12878_apa.bed'))
    counts = pd.Series(Counter(df_apa['Feature']))

    assert all(df_apa.columns == sample_col_order)

    counts = pd.Series(Counter(df_apa['Feature']))
    assert counts.idxmax() == 'three_prime_utr'
