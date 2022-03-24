import pytest
import pandas as pd
from lapa.cluster import PolyAClustering, Cluster
from conftest import fasta


@pytest.fixture
def df_tes():
    return pd.DataFrame({
        'Chromosome': ['chr17', 'chr17', 'chrM', 'chr17', 'chr17'],
        'Start': [100099, 100100, 1100, 100150, 100199],
        'End': [100100, 100101, 1101, 100151, 100200],
        'Strand': ['+', '+', '+', '-', '-'],
        'count': [10, 5, 11, 3, 7]
    })


def test_PolyAClustering_cluster(df_tes):
    clusters = list(PolyAClustering(fasta, extent_cutoff=5).cluster(df_tes))

    assert len(clusters) == 3

    cluster = clusters[0]
    assert cluster.Chromosome == 'chr17'
    assert cluster.Start == 100099
    assert cluster.End == 100101
    assert cluster.Strand == '+'
    assert cluster.counts == [(100100, 10), (100101, 5)]

    cluster = clusters[1]
    assert str(cluster) == 'chr17:100199-100200:-'
    assert cluster.counts == [(100200, 7)]

    cluster = clusters[2]
    assert cluster.Chromosome == 'chrM'
    assert cluster.Start == 1100
    assert cluster.End == 1101
    assert cluster.Strand == '+'
    assert cluster.counts == [(1101, 11)]


def test_TesCluster_to_df(df_tes):
    df_clusters = PolyAClustering(fasta, extent_cutoff=5).to_df(df_tes)

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

    cols = ['Chromosome', 'Start', 'End', 'polyA_site',
            'count', 'Strand', 'fracA', 'signal']
    pd.testing.assert_frame_equal(df_clusters[cols], df_expected)


def test_tes_cluster():

    clustering = PolyAClustering(fasta, groupby='gene_id', fields=['reads'])

    df = pd.DataFrame({
        'Chromosome': ['chr1', 'chr1', 'chr1', 'chr1'],
        'Start': [10, 11, 50, 55],
        'End': [11, 12, 51, 56],
        'Strand': ['+', '+', '-', '-'],
        'count': [2, 3, 3, 5],
        'gene_id': ['gene_a', 'gene_a', 'gene_b', 'gene_b'],
        'reads': [
            ('r1', 'r2'),
            ('r3', '4', 'r5'),
            ('r9', 'r10', 'r11'),
            ('r4', 'r5', 'r6', 'r7', 'r8')
        ]
    })

    clusters = list(clustering.cluster(df))

    assert len(clusters) == 2
    assert str(clusters[0]) == 'chr1:11-12:+'
    assert str(clusters[1]) == 'chr1:50-56:-'

    assert clusters[0].fields == {'reads': [('r3', '4', 'r5')]}
    assert clusters[1].fields == {
        'reads': [
            ('r9', 'r10', 'r11'),
            ('r4', 'r5', 'r6', 'r7', 'r8')
        ]
    }

    polya_site = clusters[1].polyA_site()
    assert polya_site == 56


def test_Cluster__count_arr():
    counts = [
        (5, 10),
        (6, 15),
        (7, 16),
        (8, 3),
        (18, 23),
        (19, 24),
        (21, 18)
    ]
    count_arr = Cluster._count_arr(counts)
    assert count_arr == [
        10, 15, 16, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 23, 24, 0, 18]
