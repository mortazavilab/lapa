import pytest
import pyBigWig
import numpy as np
import pandas as pd
from lapa import lapa, read_polyA_cluster, read_apa_sample
from lapa.lapa import count_tes, TesCluster, tes_cluster_annotate, tes_sample
from lapa.utils.io import read_talon_read_annot
from conftest import fasta, gtf, read_annot, chrom_sizes, nbr1_hepg2_read_annot


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
        'Chromosome': ['chr17', 'chr17', 'chr17', 'chr17', 'chr17'],
        'End': [100100, 100101, 100101, 100101, 100200],
        'Strand': ['+', '+', '+', '-', '-'],
        'gene_id': ['gene_a', 'gene_a', 'gene_b', 'gene_c', 'gene_c'],
        'count': [10, 5, 11, 3, 7]
    })


def test_TesCluster_cluster(df_tes):
    clusters = list(TesCluster(fasta).cluster(df_tes))
    assert len(clusters) == 2

    cluster = clusters[0]
    assert cluster.Chromosome == 'chr1'
    assert cluster.Start == 100099
    assert cluster.End == 100101
    assert cluster.Strand == '+'
    assert cluster.counts == [10, 5]

    cluster = clusters[1]
    assert cluster.Chromosome == 'chr1'
    assert cluster.Start == 100100
    assert cluster.End == 100101
    assert cluster.Strand == '+'
    assert cluster.counts == [11]


def test_TesCluster_to_df(df_tes):
    df_clusters = TesCluster(fasta).to_df(df_tes)

    pd.testing.assert_frame_equal(
        df_clusters,
        pd.DataFrame({
            'Chromosome': ['chr17', 'chr17'],
            'Start': [100099, 100100],
            'End': [100101, 100101],
            'polyA_site': [100100, 100101],
            'count': [15, 11],
            'Strand': ['+', '+'],
            'fracA': [6, 6],
            'singal': ['100157@GATAAA', '100157@GATAAA']
        })
    )


def test_tes_cluster_annotation():

    df_cluster = pd.DataFrame({
        'Chromosome': ['chr17', 'chr17'],
        'Start': [100099, 100100],
        'End': [100101, 100101],
        'polyA_site': [100100, 100101],
        'count': [15, 11],
        'Strand': ['+', '+'],
        'fracA': [6, 6],
        'singal': ['100157@GATAAA', '100157@GATAAA']
    })

    df_annotate = tes_cluster_annotate(df_cluster, gtf)
    df_expected = pd.DataFrame({
        'Chromosome': ['chr17', 'chr17'],
        'Start': [100099, 100100],
        'End': [100101, 100101],
        'Strand': ['+', '+'],
        'polyA_site': [100100, 100101],
        'count': [15, 11],
        'fracA': [6, 6],
        'singal': ['100157@GATAAA', '100157@GATAAA'],
        'Feature': ['intergenic', 'intergenic'],
        'canonical_site': [-1, -1],
        'canonical': [False, False],
        'tpm': [15 * 1000000 / (15 + 11), 11 * 1000000 / (15 + 11)]
    })
    pd.testing.assert_frame_equal(df_annotate, df_expected,
                                  check_dtype=False, check_categorical=False)


def test_tes_sample():
    df_cluster = pd.DataFrame({
        'Chromosome': ['chr17', 'chr17'],
        'Start': [100099, 100100],
        'End': [100101, 100101],
        'Strand': ['+', '+'],
        'polyA_site': [100100, 100101],
        'count': [15, 11],
        'fracA': [6, 6],
        'singal': ['100157@GATAAA', '100157@GATAAA'],
        'Feature': ['intergenic', 'intergenic'],
        'canonical_site': [-1, -1],
        'canonical': [False, False],
        'tpm': [15 * 1000000 / (15 + 11), 11 * 1000000 / (15 + 11)]
    })
    df_apa = tes_sample(df_cluster, df_tes)

    import pdb
    pdb.set_trace()


def test_lapa(tmp_path):
    output_dir = tmp_path / 'lapa'

    lapa(read_annot, fasta, gtf, chrom_sizes, output_dir)

    df_cluster = read_polyA_cluster(str(output_dir / 'polyA_clusters.bed'))

    # TODO: write better test cases!
    assert all(df_cluster.columns == [
        'Chromosome', 'Start', 'End', 'polyA_site', 'count', 'Strand',
        'fracA', 'singal', 'Feature', 'canonical_site', 'canonical', 'tpm'
    ])

    assert set(df_cluster['Chromosome']) == {'chr17', 'ERCC-00060'}

    df_apa = read_apa_sample(str(output_dir / 'gm12878_apa.bed'))
    assert all(df_apa == [
        'Chromosome', 'Start', 'End', 'count', 'polyA_site', 'Strand',
        'gene_id', 'gene_count', 'usage', 'count_cluster',
        'fracA', 'singal', 'Feature', 'canonical_site', 'canonical', 'tpm'
    ])
