import pandas as pd
from lapa.cluster import TesClustering
from conftest import fasta


def test_tes_cluster():

    clustering = TesClustering(fasta, groupby='gene_id', fields=['reads'])

    df = pd.DataFrame({
        'Chromosome': ['chr1', 'chr1', 'chr1', 'chr1'],
        'Start': [10, 11, 50, 55],
        'End': [11, 12, 51, 56],
        'Strand': ['+', '+', '-', '-'],
        'count': [2, 3, 5, 3],
        'gene_id': ['gene_a', 'gene_a', 'gene_b', 'gene_b'],
        'reads': [
            ('r1', 'r2'),
            ('r3', '4', 'r5'),
            ('r4', 'r5', 'r6', 'r7', 'r8'),
            ('r9', 'r10', 'r11')
        ]
    })

    clusters = list(clustering.cluster(df))

    assert len(clusters) == 2
    assert str(clusters[0]) == 'chr1:11-12:+'
    assert str(clusters[1]) == 'chr1:50-56:-'

    assert clusters[0].fields == {'reads': [('r3', '4', 'r5')]}
    assert clusters[1].fields == {
        'reads': [
            ('r4', 'r5', 'r6', 'r7', 'r8'),
            ('r9', 'r10', 'r11')
        ]
    }
