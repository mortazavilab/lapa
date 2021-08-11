import pytest
from conftest import gtf_brca1, agg_annotation_gene
import numpy as np
import pandas as pd
import pyranges as pr
from lapa.genomic_regions import GenomicRegions


@pytest.fixture
def genomic_regions():
    return GenomicRegions(gtf_brca1)


def test_GenomicRegions_features(genomic_regions):
    pr = genomic_regions.features()

    assert set(pr.Feature) == {'five_prime_utr', 'exon',
                               'gene', 'three_prime_utr', 'intron'}


def test_GenomicRegions_annotate(genomic_regions):
    gr = pr.PyRanges(chromosomes="17",
                     starts=(43124097, 43115726, 43104960, 43042754, 43044293),
                     ends=(43124115, 43115779, 43104990, 43044056, 43044300),
                     strands=('-', '-', '-', '-', '-'))
    gr.polyA_site = [43124106, 43115752, 43104975, 43043405, 43044296]
    df_ann = genomic_regions.annotate(gr)

    df_expected = pd.DataFrame({
        'Chromosome': ['17', '17', '17', '17', '17'],
        'Start': [43042754, 43044293, 43104960, 43115726, 43124097],
        'End': [43044056, 43044300, 43104990, 43115779, 43124115],
        'Strand': ['-', '-', '-', '-', '-'],
        'polyA_site': [43043405, 43044296, 43104975, 43115752, 43124106],
        'Feature': ['intergenic', 'three_prime_utr', 'intron', 'exon', 'exon'],
        'gene_id': ['', 'ENSG00000012048', 'ENSG00000012048', 'ENSG00000012048', 'ENSG00000012048'],
        'gene_name': ['', 'BRCA1', 'BRCA1', 'BRCA1', 'BRCA1'],
        'canonical_site': [-1, 43044294, -1, -1, -1]
    })
    pd.testing.assert_frame_equal(df_ann, df_expected)


def test_GenomicRegions__agg_annotation_gene():
    df = pd.read_csv(agg_annotation_gene[0])

    df_ann = GenomicRegions._agg_annotation_gene(
        df).reset_index(drop=True)
    df_ann.columns.name = None

    pd.testing.assert_frame_equal(
        df_ann,
        pd.DataFrame({
            'Chromosome': ['chr17'],
            'Start': [3510502],
            'End': [3510505],
            'polyA_site': [3510503],
            'count': [27],
            'Strand': ['-'],
            'fracA': [2],
            'signal': ['3510539@AATAAA'],
            'Feature': ['three_prime_utr'],
            'gene_id': ['ENSG00000167723.15'],
            'gene_name': ['TRPV3'],
            'canonical_site': [3510502]
        })
    )

    df = pd.read_csv(agg_annotation_gene[1])

    df_ann = GenomicRegions._agg_annotation_gene(
        df).reset_index(drop=True)
    df_ann.columns.name = None

    pd.testing.assert_frame_equal(
        df_ann,
        pd.DataFrame({
            'Chromosome': ['chr17', 'chr17'],
            'Start': [2043419, 2043419],
            'End': [2043430, 2043430],
            'polyA_site': [2043424, 2043424],
            'count': [761, 761],
            'Strand': ['+', '+'],
            'fracA': [3, 3],
            'signal': ['2043461@AATAAA', '2043461@AATAAA'],
            'Feature': ['three_prime_utr', 'three_prime_utr'],
            'gene_id': ['ENSG00000108963.19', 'ENSG00000262664.3'],
            'gene_name': ['DPH1', 'OVCA2'],
            'canonical_site': [2043430, 2043430]
        })
    )
