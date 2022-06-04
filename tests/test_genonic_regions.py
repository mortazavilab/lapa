import pytest
import numpy as np
import pandas as pd
import pyranges as pr
from lapa.genomic_regions import TssGenomicRegions, PolyAGenomicRegions
from conftest import gtf_brca1, agg_annotation_gene, gtf


@pytest.fixture
def polya_genomic_regions():
    return PolyAGenomicRegions(gtf_brca1)


@pytest.fixture
def tss_genomic_regions():
    return TssGenomicRegions(gtf)


def test_GenomicRegions_features(polya_genomic_regions):
    pr = polya_genomic_regions.features()

    assert set(pr.Feature) == {'five_prime_utr', 'exon',
                               'gene', 'three_prime_utr', 'intron'}


def test_GenomicRegions_annotate(polya_genomic_regions):
    gr = pr.PyRanges(chromosomes="chr17",
                     starts=(4303000, 43044826, 43046541, 43115728, 43093458),
                     ends=(4303100, 43045289, 43046997, 43115767, 43093573),
                     strands=('-', '-', '-', '-', '-'))

    gr.polyA_site = [4303050, 43045057, 43046769, 43057094, 43093515]
    df_ann = polya_genomic_regions.annotate(gr)

    df_ann['Chromosome'] = df_ann['Chromosome'].astype('str')
    df_ann['Strand'] = df_ann['Strand'].astype('str')

    df_expected = pd.DataFrame({
        'Chromosome': ['chr17', 'chr17', 'chr17', 'chr17', 'chr17'],
        'Start': [4303000, 43044826, 43046541, 43093458, 43115728],
        'End': [4303100, 43045289, 43046997, 43093573, 43115767],
        'Strand': ['-', '-', '-', '-', '-'],
        'polyA_site': [4303050, 43045057, 43046769, 43093515, 43057094],
        'Feature': ['intergenic', 'three_prime_utr', 'intron', 'three_prime_utr', 'exon'],
        'gene_id': ['intergenic_0', 'ENSG00000012048.23', 'ENSG00000012048.23',
                    'ENSG00000012048.23', 'ENSG00000012048.23'],
        'gene_name': ['intergenic_0', 'BRCA1', 'BRCA1', 'BRCA1', 'BRCA1'],
        'annotated_site': [-1, 43044294, -1, 43091434, -1]
    })

    pd.testing.assert_frame_equal(df_ann, df_expected)


def test_TssGenomicRegions_annotate(tss_genomic_regions):
    gr = pr.PyRanges(chromosomes="chr17",
                     starts=(50201626,),
                     ends=(50201686,),
                     strands=('-',))

    gr.tss_site = [50201631]
    df_ann = tss_genomic_regions.annotate(gr)

    df_ann['Chromosome'] = df_ann['Chromosome'].astype('str')
    df_ann['Strand'] = df_ann['Strand'].astype('str')

    df_expected = pd.DataFrame({
        'Chromosome': ['chr17'],
        'Start': [50201626],
        'End': [50201686],
        'Strand': ['-'],
        'tss_site': [50201631],
        'Feature': ['five_prime_utr'],
        'gene_id': ['ENSG00000108821.14'],
        'gene_name': ['COL1A1'],
        'annotated_site': [50201631]
    })

    pd.testing.assert_frame_equal(df_ann, df_expected)


def test_GenomicRegions__agg_annotation_gene(polya_genomic_regions):
    df = pd.read_csv(agg_annotation_gene[0])

    df_ann = polya_genomic_regions._agg_annotation_gene(
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
            'annotated_site': [3510502]
        })
    )

    df = pd.read_csv(agg_annotation_gene[1])

    df_ann = polya_genomic_regions._agg_annotation_gene(
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
            'annotated_site': [2043430, 2043430]
        })
    )
