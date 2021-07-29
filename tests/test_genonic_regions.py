import pytest
from conftest import gtf_brca1
import pyranges as pr
from longread_postprocessing.genomic_regions import GenomicRegions


@pytest.fixture
def genomic_regions():
    return GenomicRegions(gtf_brca1)


def test_genomic_regions_features(genomic_regions):
    pr = genomic_regions.features()

    assert set(pr.Feature) == {'five_prime_utr', 'exon',
                               'gene', 'three_prime_utr', 'intron'}


def test_genomics_regions_annotate(genomic_regions):
    gr = pr.PyRanges(chromosomes="17",
                     starts=(43124097, 43115726, 43104960, 43042754),
                     ends=(43124115, 43115779, 43104990, 43044056),
                     strands=('-', '-', '-', '-'))
    gr.gene = 'brca1'

    gr_ann = genomic_regions.annotate(gr)
    assert gr_ann.df.shape == (4, 6)

    gr_ann = genomic_regions.annotate(gr, single=True)
    assert gr_ann.df['Feature'].tolist() == [
        'five_prime_utr', 'five_prime_utr', 'intron', 'intergenic']
    assert gr_ann.df.shape == (4, 6)
