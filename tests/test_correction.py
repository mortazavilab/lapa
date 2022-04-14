import pytest
import pandas as pd
import pyranges as pr
from lapa.correction import correct_talon, \
    _links_transcript_agg, _transcript_tss_tes, \
    correct_talon, TranscriptModifier
from conftest import gtf_gm12_pb, read_annot_gm12_pb, \
    lapa_links_chr17, gtf_gm12_pb, chr17_abundance


@pytest.fixture
def read_annot_link_path(tmp_path):
    read_annot_path = tmp_path / 'test_read_annot.tsv'

    pd.DataFrame({
        'read_name': ['r1', 'r2', 'r3', 'r6', 'r7'],
        'chrom': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
        'read_start': [5000, 6000, 90000, 5000, 4000],
        'read_end': [10001, 9750, 100000, 9000, 10000],
        'strand': ['+', '+', '+', '-', '+'],
        'annot_transcript_id': ['t1', 't1', 't2', 't1', 't1'],
        'annot_gene_id': ['g1', 'g1', 'g2', 'g1', 'g1'],
        'dataset': ['test', 'test', 'test', 'test', 'test']
    }).to_csv(read_annot_path, sep='\t', index=False)

    return read_annot_path


@pytest.fixture
def link_path(tmp_path):
    link_path = tmp_path / 'test_links.csv'

    pd.DataFrame({
        'read_name': ['r1', 'r2', 'r6', 'r3', 'r4', 'r5', 'r7'],
        'Chromosome': pd.Categorical(['chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr5', 'chr1']),
        'read_Start': [5000, 6000, 5000, 90000, 10, 100, 4000],
        'read_End': [10001, 9750, 9000, 100000, 20, 200, 10000],
        'Strand': pd.Categorical(['+', '+', '+', '+', '+', '-', '+']),
        'start_site': [5000, 5000, 5000, -1, 10, 200, 4000],
        'polyA_site': [10005, 10005, -1, -1, 20, -1, 10000]
    }).to_csv(link_path, index=False)

    return link_path


def test__transcript_tss_tes(link_path, read_annot_link_path):
    df = _links_transcript_agg(link_path, read_annot_link_path)
    pd.testing.assert_frame_equal(df, pd.DataFrame({
        'transcript_id': ['t1', 't1', 't1', 't2'],
        'Strand': ['+', '+', '+', '+'],
        'start_site': [4000, 5000, 5000, -1],
        'polyA_site': [10000, -1, 10005, -1],
        'read_Start': [4000, 5000, 5000, 90000],
        'read_End': [10000, 9000, 10001, 100000],
        'count': [1, 1, 2, 1]
    }))

    df = _transcript_tss_tes(df)
    pd.testing.assert_frame_equal(df, pd.DataFrame({
        'transcript_id': ['t1_0', 't1_1'],
        'start_site': [4000, 5000],
        'polyA_site': [10000, 10005],
        'count': [1, 2]
    }))


@pytest.fixture
def modifier():
    return TranscriptModifier(gtf_gm12_pb)


def test_TranscriptModifier_fetch_transcript(modifier):
    transcript = modifier.fetch_transcript('ENST00000308278.12')

    assert transcript.transcript_id == 'ENST00000308278.12'
    assert transcript.df.Feature.tolist() == [
        'transcript', 'exon', 'exon', 'exon', 'exon', 'exon']

    assert transcript


def test_Transcript_copy(modifier):
    transcript = modifier.fetch_transcript('ENST00000308278.12') \
                         .copy('ENST00000308278.12_0')

    assert transcript.transcript_id == 'ENST00000308278.12_0'
    df_exons = transcript.df[1:]
    assert all(df_exons['exon_id'].str.endswith('_0'))


def test_Transcript_update_start_site(modifier):
    transcript = modifier.fetch_transcript('ENST00000308278.12') \
                         .copy('ENST00000308278.12_0')

    transcript.update_start_site(732411 - 100)
    transcript.df['Start'].iloc[0] == 732411 - 100
    transcript.df['Start'].iloc[1] == 732411 - 100

    transcript = modifier.fetch_transcript('ENST00000579788.5') \
                         .copy('ENST00000579788.5_0')

    transcript.update_start_site(64020602 + 100)
    transcript.df['End'].iloc[0] == 64020602 + 100
    transcript.df['End'].iloc[0] == 64020602 + 100


def test_Transcript_update_polya_site(modifier):
    transcript = modifier.fetch_transcript('ENST00000308278.12') \
                         .copy('ENST00000308278.12_0')

    transcript.update_polyA_site(742972 + 100)
    transcript.df['Start'].iloc[0] == 742972 + 100
    transcript.df['Start'].iloc[1] == 742972 + 100

    transcript = modifier.fetch_transcript('ENST00000579788.5') \
                         .copy('ENST00000579788.5_0')

    transcript.update_polyA_site(64002595 - 100)
    transcript.df['Start'].iloc[0] == 64002595 - 100
    transcript.df['Start'].iloc[0] == 64002595 - 100


def test_GtfModifier_to_gtf(tmp_path, modifier):
    gtf_cor_path = tmp_path / 'corr.gtf'

    t1 = modifier.fetch_transcript('ENST00000308278.12').copy(
        'ENST00000308278.12_0')
    t2 = modifier.fetch_transcript('ENST00000308278.12').copy(
        'ENST00000308278.12_1')

    t1.update_start_site(732411 - 100)
    t1.update_polyA_site(742972 + 100)

    t2.update_start_site(732411 + 100)
    t2.update_polyA_site(742972 - 100)

    modifier.add_transcript(t1)
    modifier.add_transcript(t2)
    modifier.to_gtf(gtf_cor_path)

    df_gtf = pr.read_gtf(gtf_cor_path).df
    assert df_gtf.shape[0] == 13


def test_correct_talon(tmp_path, lapa_links_chr17):

    output_gtf = tmp_path / 'chr17_corrected.gtf'

    correct_talon(lapa_links_chr17, read_annot_gm12_pb,
                  gtf_gm12_pb, output_gtf, chr17_abundance)

    df_gtf_input = pr.read_gtf(gtf_gm12_pb).df
    df_gtf = pr.read_gtf(output_gtf).df

    assert len(df_gtf_input['gene_id'].unique()) >= len(
        df_gtf['gene_id'].unique())
