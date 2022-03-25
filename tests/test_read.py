from pathlib import Path
import pandas as pd
import pyranges as pr
from tqdm import tqdm
from conftest import read_annot, fasta, gtf, chrom_sizes, \
    read_annot_gm12_pb, gtf_gm12_pb
from lapa import lapa, lapa_tss
from lapa.utils.io import read_talon_read_annot, read_tss_cluster
from lapa.read import read_tes_mapping, correct_gtf, sort_gtf, \
    read_tss_mapping, _transcript_mapping, read_tss_mapping, \
    tes_transcript_mapping
from lapa.result import LapaResult

tqdm.pandas()


def test_read_tes_mapping(tmp_path):
    output_dir = tmp_path / 'lapa'
    lapa(read_annot, fasta, gtf, chrom_sizes, output_dir)

    df_read_annot = read_talon_read_annot(read_annot)
    df_cluster = LapaResult(output_dir).read_cluster()
    df_mapping = read_tes_mapping(df_cluster, df_read_annot)

    assert (df_mapping.shape[0] / df_read_annot.shape[0]) > 0.9


def test_read_tes_mapping_uncluster():
    df_cluster = pd.DataFrame({
        'Chromosome': ['chr1', 'chr1'],
        'Start': [10000, 10000000],
        'End': [10010, 10000001],
        'Strand': ['+', '+'],
        'polyA_site': [10005, 10000000]
    })

    df_read = pd.DataFrame({
        'read_name': ['r1', 'r2', 'r3'],
        'Chromosome': ['chr1', 'chr1', 'chr1'],
        'Start': [5000, 6000, 90000],
        'End': [10001, 9750, 100000],
        'Strand': ['+', '+', '+']
    })

    df = read_tes_mapping(df_cluster, df_read)

    pd.testing.assert_frame_equal(df, pd.DataFrame({
        'read_name': ['r1', 'r2', 'r3'],
        'Chromosome': ['chr1', 'chr1', 'chr1'],
        'polyA_site': [10005, 10005, 100000],
        'Strand': ['+', '+', '+']
    }))


def test_read_tss_mapping(tmp_path):
    df_cluster = pd.DataFrame({
        'Chromosome': ['chr1', 'chr1'],
        'Start': [4990, 100],
        'End': [5010, 110],
        'Strand': ['+', '+'],
        'start_site': [5000, 105]
    })

    df_read = pd.DataFrame({
        'read_name': ['r1', 'r2', 'r3'],
        'Chromosome': ['chr1', 'chr1', 'chr1'],
        'Start': [5000, 6000, 90000],
        'End': [10001, 9750, 100000],
        'Strand': ['+', '+', '+'],
        'gene_id': ['g1', 'g1', 'g2'],
        'transcript_id': ['t1', 't1', 't2'],
        'sample': ['a', 'a', 'a']
    })

    df = read_tss_mapping(df_cluster, df_read)
    pd.testing.assert_frame_equal(
        df,
        pd.DataFrame({
            'read_name': ['r1', 'r2', 'r3'],
            'Chromosome': pd.Categorical(['chr1', 'chr1', 'chr1']),
            'start_site': [5000, 5000, 90000],
            'Strand': pd.Categorical(['+', '+', '+'])
        }))


def test__transcript_mapping():
    df_read_tes = pd.DataFrame({
        'read_name': ['r1', 'r2', 'r3', 'r4', 'r5', 'r6'],
        'Chromosome': 'chr1',
        'polyA_site': [10005, 10005, 11010, 11010, 10005, 100000],
        'Strand': '+'
    })

    df_read_transcript = pd.DataFrame({
        'read_name': ['r1', 'r2', 'r3', 'r4', 'r5', 'r6'],
        'transcript_id': ['t1', 't1', 't1', 't1', 't1', 't2']
    })

    df = _transcript_mapping(df_read_tes, df_read_transcript, 'polyA_site')

    pd.testing.assert_frame_equal(
        pd.DataFrame({
            'transcript_id': ['t1', 't1', 't2'],
            'Strand': ['+', '+', '+'],
            'polyA_site': [10005, 11010, 100000],
            'count': [3, 2, 1],
            'threshold': [1.5, 1.5, 0.5]
        }), df
    )


def test_tes_transcript_mapping():
    df_read_tes = pd.DataFrame({
        'read_name': ['r1', 'r2', 'r3', 'r4', 'r5', 'r6'],
        'Chromosome': 'chr1',
        'polyA_site': [10005, 10005, 11010, 11010, 10005, 100000],
        'Strand': '+'
    })

    df_read_transcript = pd.DataFrame({
        'read_name': ['r1', 'r2', 'r3', 'r4', 'r5', 'r6'],
        'transcript_id': ['t1', 't1', 't1', 't1', 't1', 't2']
    })

    df = tes_transcript_mapping(df_read_tes, df_read_transcript)
    pd.testing.assert_frame_equal(
        pd.DataFrame({
            'transcript_id': ['t1', 't2'],
            'polyA_site': [11010, 100000]
        }).set_index('transcript_id'), df
    )


def test_tss_transcript_mapping():
    pass


def test__correct_trancript():
    pass


def test__correct_gene():
    pass


def test_sort_gtf():
    pass


def test_correct_gtf_tes():
    pass


def test_correct_gtf_tes_integration(tmp_path):
    output_dir = tmp_path / 'lapa'
    lapa(read_annot_gm12_pb, fasta, gtf, chrom_sizes, output_dir)

    output_tss_dir = tmp_path / 'lapa_tss'

    lapa_tss(read_annot_gm12_pb, fasta, gtf, chrom_sizes, output_tss_dir)

    gtf_corrected = tmp_path / 'corrected.gtf'

    correct_gtf(gtf_gm12_pb, gtf_corrected, output_dir,
                output_tss_dir, read_annot_gm12_pb, fasta)

    df_gtf = pr.read_gtf(gtf_gm12_pb).df
    df_gtf_corrected = pr.read_gtf(gtf_corrected).df

    assert set(df_gtf['gene_id']) == set(df_gtf_corrected['gene_id'])
    assert set(df_gtf['exon_id']) == set(df_gtf_corrected['exon_id'])
    assert set(df_gtf['transcript_id']) == set(
        df_gtf_corrected['transcript_id'].str.split('_').str.get(0))

    # Gene bounders should be larger than transcripts
    _df_gene = df_gtf_corrected[df_gtf_corrected['Feature'] == 'gene']
    _df_transcript = df_gtf_corrected[
        df_gtf_corrected['Feature'] == 'transcript']
    _df_exon = df_gtf_corrected[df_gtf_corrected['Feature'] == 'exon']

    _df = _df_exon.set_index('gene_id').join(
        _df_gene.set_index('gene_id'), rsuffix='_gene')

    assert all(_df['Start'] >= _df['Start_gene'])
    assert all(_df['End'] <= _df['End_gene'])

    _df = _df_exon.set_index('transcript_id').join(
        _df_transcript.set_index('transcript_id'), rsuffix='_transcript')

    assert all(_df['Start'] >= _df['Start_transcript'])
    assert all(_df['End'] <= _df['End_transcript'])

    # Introns do not change after updates
    df_introns = pr.PyRanges(df_gtf).features.introns(
        by="transcript").df.sort_values('Start')
    df_introns_corr = pr.PyRanges(df_gtf_corrected).features.introns(
        by="transcript").df.sort_values('Start')

    df_introns = df_introns[['Chromosome', 'Start', 'End', 'Strand']] \
        .drop_duplicates()

    df_introns = df_introns.sort_values(
        df_introns.columns.tolist()).reset_index(drop=True)

    df_introns_corr = df_introns_corr[['Chromosome', 'Start',
                                       'End', 'Strand']] \
        .drop_duplicates()
    df_introns_corr = df_introns_corr.sort_values(
        df_introns_corr.columns.tolist()).reset_index(drop=True)

    pd.testing.assert_frame_equal(df_introns, df_introns_corr)
