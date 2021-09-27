import pandas as pd
import pyranges as pr
from tqdm import tqdm
from conftest import read_annot, fasta, gtf, chrom_sizes, \
    read_annot_gm12_pb, gtf_gm12_pb
from lapa import lapa
from lapa.utils.io import read_talon_read_annot
from lapa.read import read_tes_mapping, correct_gtf_tes, sort_gtf


tqdm.pandas()


def test_read_tes_mapping(tmp_path):
    output_dir = tmp_path / 'lapa'
    lapa(read_annot, fasta, gtf, chrom_sizes, output_dir)

    df_read_annot = read_talon_read_annot(read_annot)
    df_mapping = read_tes_mapping(output_dir, read_annot)

    assert (df_mapping.shape[0] / df_read_annot.shape[0]) > 0.9


def test_correct_gtf_tes(tmp_path):
    df_read_annot = read_talon_read_annot(read_annot_gm12_pb)

    output_dir = tmp_path / 'lapa'
    lapa(read_annot_gm12_pb, fasta, gtf, chrom_sizes, output_dir)

    df_mapping = read_tes_mapping(output_dir, read_annot_gm12_pb)
    gtf_corrected = tmp_path / 'corrected.gtf'

    df_read_transcript = df_read_annot.rename(
        columns={'annot_transcript_id': 'transcript_id'})[
        ['read_name', 'transcript_id']]
    correct_gtf_tes(df_mapping, df_read_transcript, gtf_gm12_pb, gtf_corrected)

    df_gtf = pr.read_gtf(gtf_gm12_pb).df
    df_gtf_corrected = pr.read_gtf(gtf_corrected).df

    assert set(df_gtf['gene_id']) == set(df_gtf_corrected['gene_id'])
    assert set(df_gtf['exon_id']) == set(df_gtf_corrected['exon_id'])
    assert set(df_gtf['transcript_id']) == set(
        df_gtf_corrected['transcript_id'].str.split('_').str.get(0))

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


def test_correct_gtf_tes_NBR1(tmp_path):

    transripts_id = 'ENCODEHT000443413'

    df_read_annot = read_talon_read_annot(read_annot_gm12_pb)
    df_read_annot = df_read_annot[
        df_read_annot['annot_transcript_name'] == transripts_id]

    read_annot_nbr1 = str(tmp_path / 'nbr1_read_annot.tsv')
    df_read_annot.to_csv(read_annot_nbr1, sep='\t', index=False)

    output_dir = tmp_path / 'lapa'
    lapa(read_annot_nbr1, fasta, gtf, chrom_sizes, output_dir)

    df_mapping = read_tes_mapping(output_dir, read_annot_gm12_pb)

    gtf_corrected = tmp_path / 'corrected.gtf'

    df_read_transcript = df_read_annot.rename(
        columns={'annot_transcript_id': 'transcript_id'})[
        ['read_name', 'transcript_id']]
    correct_gtf_tes(df_mapping, df_read_transcript, gtf_gm12_pb, gtf_corrected)

    df_gtf = pr.read_gtf(gtf_gm12_pb).df
    df_gtf_corrected = pr.read_gtf(gtf_corrected).df

    df_gtf_sort = sort_gtf(
        df_gtf[df_gtf['gene_name'] != 'NBR1']).reset_index(drop=True)
    df_gtf_sort_corr = sort_gtf(df_gtf_corrected[
        df_gtf_corrected['gene_name'] != 'NBR1'][
            df_gtf.columns]).reset_index(drop=True)

    assert all(df_gtf_sort_corr['Start'] == df_gtf_sort['Start'])
    # TALON have transcripts longer than genes
    assert all(
        df_gtf_sort_corr[df_gtf_sort_corr['Feature'] == 'transcript']['End']
        == df_gtf_sort[df_gtf_sort['Feature'] == 'transcript']['End'])

    _df = df_gtf_corrected[df_gtf_corrected['transcript_id'] == transripts_id]

    # last exon end and transcript end should be same
    transcript_end = _df.iloc[0]['End']
    assert all(_df[_df['exon_number'] == '17']['End'] == transcript_end)

    polyA_site = df_mapping['polyA_site'].unique()
    assert len(polyA_site) == 1
    assert transcript_end == polyA_site[0]
