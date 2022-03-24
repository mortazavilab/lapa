import numpy as np
import pandas as pd
import pyranges as pr
from tqdm import tqdm
from lapa.utils.io import read_talon_read_annot
from lapa.result import LapaResult


def read_tes_mapping(df_cluster, read_annot, distance=1000):
    if type(read_annot) == str:
        df_reads = read_talon_read_annot(read_annot)
    elif type(read_annot) == pd.DataFrame:
        df_reads = read_annot
    else:
        raise ValueError('`read_annot` should be `str` or `pd.DataFrame`')

    df_reads['End'] = np.where(df_reads['Strand'] == '-',
                               df_reads['Start'],
                               df_reads['End'])
    df_reads['Start'] = df_reads['End'] - 1
    gr_reads = pr.PyRanges(df_reads)

    df = gr_reads.nearest(pr.PyRanges(df_cluster), how='downstream',
                          strandedness='same').df
    df = df[df.Distance < distance]

    df = df_reads.set_index('read_name').join(
        df.set_index('read_name')[['polyA_site']])

    unclustered = df['polyA_site'].isna()
    df.loc[unclustered, 'polyA_site'] = df.loc[unclustered, 'End']
    # df.loc[unclustered, 'count'] = 1

    df['polyA_site'] = df['polyA_site'].astype(int)

    return df.reset_index()[['read_name', 'Chromosome',
                             'polyA_site', 'Strand']]


def tss_mapping(df_tss_cluster, read_annot, distance=1000):
    df_reads = read_tss_read_annot(read_annot)
    gr_reads = pr.PyRanges(df_reads)

    df = gr_reads.nearest(pr.PyRanges(df_tss_cluster), how='upstream').df
    unclustered = df.Distance > distance
    df.loc[unclustered, 'start_site'] = df.loc[unclustered, 'End']
    df.loc[unclustered, 'count'] = 1

    df['start_site'] = df['start_site'].astype(int)

    return df[['read_name', 'Chromosome', 'start_site', 'Strand']]


def _correct_transcript(df):
    df = df.reset_index(drop=True)

    if all(df['polyA_site'].isna()):
        if 'start_site' in df.columns:
            if all(df['start_site'].isna()):
                return df
        else:
            return df

    strand = df.iloc[0]['Strand']
    _df_exon = df[df['Feature'] == 'exon']
    first_exon = _df_exon['Start'].idxmin()
    last_exon = _df_exon['End'].idxmax()

    transcript = df['Feature'] == 'transcript'
    min_exon_len = 25

    # Update transcripts
    polyA_site = df.loc[transcript, 'polyA_site']
    if not all(polyA_site.isna()):
        polyA_site = polyA_site.astype(int)[0]
        if strand == '+':
            if df.loc[last_exon, 'Start'] < polyA_site - min_exon_len:
                df.loc[last_exon, 'End'] = polyA_site
                df.loc[transcript, 'End'] = polyA_site
        elif strand == '-':
            if polyA_site < df.loc[first_exon, 'End'] - min_exon_len:
                df.loc[first_exon, 'Start'] = polyA_site
                df.loc[transcript, 'Start'] = polyA_site

    if 'start_site' in df.columns:
        start_site = df.loc[transcript, 'start_site']
        if not all(start_site.isna()):
            start_site = start_site.astype(int)[0]
            if strand == '+':
                if start_site < df.loc[first_exon, 'End'] - min_exon_len:
                    df.loc[first_exon, 'Start'] = start_site
                    df.loc[transcript, 'Start'] = start_site
            elif strand == '-':
                if df.loc[last_exon, 'Start'] < start_site - min_exon_len:
                    df.loc[last_exon, 'End'] = start_site
                    df.loc[transcript, 'End'] = start_site

    return df


def _correct_gene(df_transcript, df_gene):
    df = df_transcript.groupby('gene_id').agg({'Start': 'min', 'End': 'max'})
    df_join = df_gene.set_index('gene_id').join(df, rsuffix='_b')

    df_join['Start'] = np.minimum(df_join['Start'], df_join['Start_b'])
    df_join['End'] = np.maximum(df_join['End'], df_join['End_b'])

    del df_join['Start_b']
    del df_join['End_b']

    return df_join.reset_index()


def _sort_gtf_key(col):
    if col.name == 'End':
        return -col
    elif col.name == 'exon_number':
        return col.astype(float)
    else:
        return col


def sort_gtf(df):
    return df.sort_values([
        'Chromosome', 'gene_id', 'transcript_id', 'exon_number',
        'Start', 'End'
    ], na_position='first', key=_sort_gtf_key)


def _transcript_mapping(df_read_tes, df_read_transcript, site, minor_threshold=0.5):

    df = df_read_tes.set_index('read_name').join(
        df_read_transcript.set_index('read_name')
    ).reset_index()
    df['Strand'] = df['Strand'].astype('str')

    df = df.groupby(['transcript_id', 'Strand', site]).agg(
        count=('read_name', 'count')).reset_index(site)

    df = df.join(df.groupby(['transcript_id', 'Strand'])[['count']].agg('max')
                 .rename(columns={'count': 'threshold'}) * minor_threshold)

    return df[df['count'] > df['threshold']].reset_index()


def tes_transcript_mapping(df_read_tes, df_read_transcript,
                           minor_threshold=0.5):
    site = 'polyA_site'
    df = _transcript_mapping(df_read_tes, df_read_transcript, site,
                             minor_threshold=minor_threshold)
    return pd.concat([
        df[df['Strand'] == '-'].groupby('transcript_id')[[site]].agg('min'),
        df[df['Strand'] == '+'].groupby('transcript_id')[[site]].agg('max')
    ])


def tss_transcript_mapping(df_read_tss, df_read_transcript,
                           minor_threshold=0.5):
    site = 'start_site'
    df = _transcript_mapping(df_read_tss, df_read_transcript, site,
                             minor_threshold=minor_threshold)
    return pd.concat([
        df[df['Strand'] == '-'].groupby('transcript_id')[[site]].agg('max'),
        df[df['Strand'] == '+'].groupby('transcript_id')[[site]].agg('min')
    ])


def correct_gtf_tes(df_read_tes, df_read_transcript, gtf, gtf_output,
                    df_read_tss=None):

    df_tes = tes_transcript_mapping(df_read_tes, df_read_transcript)

    if df_read_tss is not None:
        df_tss = tss_transcript_mapping(df_read_tss, df_read_transcript)

    df_gtf = pr.read_gtf(gtf).df
    df_transcript = df_gtf[df_gtf['Feature'].isin({'transcript', 'exon'})]

    df_transcript = df_transcript.set_index('transcript_id') \
                                 .join(df_tes, how='left')

    if df_read_tss is not None:
        df_transcript = df_transcript.join(df_tss, how='left')

    tqdm.pandas()
    df_transcript_cor = df_transcript.groupby(level=0).progress_apply(
        _correct_transcript).reset_index()

    del df_transcript_cor['polyA_site']

    if df_read_tss is not None:
        del df_transcript_cor['start_site']

    df_gene_cor = _correct_gene(
        df_transcript_cor, df_gtf[df_gtf['Feature'] == 'gene'])

    df_gtf_cor = sort_gtf(pd.concat([df_gene_cor, df_transcript_cor]))

    pr.PyRanges(df_gtf_cor).to_gtf(gtf_output)


def correct_gtf(gtf, gtf_output, lapa_dir, lapa_tss_dir,
                read_annot, fasta, tss_correct=True):
    print('TES read mapping (1 / 3)...')
    df_cluster = LapaResult(lapa_dir).read_cluster()
    df_mapping = read_tes_mapping(df_cluster, read_annot)

    print('TSS read mapping (2 / 3)...')  
    df_tss_cluster = pd.read_csv(Path(lapa_tss_dir) / 'tss_clusters.bed')
    df_tss_mapping = tss_mapping(df_tss_cluster, read_annot)

    df_reads = read_talon_read_annot(read_annot).rename(
        columns={'annot_transcript_id': 'transcript_id'})

    print('Correcting gtf (3 / 3)...')
    correct_gtf_tes(
        df_mapping,
        df_reads[['read_name', 'transcript_id']],
        gtf,
        gtf_output,
        df_tss_mapping
    )
