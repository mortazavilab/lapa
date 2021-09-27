import numpy as np
import pandas as pd
import pyranges as pr
from tqdm import tqdm
from lapa.utils.io import read_talon_read_annot
from lapa.cluster import TesClustering
from lapa.result import LapaResult


tqdm.pandas()


def read_tes_mapping(df_cluster, read_annot, filter_internal_priming=True,
                     distance=500):
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

    gr = gr_reads.nearest(pr.PyRanges(df_cluster), how='downstream')
    df = gr[gr.Distance < distance].df
    df['polyA_site'] = df['polyA_site'].astype(int)

    return df[['read_name', 'Chromosome', 'polyA_site', 'Strand', 'Distance']]


def read_tss_read_annot(read_annot):
    df = read_talon_read_annot(read_annot)
    df['End'] = np.where(df['Strand'] == '-', df['End'], df['Start'])
    del df['Start']
    df['Start'] = df['End'] - 1
    return df


def read_tss_read_annot_count(read_annot):
    df = read_tss_read_annot(read_annot)
    df['count'] = 1

    columns = ['Chromosome', 'Start', 'End', 'Strand']
    df = df.groupby(columns).agg('sum').reset_index()
    return df


def tss_cluster(read_annot, fasta):
    df_reads = read_tss_read_annot_count(read_annot)

    tss_clusters = list(TesClustering(fasta).cluster(df_reads))
    return pd.DataFrame({
        'Chromosome': [c.Chromosome for c in tss_clusters],
        'Start': [c.Start for c in tss_clusters],
        'End': [c.End for c in tss_clusters],
        'Strand': [c.Strand for c in tss_clusters],
        'count': [c.total_count for c in tss_clusters],
        'start_site': [c.polyA_site() for c in tss_clusters]
    })


def tss_mapping(df_tss_cluster, read_annot, distance=500):
    df_reads = read_tss_read_annot(read_annot)
    gr_reads = pr.PyRanges(df_reads)

    gr = gr_reads.nearest(pr.PyRanges(df_tss_cluster), how='upstream')
    df = gr[gr.Distance < distance].df
    df['start_site'] = df['start_site'].astype(int)

    return df[['read_name', 'Chromosome', 'start_site', 'Strand', 'Distance']]


def _correct_transcript(df):
    df = df.reset_index()

    if all(df['polyA_site'].isna()) and all(df['start_site'].isna()):
        return df

    if df.iloc[0]['Strand'] == '+':
        tss_pos = 'Start'
        tes_pos = 'End'
    else:
        tss_pos = 'End'
        tes_pos = 'Start'

    # Update transcripts
    transcript = df['Feature'] == 'transcript'

    polyA_site = df.loc[transcript, 'polyA_site']
    if not all(polyA_site.isna()):
        df.loc[transcript, tes_pos] = polyA_site.astype(int)
    start_site = df.loc[transcript, 'start_site']
    if not all(start_site.isna()):
        df.loc[transcript, tss_pos] = start_site.astype(int)

    # Update exons boundries
    exons = df['Feature'] == 'exon'
    _df_exon = df[exons]
    _df_transcript = df[df['Feature'] == 'transcript']

    first_exon = _df_exon['Start'].idxmin()
    last_exon = _df_exon['End'].idxmax()

    df.loc[first_exon, 'Start'] = _df_transcript.iloc[0]['Start']
    df.loc[last_exon, 'End'] = _df_transcript.iloc[0]['End']

    return df


def _correct_gene(df_transcript, df_gene):
    df = df_transcript.groupby('gene_id').agg({'Start': 'min', 'End': 'max'})
    df_join = df_gene.set_index('gene_id').join(df, rsuffix='_b')

    df_join['Start'] = np.minimum(df_join['Start'], df_join['Start_b'])
    df_join['End'] = np.maximum(df_join['End'], df_join['End_b'])

    del df_join['Start_b']
    del df_join['End_b']

    return df_join.reset_index()


def sort_gtf(df):
    return df.sort_values([
        'Chromosome', 'gene_id', 'transcript_id',  'Start', 'End'
    ], na_position='first', key=lambda x: x if x.name != 'End' else -x)


def _transcript_mapping(df_read_tes, df_read_transcript, site):
    df = df_read_tes.set_index('read_name').join(
        df_read_transcript.set_index('read_name')
    ).reset_index()
    df['Strand'] = df['Strand'].astype('str')

    df = df.groupby(['transcript_id', 'Strand', site]).agg(
        count=('read_name', 'count')).reset_index(site)

    df = df.join(df.groupby(['transcript_id', 'Strand'])[['count']].agg('max')
                 .rename(columns={'count': 'threshold'}) * 0.5)

    df = df[df['count'] > df['threshold']].reset_index()
    return df


def tes_transcript_mapping(df_read_tes, df_read_transcript):
    site = 'polyA_site'
    df = _transcript_mapping(df_read_tes, df_read_transcript, site)
    return pd.concat([
        df[df['Strand'] == '-'].groupby('transcript_id')[[site]].agg('min'),
        df[df['Strand'] == '+'].groupby('transcript_id')[[site]].agg('max')
    ])


def tss_transcript_mapping(df_read_tss, df_read_transcript):
    site = 'start_site'
    df = _transcript_mapping(df_read_tss, df_read_transcript, site)
    return pd.concat([
        df[df['Strand'] == '-'].groupby('transcript_id')[[site]].agg('max'),
        df[df['Strand'] == '+'].groupby('transcript_id')[[site]].agg('min')
    ])


def correct_gtf_tes(df_read_tes, df_read_transcript, gtf, gtf_output,
                    df_read_tss=None):

    df_tes = tes_transcript_mapping(df_read_tes, df_read_transcript)
    df_tss = tss_transcript_mapping(df_read_tss, df_read_transcript)

    df_gtf = pr.read_gtf(gtf).df
    df_transcript = df_gtf[df_gtf['Feature'].isin({'transcript', 'exon'})]
    df_transcript = df_transcript.set_index('transcript_id') \
                                 .join(df_tes, how='left') \
                                 .join(df_tss, how='left')

    df_transcript_cor = df_transcript.groupby('transcript_id').progress_apply(
        _correct_transcript).reset_index(drop=True)

    del df_transcript_cor['polyA_site']

    df_gene_cor = _correct_gene(df_transcript_cor,
                                df_gtf[df_gtf['Feature'] == 'gene'])

    df_gtf_cor = sort_gtf(pd.concat([df_gene_cor, df_transcript_cor]))
    pr.PyRanges(df_gtf_cor).to_gtf(gtf_output)
