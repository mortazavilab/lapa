import numpy as np
import pandas as pd
import pyranges as pr
from tqdm import tqdm
from lapa.utils.io import read_talon_read_annot
from lapa.result import LapaResult


tqdm.pandas()


def read_tes_mapping(lapa_dir, read_annot, filter_internal_priming=True,
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

    df_cluster = LapaResult(lapa_dir).read_cluster()

    if filter_internal_priming:
        df_cluster = df_cluster[(~(
            (df_cluster['fracA'] > 7) &
            (df_cluster['signal'] == 'None@None')))
            | (df_cluster['canonical_site'] != -1)
        ]

    gr = gr_reads.nearest(pr.PyRanges(df_cluster), how='downstream')
    df = gr[gr.Distance < distance].df
    df['polyA_site'] = df['polyA_site'].astype(int)

    return df[['read_name', 'Chromosome', 'polyA_site', 'Strand', 'Distance']]


def _correct_transcript(df):
    df = df.reset_index()

    if all(df['polyA_site'].isna()):
        return df

    pos = 'End' if df.iloc[0]['Strand'] == '+' else 'Start'
    transcript = df['Feature'] == 'transcript'
    df.loc[transcript, pos] = df.loc[transcript, 'polyA_site'].astype(int)

    # TODO: update start based on tss_site

    # dfs = list()

    # for i, (_, _df) in enumerate(df.groupby('polyA_site')):
    #     _df['transcript_id'] = _df['transcript_id'] + f'_{i}'

    # exons = _df['Feature'] == 'exon'
    # _df_exon = _df[exons]
    # _df_transcript = _df[_df['Feature'] == 'transcript']

    # first_exon = _df_exon['Start'].idxmin()
    # last_exon = _df_exon['End'].idxmax()

    # _df.loc[first_exon, 'Start'] = _df_transcript.iloc[0]['Start']
    # _df.loc[last_exon, 'End'] = _df_transcript.iloc[0]['End']

    # dfs.append(_df)

    # return pd.concat(dfs)

    exons = df['Feature'] == 'exon'
    _df_exon = df[exons]
    _df_transcript = df[df['Feature'] == 'transcript']

    first_exon = _df_exon['Start'].idxmin()
    last_exon = _df_exon['End'].idxmax()

    df.loc[first_exon, 'Start'] = _df_transcript.iloc[0]['Start']
    df.loc[last_exon, 'End'] = _df_transcript.iloc[0]['End']

    return df

# def _correct_exon(df_transcript, df_exon):

#     df_transcript['transcript_id_uncor'] = df_transcript[
#         'transcript_id'].str.split('_').str.get(0)

#     import pdb
#     pdb.set_trace()

#     df_exon = df_exon.set_index('transcript_id') \
#         .join(df_exon.groupby('transcript_id')
#               .agg({'exon_number': 'max'})
#               .rename(columns={'exon_number': 'last_exon'})) \
#         .join(df_transcript.set_index('transcript_id')[['Start', 'End']],
#               rsuffix='_t') \
#         .reset_index()

#     df_exon['last_exon'] = df['last_exon'] == df['exon_number']
#     df_exon['first_exon'] = df['exon_number'] == 1

#     df_exon.loc[df['first_exon'], 'Start'] = df_exon.loc[
#         df_exon['first_exon'], 'Start_t']
#     df_exon.loc[df['last_exon'], 'End'] = df_exon.loc[
#         df_exon['first_exon'], 'End_t']

#     del df['Start_t']
#     del df['End_t']

#     return df_exon


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


def correct_gtf_tes(df_read_tes, df_read_transcript, gtf, gtf_output):

    # cols = ['Chromosome', 'polyA_site', 'Strand', 'transcript_id']
    df = df_read_tes.set_index('read_name').join(
        df_read_transcript.set_index('read_name')
    )

    indexes = df.groupby(['transcript_id', 'polyA_site'])['Strand'] \
        .count().groupby(level=0).idxmax()
    df = pd.DataFrame({
        'transcript_id': indexes.str.get(0).tolist(),
        'polyA_site': indexes.str.get(1).tolist()
    }).set_index('transcript_id')

    # df = df[cols].drop_duplicates().set_index('transcript_id')

    df_gtf = pr.read_gtf(gtf).df
    df_transcript = df_gtf[df_gtf['Feature'].isin({'transcript', 'exon'})]
    df_transcript = df_transcript.set_index('transcript_id').join(
        df, rsuffix='_b', how='left')

    df_transcript_cor = df_transcript.groupby(
        'transcript_id').progress_apply(_correct_transcript).reset_index(drop=True)

    del df_transcript_cor['polyA_site']

    df_gene_cor = _correct_gene(df_transcript_cor,
                                df_gtf[df_gtf['Feature'] == 'gene'])

    df_gtf_cor = sort_gtf(pd.concat([df_gene_cor, df_transcript_cor]))
    pr.PyRanges(df_gtf_cor).to_gtf(gtf_output)
