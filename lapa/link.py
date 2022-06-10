from pathlib import Path
import numpy as np
import pyranges as pr
import pandas as pd
from lapa.utils.io import read_talon_read_annot
from lapa.result import LapaResult, LapaTssResult


_reads_cols = ['read_name', 'Chromosome', 'read_Start', 'read_End', 'Strand']


def _read_alignment_reads(alignment, mapq=10, min_read_length=100):
    '''
    Read `read_name` and locus of reads from alignment file.
    '''
    if str(alignment).endswith('_read_annot.tsv'):
        df = read_talon_read_annot(alignment).rename(
            columns={'Start': 'read_Start', 'End': 'read_End'}
        )[_reads_cols]

    elif str(alignment).endswith('.bam'):
        df = pd.concat([
            pr.read_bam(i, sparse=False, as_df=True, mapq=mapq)
            .rename(
                columns={'Start': 'read_Start',
                         'End': 'read_End',
                         'Name': 'read_name'}
            )[_reads_cols]
            for i in str(alignment).split(',')
        ])
    else:
        raise ValueError('Only bam and TALON read_annot'
                         ' file formats are supported.')

    return df[(df['read_End'] - df['read_Start']) > min_read_length]


def _link_reads_to_tes(df_tes_cluster, df_reads, distance=50):
    '''
    Link reads to transcript end sites
    '''
    df_reads['End'] = np.where(df_reads['Strand'] == '-',
                               df_reads['read_Start'],
                               df_reads['read_End'])
    df_reads['Start'] = df_reads['End'] - 1
    gr_reads = pr.PyRanges(df_reads)

    df = gr_reads.nearest(pr.PyRanges(df_tes_cluster),
                          strandedness='same').df
    df.loc[df['Distance'] > distance, 'polyA_site'] = -1
    return df[[*_reads_cols, 'polyA_site']]


def _link_reads_to_tss(df_tss_cluster, df_reads, distance=50):
    '''
    Link reads to transcript start sites
    '''
    df_reads['End'] = np.where(df_reads['Strand'] == '-',
                               df_reads['read_End'],
                               df_reads['read_Start'])
    df_reads['Start'] = df_reads['End'] - 1
    gr_reads = pr.PyRanges(df_reads)

    df = gr_reads.nearest(pr.PyRanges(df_tss_cluster),
                          strandedness='same').df
    df.loc[df['Distance'] > distance, 'tss_site'] = -1
    return df[[*_reads_cols, 'tss_site', 'polyA_site']]


def link_tss_to_tes(alignment, lapa_dir, lapa_tss_dir, distance=50,
                    mapq=10, min_read_length=100, dataset='all'):
    '''
    Link transcript site sites to transcript end sites using
    long-read from the alignment file.

    Args:
      alignment (str): Path to bam file or TALON read_annot file.
      lapa_dir (str): Path to lapa output directory generated
        with `lapa` command
      lapa_tss_dir (str): Path to lapa tss directory
        with `lapa_tss` command
    '''
    print('Reading alignment file...')
    df_reads = _read_alignment_reads(alignment, mapq=mapq,
                                     min_read_length=min_read_length)

    lapa = LapaResult(lapa_dir)
    lapa_tss = LapaTssResult(lapa_tss_dir)
    
    print('TES read mapping...')
    if dataset == 'raw_all':
        lapa.replicated = False
        df_cluster = lapa.read_clusters()
    elif dataset == 'all':
        df_cluster = lapa.read_clusters()
    elif dataset in lapa.datasets:
        df_cluster = lapa.read_dataset(dataset)
    else:
        raise ValueError(f'{dataset} is not in datasets.'
                         'Valid options are `all`, `raw_all`, or dataset name')

    core_cols = ['Chromosome', 'Start', 'End', 'Strand']
    df_cluster = df_cluster.drop_duplicates(core_cols)
    df_reads = _link_reads_to_tes(df_cluster, df_reads, distance=distance)

    print('TSS read mapping...')
    if dataset == 'raw_all':
        lapa_tss.replicated = False
        df_tss_cluster = lapa_tss.read_clusters()
    elif dataset == 'all':
        df_tss_cluster = lapa_tss.read_clusters()
    elif dataset in lapa_tss.datasets:
        df_tss_cluster = lapa_tss.read_dataset(dataset)
    else:
        raise ValueError(f'{dataset} is not in datasets.'
                         'Valid options are `all`, `raw_all`, or dataset name')

    df_tss_cluster = df_tss_cluster.drop_duplicates(core_cols)
    df_reads = _link_reads_to_tss(df_tss_cluster, df_reads, distance=distance)

    valid = np.where(df_reads['Strand'] == '+',
                     df_reads['tss_site'] < df_reads['polyA_site'],
                     df_reads['polyA_site'] < df_reads['tss_site'])
    valid = valid | (df_reads['tss_site'] == -1) \
        | (df_reads['polyA_site'] == -1)

    return df_reads[valid]
