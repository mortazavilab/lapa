from pathlib import Path
import numpy as np
import pyranges as pr
import pandas as pd
from lapa.result import LapaResult
from lapa.utils.io import read_talon_read_annot, \
    read_tss_cluster


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


def _link_reads_to_tes(df_tes_cluster, df_reads, distance=1000):
    '''
    Link reads to transcript end sites
    '''
    df_reads['End'] = np.where(df_reads['Strand'] == '-',
                               df_reads['read_Start'],
                               df_reads['read_End'])
    df_reads['Start'] = df_reads['End'] - 1
    gr_reads = pr.PyRanges(df_reads)

    df = gr_reads.nearest(pr.PyRanges(df_tes_cluster), how='downstream',
                          strandedness='same').df
    df.loc[df['Distance'] > distance, 'polyA_site'] = -1
    return df[[*_reads_cols, 'polyA_site']]


def _link_reads_to_tss(df_tss_cluster, df_reads, distance=1000):
    '''
    Link reads to transcript start sites
    '''
    df_reads['End'] = np.where(df_reads['Strand'] == '-',
                               df_reads['read_End'],
                               df_reads['read_Start'])
    df_reads['Start'] = df_reads['End'] - 1
    gr_reads = pr.PyRanges(df_reads)

    df = gr_reads.nearest(pr.PyRanges(df_tss_cluster), how='upstream',
                          strandedness='same').df
    df.loc[df['Distance'] > distance, 'start_site'] = -1
    return df[[*_reads_cols, 'start_site', 'polyA_site']]


def link_tss_to_tes(alignment, lapa_dir, lapa_tss_dir, distance=1000,
                    mapq=10, min_read_length=100):
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

    print('TES read mapping...')
    df_cluster = LapaResult(lapa_dir, tpm_cutoff=1).read_cluster()
    df_reads = _link_reads_to_tes(df_cluster, df_reads, distance=distance)

    print('TSS read mapping...')
    df_tss_cluster = read_tss_cluster(
        Path(lapa_tss_dir) / 'tss_clusters.bed')
    df_reads = _link_reads_to_tss(df_tss_cluster, df_reads, distance=distance)

    return df_reads
