from pathlib import Path
import numpy as np
import pandas as pd
from bamread import read_bam
from lapa.utils.common import chroms_chr


cluster_col_order = [
    'Chromosome', 'Start', 'End', 'polyA_site', 'tpm', 'Strand',
    'Feature', 'count', 'fracA', 'signal', 'canonical_site'
]

sample_col_order = [
    'Chromosome', 'Start', 'End', 'polyA_site', 'tpm', 'Strand',
    'Feature', 'gene_id', 'count', 'gene_count', 'usage',
    'fracA', 'signal', 'canonical_site'
]


def remove_chr(df):
    df['Chromosome'] = df['Chromosome'].str.replace('chr', '')
    return df


def read_talon_read_annot(path, usecols=(
        'read_name', 'chrom', 'read_start', 'read_end',
        'strand', 'annot_gene_id', 'annot_transcript_id', 'dataset')):

    df = pd.read_csv(path, sep='\t', usecols=usecols)
    df = df.rename(columns={
        'chrom': 'Chromosome',
        'read_start': 'Start',
        'read_end': 'End',
        'strand': 'Strand',
        'annot_gene_id': 'gene_id',
        'dataset': 'sample'
    })

    start = np.where(df['Start'] < df['End'], df['Start'], df['End'])
    end = np.where(df['Start'] > df['End'], df['Start'], df['End'])
    df['Start'] = start.copy()
    df['End'] = end.copy()

    return df


def read_talon_read_annot_count(path):
    df = read_talon_read_annot(path)

    df['End'] = np.where(df['Strand'] == '-', df['Start'], df['End'])
    del df['Start']
    df['Start'] = df['End'] - 1
    df['count'] = 1

    return df[['Chromosome', 'Start', 'End', 'Strand', 'count', 'sample']]


def read_bam_ends(path, mapq=10, sample=None):
    df = read_bam(path, mapq=mapq)

    df['End'] = np.where(df['Strand'] == '-', df['Start'], df['End'])
    del df['Start']

    if sample is None:
        df['sample'] = Path(path).stem
    else:
        df['sample'] = sample

    return df


def read_sample_csv(path, mapq=10):
    df_sample = pd.read_csv(path)

    df = list()
    for i, row in df_sample.iterrows():
        df.append(read_bam_ends(row['path'], sample=row['sample']))

    return pd.concat(df)


def read_chrom_sizes(chrom_size_file):
    df_chrom_sizes = pd.read_csv(chrom_size_file, sep='\t', header=None)
    df_chrom_sizes = df_chrom_sizes[df_chrom_sizes[0].isin(chroms_chr)]
    chrom_sizes = df_chrom_sizes.set_index(0)[1].to_dict()
    # return {k.replace('chr', ''): v for k, v in chrom_sizes.items()}
    return chrom_sizes


def bw_from_pyranges(gr, value_col, chrom_size_file, bw_pos_file, bw_neg_file):
    chrom_sizes = read_chrom_sizes(chrom_size_file)

    gr = gr[gr.Chromosome.isin(chrom_sizes.keys())]

    if gr['+'].length > 0:
        gr['+'].to_bigwig(bw_pos_file, chromosome_sizes=chrom_sizes,
                          rpm=False, value_col=value_col)
    if gr['-'].length > 0:
        gr['-'].to_bigwig(bw_neg_file, chromosome_sizes=chrom_sizes,
                          rpm=False, value_col=value_col)


def read_polyA_cluster(path):
    df = pd.read_csv(path, header=None, sep='\t')
    df.columns = cluster_col_order
    return df


def read_apa_sample(path):
    df = pd.read_csv(path, header=None, sep='\t')
    df.columns = sample_col_order
    return df
