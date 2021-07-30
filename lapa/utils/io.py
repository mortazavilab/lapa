import numpy as np
import pandas as pd
import pyranges as pr
from lapa.utils.common import filter_main_chroms, \
    chroms, chroms_chr


def remove_chr(df):
    df['Chromosome'] = df['Chromosome'].str.replace('chr', '')
    return df


def read_talon_read_annot(path):
    df = pd.read_csv(path, sep='\t')
    df = df.rename(columns={
        'chrom': 'Chromosome',
        'read_start': 'Start',
        'read_end': 'End',
        'strand': 'Strand',
        'annot_gene_id': 'gene_id',
        'dataset': 'sample'
    })
    # TO FIX: start end strand may not be correct
    return df


def read_bam_ends(path):
    df = pr.read_bam(path)

    import pdb
    pdb.set_trace()


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
    cols = [
        'Chromosome', 'Start', 'End', 'polyA_site', 'count', 'Strand',
        'fracA', 'singal', 'Feature', 'canonical_site', 'canonical', 'tpm'
    ]
    df = pd.read_csv(path, header=None, sep='\t')
    df.columns = cols
    return df


def read_apa_sample(path):
    cols = [
        'Chromosome', 'Start', 'End', 'count', 'polyA_site', 'Strand',
        'gene_id', 'gene_count', 'usage', 'count_cluster',
        'fracA', 'singal', 'Feature', 'canonical_site', 'canonical', 'tpm'
    ]
    df = pd.read_csv(path, header=None, sep='\t')
    df.columns = cols
    return df
