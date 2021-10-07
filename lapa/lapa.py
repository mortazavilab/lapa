from pathlib import Path
import numpy as np
import pandas as pd
import pyranges as pr
from tqdm import tqdm
from kipoiseq import Interval
from kipoiseq.extractors import FastaStringExtractor
from lapa.utils.io import read_talon_read_annot_count, bw_from_pyranges, \
    read_bam_ends, cluster_col_order, sample_col_order, read_sample_csv
from lapa.utils.common import pad_series, polyA_signal_seqs
from lapa.genomic_regions import GenomicRegions
from lapa.count import count_tes_bam_samples, agg_tes_samples
from lapa.cluster import TesClustering


def tes_cluster_annotate(df_cluster, annotation):
    gr_apa = pr.PyRanges(df_cluster)

    total = gr_apa.count.sum()

    gr_gtf = pr.read_gtf(annotation)
    greg = GenomicRegions(gr_gtf)

    df = greg.annotate(gr_apa)

    df['tpm'] = df['count'] * 1000000 / total

    return df.sort_values(['Chromosome', 'End'])


def tes_sample(df_cluster, df_tes_sample,
               filter_intergenic=True,
               filter_internal_priming=True):

    columns = ['Chromosome', 'Start', 'End', 'Strand', 'gene_id']
    gr = pr.PyRanges(df_tes_sample)

    if filter_intergenic:
        df_cluster = df_cluster[df_cluster['Feature'] != 'intergenic']

    if filter_internal_priming:
        df_cluster = df_cluster[~((df_cluster['fracA'] > 7)
                                  & (df_cluster['signal'] == 'None@None'))]

    gr_cluster = pr.PyRanges(df_cluster[columns])
    df_join = gr_cluster.join(gr, suffix='_sample').df

    df_join = df_join.groupby(columns, observed=True) \
                     .agg({'count': 'sum'}).reset_index()
    df_join = df_join[df_join['count'] > 0]

    # get number of tes in the gene
    df_join = df_join.set_index('gene_id').join(
        df_join[['gene_id', 'count']]
        .groupby('gene_id').agg('sum').rename(
            columns={'count': 'gene_count'}))

    # calculate apa usage
    df_join['usage'] = df_join['count'] / df_join['gene_count']

    df_apa = df_join.reset_index().set_index(columns).join(
        df_cluster.set_index(columns), rsuffix='_cluster')

    df_apa['tpm'] = (df_apa['count'] * 1000000 /
                     df_apa['count'].sum()).round(2)
    return df_apa.reset_index()


def prepare_alignment(alignment):
    if alignment.endswith('.bam'):
        alignments = alignment.split(',')
        return pd.DataFrame({
            'sample': ['all'] * len(alignments),
            'path': alignments
        })
    elif alignment.endswith('.csv'):
        df = pd.read_csv(alignment)
        assert all(df.columns == ['sample', 'path']), \
            'provided csv file should be consist of columns `sample` and `path`'
        return df


def lapa(alignment, fasta, annotation, chrom_sizes, output_dir, method=None,
         min_tail_len=10, min_percent_a=0.9, mapq=10,
         cluster_extent_cutoff=3, cluster_window=25):

    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    print('Reading the alignment file...')
    if alignment.endswith('_read_annot.tsv'):
        df_count = read_talon_read_annot_count(alignment)
    elif alignment.endswith('.bam') or alignment.endswith('.csv'):
        if method not in {'tail', 'end'}:
            raise ValueError('method need to be either `tail` or `end`')
        df_alignment = prepare_alignment(alignment)
        df_count = count_tes_bam_samples(df_alignment, method, min_tail_len,
                                         min_percent_a, mapq)
    else:
        raise ValueError(
            'Unknown file alignment format: supported '
            'file formats are `bam` and `sample.csv`')

    print('Counting TES...')
    df_tes, tes = agg_tes_samples(df_count, chrom_sizes, output_dir)
    del df_count

    print('Clustering TES and calculating polyA_sites...')
    df_cluster = TesClustering(fasta,
                               extent_cutoff=cluster_extent_cutoff,
                               window=cluster_window).to_df(df_tes)
    del df_tes

    print('Annotating TES cluster...')
    df_cluster = tes_cluster_annotate(df_cluster, annotation)
    df_cluster[cluster_col_order].to_csv(output_dir / 'polyA_clusters.bed',
                                         index=False, sep='\t', header=False)

    print('Calculationg APA per samples...')
    for sample, df_tes_sample in tes.items():
        df_apa = tes_sample(df_cluster, df_tes_sample)
        df_apa[sample_col_order].to_csv(output_dir / f'{sample}_apa.bed',
                                        index=False, sep='\t', header=False)
