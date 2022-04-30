from pathlib import Path
import pandas as pd
import pyranges as pr
from lapa.cluster import PolyAClustering, TssClustering
from lapa.genomic_regions import GenomicRegions
from lapa.count import TesMultiCounter, TssMultiCounter
from lapa.utils.io import cluster_col_order, sample_col_order


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


def tss_sample(df_cluster, df_tss_sample):
    columns = ['Chromosome', 'Start', 'End', 'Strand']
    gr = pr.PyRanges(df_tss_sample)

    gr_cluster = pr.PyRanges(df_cluster[columns])
    df_join = gr_cluster.join(gr, suffix='_sample').df

    df_join = df_join.groupby(columns, observed=True) \
                     .agg({'count': 'sum'}).reset_index()
    df_join = df_join[df_join['count'] > 0]

    df_tss = df_join.set_index(columns).join(
        df_cluster.set_index(columns), rsuffix='_cluster').reset_index()

    cols = ['Chromosome', 'Start', 'End', 'peak', 'count', 'Strand']
    return df_tss[cols]


def lapa(alignment, fasta, annotation, chrom_sizes, output_dir, method='end',
         min_tail_len=10, min_percent_a=0.9, mapq=10,
         cluster_extent_cutoff=3, cluster_window=25, cluster_ratio_cutoff=0.05):

    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    print('Counting TES (1 / 4)...')
    counter = TesMultiCounter(alignment, method, mapq,
                              min_tail_len, min_percent_a)
    df_tes, tes = counter.to_df()
    counter._to_bigwig(df_tes, tes, chrom_sizes,
                       output_dir, prefix='tes_counts')

    print('Clustering TES and calculating polyA_sites (2 / 4)...')
    df_cluster = PolyAClustering(fasta,
                                 extent_cutoff=cluster_extent_cutoff,
                                 ratio_cutoff=cluster_ratio_cutoff,
                                 window=cluster_window).to_df(df_tes)
    del df_tes

    print('Annotating TES cluster (3 / 4)...')
    df_cluster = tes_cluster_annotate(df_cluster, annotation)
    df_cluster[cluster_col_order].to_csv(output_dir / 'polyA_clusters.bed',
                                         index=False, sep='\t', header=False)

    print('Calculationg APA per samples (4 / 4)...')
    for sample, df_tes_sample in tes.items():
        df_apa = tes_sample(df_cluster, df_tes_sample)
        df_apa[sample_col_order].to_csv(output_dir / f'{sample}_apa.bed',
                                        index=False, sep='\t', header=False)


def lapa_tss(alignment, fasta, annotation, chrom_sizes, output_dir,
             method='start', mapq=10, cluster_extent_cutoff=3, cluster_window=25,
             cluster_ratio_cutoff=0.05):

    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    print('Counting TSS (1 / 4)...')
    counter = TssMultiCounter(alignment, method, mapq)
    df_tss, tss = counter.to_df()
    counter._to_bigwig(df_tss, tss, chrom_sizes,
                       output_dir, prefix='tss_counts')

    print('Clustering TES and calculating polyA_sites (2 / 4)...')
    df_cluster = TssClustering(fasta,
                               extent_cutoff=cluster_extent_cutoff,
                               window=cluster_window).to_df(df_tss)
    del df_tss

    cluster_cols = ['Chromosome', 'Start', 'End', 'peak', 'count', 'Strand']
    df_cluster[cluster_cols].to_csv(output_dir / 'tss_clusters.bed',
                                    index=False, sep='\t', header=False)

    for sample, df_tss_sample in tss.items():
        df_tss = tss_sample(df_cluster, df_tss_sample)
        df_tss[cluster_cols].to_csv(output_dir / f'{sample}_tss.bed',
                                    index=False, sep='\t', header=False)
