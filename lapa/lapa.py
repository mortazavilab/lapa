from pathlib import Path
import numpy as np
import pandas as pd
import pyranges as pr
from tqdm import tqdm
from kipoiseq import Interval
from kipoiseq.extractors import FastaStringExtractor
from lapa.utils.io import read_talon_read_annot, bw_from_pyranges, \
    read_bam_ends, cluster_col_order, sample_col_order, read_sample_csv
from lapa.utils.common import pad_series, polyA_signal_seqs
from lapa.genomic_regions import GenomicRegions


def _count_tes(df_alignment, chrom_sizes, sample, output_dir):
    # count TES sites
    columns = ['Chromosome', 'End', 'Strand']
    df = df_alignment[columns]

    df['count'] = 1
    df = df.groupby(columns).agg('sum').reset_index()

    df['Start'] = df['End'] - 1

    bw_from_pyranges(
        pr.PyRanges(df), 'count',
        chrom_sizes,
        str(output_dir / f'{sample}_tes_counts_pos.bw'),
        str(output_dir / f'{sample}_tes_counts_neg.bw')
    )

    return df


def count_tes(df_alignment, chrom_sizes, output_dir):
    tes = {
        sample: _count_tes(df, chrom_sizes, sample, output_dir)
        for sample, df in df_alignment.groupby('sample')
    }
    # TODO: combine tes samples counts pd.concat |> agg(sum)
    return _count_tes(df_alignment, chrom_sizes, 'all', output_dir), tes


class Cluster:

    def __init__(self, Chromosome: str, Start: int, End: int,
                 Strand: str, counts=None):
        self.Chromosome = str(Chromosome)
        self.Start = Start
        self.End = End
        self.Strand = Strand
        self.counts = counts or list()

    def extend(self, end, count):
        self.End = end
        self.counts.append(count)

    @property
    def total_count(self):
        return sum(self.counts)

    def __len__(self):
        return self.End - self.Start

    def polyA_site(self, window=5, std=1):
        pad_size = window // 2
        counts = pad_series(self.counts, pad_size=pad_size)
        moving_sum = counts.rolling(window, center=True,
                                    win_type='gaussian').sum(std=std)
        return self.Start + moving_sum.idxmax() + 1

    def polyA_signal_sequence(self, fasta, polyA_site):
        # the list of poly(A) signals
        # that are annotated for single 3' end processing sites
        # in the region -60 to +10 nt around them
        # According to Gruber et al., 2016, Genome Research.
        if self.Strand == '+':
            start = polyA_site - 60
            end = polyA_site + 10
        elif self.Strand == '-':
            start = polyA_site - 10
            end = polyA_site + 60

        interval = Interval(self.Chromosome, start, end,
                            strand=self.Strand)
        if self.Chromosome not in fasta.fasta or (interval.start < 0):
            return None, None

        seq = fasta.extract(interval)

        for signal_seq in polyA_signal_seqs:
            index = seq.find(signal_seq)
            if index == -1:
                continue
            else:
                return self.Start + index, signal_seq
        return None, None

    def fraction_A(self, fasta, polyA_site):
        if self.Strand == '+':
            start = polyA_site
            end = polyA_site + 10
        elif self.Strand == '-':
            start = polyA_site - 10
            end = polyA_site

        interval = Interval(self.Chromosome,
                            start, end, strand=self.Strand)
        if (self.Chromosome not in fasta.fasta) or (interval.start < 0):
            return -1

        seq = fasta.extract(interval)

        return sum('A' == i for i in seq)

    # def bed_line(self, fasta):
    #     row = self.to_dict(fasta)
    #     return f'{row["Chromosome"]}\t{row["Start"]}\t{row["End"]}\t{row["polyA_site"]}\t{row["total"]}\t{row["Strand"]}\t{row["fracA"]}\t{row["signal"]}\n'

    def to_dict(self, fasta):
        total = sum(self.counts)
        polyA = self.polyA_site()
        fracA = self.fraction_A(fasta, polyA)
        signal_seq_loc, signal_seq = self.polyA_signal_sequence(fasta, polyA)
        signal = f'{signal_seq_loc}@{signal_seq}'
        return {
            'Chromosome': self.Chromosome,
            'Start': self.Start,
            'End': self.End,
            'polyA_site': polyA,
            'count': total,
            'Strand': self.Strand,
            'fracA': fracA,
            'signal': signal
        }

    def __str__(self):
        return f'{self.Chromosome}:{self.Start}-{self.End}:{self.Strand}'

    def __repr__(self):
        return f'Cluster({str(self)})'


class TesClustering:

    def __init__(self, fasta):
        self.fasta = FastaStringExtractor(fasta, use_strand=True)

    def cluster(self, df_tes, extent_cutoff=3, window=25):
        for i, _df in tqdm(df_tes.groupby(['Chromosome', 'Strand'])):
            cluster = None

            for j, row in _df.sort_values('End').iterrows():

                if cluster is not None:
                    if (row['End'] - cluster.End) > window:
                        yield cluster
                        cluster = None

                # if enough reads supporting TES, create or extent cluster
                if row['count'] > extent_cutoff:
                    if cluster is None:
                        cluster = Cluster(row['Chromosome'], row['End'] - 1,
                                          row['End'], row['Strand'])
                    cluster.extend(row['End'], row['count'])

            if cluster is not None:
                yield cluster

    def to_bed(self, df_tes, bed_path):
        with open(bed_path, 'w') as f:
            for cluster in self.cluster(df_tes):
                f.write(cluster.bed_line(fasta))

    def to_df(self, df_tes):
        return pd.DataFrame([
            cluster.to_dict(self.fasta)
            for cluster in self.cluster(df_tes)
        ])


def tes_cluster_annotate(df_cluster, annotation):
    gr_apa = pr.PyRanges(df_cluster)

    gr_gtf = pr.read_gtf(annotation)
    greg = GenomicRegions(gr_gtf)
    gr = greg.annotate(gr_apa, single=True)
    # if genes already defined
    gr = greg.overlap_gene(gr)
    df = greg.overlap_tes(gr).df

    df.loc[df['Feature'].isna(), 'Feature'] = 'intergenic'

    df['canonical'] = df['canonical'] == 1
    df = df.rename(columns={'Start_b': 'canonical_site'})
    df['gene_id'] = df['gene_id'].replace('-1', '')
    df['gene_name'] = df['gene_name'].replace('-1', '')

    df['tpm'] = (df['count'] * 1000000 / df['count'].sum()).round(4)

    return df.sort_values(['Chromosome', 'End'])


def tes_sample(df_cluster, df_tes_sample,
               filter_intergenic=True,
               filter_internal_priming=True):
    columns = ['Chromosome', 'Start', 'End', 'Strand', 'gene_id']
    gr = pr.PyRanges(df_tes_sample)

    if filter_intergenic:
        df_cluster = df_cluster[df_cluster['Feature'] != 'intergenic']

    if filter_internal_priming:
        df_cluster = df_cluster[(df_cluster['fracA'] > 7)
                                & (df_cluster['signal'] == 'None@None')]

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
                     df_apa['count'].sum()).round(4)
    return df_apa.reset_index()


def lapa(alignment, fasta, annotation, chrom_sizes, output_dir, mapq=10):
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    print('Reading the alignment file...')
    if alignment.endswith('_read_annot.tsv'):
        df_alignment = read_talon_read_annot(alignment)
    elif alignment.endswith('.bam'):
        df_alignment = read_bam_ends(alignment, mapq=mapq)
    elif alignment.endswith('.csv'):
        df_alignment = read_sample_csv(alignment, mapq=10)
    else:
        raise ValueError(
            'Unknown file alignment format: supported '
            'file formats are `bam` and `sample.csv`')

    print('Counting TES...')
    df_tes, tes = count_tes(df_alignment, chrom_sizes, output_dir)
    del df_alignment

    print('Clustering TES and calculating polyA_sites...')
    df_cluster = TesClustering(fasta).to_df(df_tes)
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
