from pathlib import Path
import numpy as np
import pandas as pd
import pyranges as pr
from tqdm import tqdm
from kipoiseq import Interval
from kipoiseq.extractors import FastaStringExtractor
from longread_postprocessing.utils.io import read_talon_read_annot, \
    bw_from_pyranges
from longread_postprocessing.utils.common import pad_series, polyA_signal_seqs
from longread_postprocessing.genomic_regions import GenomicRegions


def count_tes(df_alignment, chrom_sizes, output_dir):
    # count TES sites
    columns = ['Chromosome', 'End', 'Strand', 'gene_id']
    df = df_alignment[columns]
    df['count'] = 1
    df = df.groupby(columns).agg('sum').reset_index()

    df['Start'] = df['End'] - 1

    bw_from_pyranges(
        pr.PyRanges(df), 'count',
        chrom_sizes,
        str(output_dir / 'tes_counts_pos.bw'),
        str(output_dir / 'tes_counts_neg.bw')
    )
    del df['Start']

    df = df.rename(columns={'End': 'tes_loc'})

    # set polyA site
    df['tes_site'] = df['Chromosome'].astype(str) + ':' + \
        df['tes_loc'].astype(str) + ':' + df['Strand']

    return df.set_index('tes_site')


class Cluster:

    def __init__(self, Chromosome: str, Start: int, End: int, Strand: str, gene_id: str, counts=None):
        self.Chromosome = str(Chromosome)
        self.Start = Start
        self.End = End
        self.Strand = Strand
        self.gene_id = gene_id
        self.counts = counts or list()

    def extend(self, end, count):
        self.End = end
        self.counts.append(count)

    def __len__(self):
        return self.End - self.Start

    def polyA_site(self, window=5, std=1):
        pad_size = window // 2
        counts = pad_series(self.counts, pad_size=pad_size)
        moving_sum = counts.rolling(window, center=True,
                                    win_type='gaussian').sum(std=std)
        return self.Start + moving_sum.idxmax()

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

        for singal_seq in polyA_signal_seqs:
            index = seq.find(singal_seq)
            if index == -1:
                continue
            else:
                return self.Start + index, singal_seq
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

    def bed_line(self, fasta):
        row = self.to_dict(fasta)
        return f'{row["Chromosome"]}\t{row["Start"]}\t{row["End"]}\t{row["polyA_site"]}\t{row["total"]}\t{row["Strand"]}\t{row["fracA"]}\t{row["signal"]}\n'

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
            'singal': signal
        }


class TesCluster:

    def __init__(self, fasta):
        self.fasta = FastaStringExtractor(fasta, use_strand=True)

    def cluster(self, df_tes, extent_cutoff=3, window=25):
        for i, _df in tqdm(df_tes.groupby('gene_id')):
            cluster = None

            for j, row in _df.iterrows():
                # if enough reads supporting TES, create or extent cluster
                if row['count'] > 3:
                    if cluster is None:
                        cluster = Cluster(row['Chromosome'], row['tes_loc'],
                                          row['tes_loc'], row['Strand'],
                                          row['gene_id'])
                    cluster.extend(row['tes_loc'], row['count'])

                if cluster is not None:
                    if (row['tes_loc'] - cluster.End) > 25:
                        yield cluster
                        cluster = None
            if cluster:
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

    gr_tes = gr_gtf.features.tes()
    df_tes = gr_tes.df[['Chromosome', 'Start',
                        'End', 'Strand']].drop_duplicates()
    df_tes['canonical'] = True

    gr_tes = pr.PyRanges(df_tes, int64=True)
    gr = gr.join(gr_tes, how='left')

    df = gr.df

    del df['End_b']
    del df['Strand_b']

    df['canonical'] = df['canonical'] == 1
    df = df.rename(columns={'Start_b': 'canonical_site'})
    # df['canonical_site'] = df['canonical_site'].replace(-1, np.nan)
    df['tpm'] = df['count'] * 1000000 / df['count'].sum()

    return df


def longread_postprocessing_tes(alignment, fasta, annotation, chrom_sizes, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    print('Reading the alignment file...')
    if alignment.endswith('_read_annot.tsv'):
        df_alignment = read_talon_read_annot(alignment)
    elif alignment.endswith('.bam') or alignment.endswith('.sam'):
        raise ValueError('Bam file is not supported at the moment')
    else:
        raise ValueError('Alignment file need to be sam or bam file.')

    print('Counting TES...')
    df_tes = count_tes(df_alignment, chrom_sizes, output_dir)
    # del df_alignment

    print('Clustering TES and calculating polyA_sites...')
    df_cluster = TesCluster(fasta).to_df(df_tes)
    # del df_tes

    print('Annotating TES cluster...')
    df_cluster = tes_cluster_annotate(df_cluster, annotation)

    df_cluster.to_csv(output_dir / 'polyA_clusters.bed',
                      index=False, sep='\t', header=False)
