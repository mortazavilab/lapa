from collections import defaultdict, Counter
import pysam
import pandas as pd
import pyranges as pr
from tqdm import tqdm
import matplotlib.pyplot as plt
from lapa.utils.io import bw_from_pyranges, read_sample_csv


class BaseTesCounter:

    def __init__(self, bam_file, mapq=10):
        self.mapq = mapq

        if type(bam_file) is pysam.AlignmentFile:
            self.bam_file = bam_file.filename
        elif type(bam_file) is str:
            self.bam_file = bam_file
        else:
            raise ValueError('Type of bam_file need to be either'
                             ' `str` (file path) or `pysam.AlignmentFile`')

    @property
    def bam(self):
        return pysam.AlignmentFile(self.bam_file, 'rb')

    def count(self):
        raise NotImplementedError()

    def to_df(self):
        rows = [
            (chrom, site, strand, count)
            for (chrom, site, strand), count in self.count().items()
        ]

        df = pd.DataFrame(
            rows, columns=['Chromosome', 'End', 'Strand', 'count'])
        df['Start'] = df['End'] - 1
        return df[['Chromosome', 'Start', 'End', 'Strand', 'count']]


class TailTesCounter(BaseTesCounter):

    def __init__(self, bam_file, min_tail_len=10, min_percent_a=0.9, mapq=10):
        self.min_tail_len = min_tail_len
        self.min_percent_a = min_percent_a
        super().__init__(bam_file, mapq)

    @staticmethod
    def detect_polyA_tail(read: pysam.AlignedSegment):
        """Detect polyA tails from a read

        Args:
          read: aligned reads

        Returns:
          Tuple of polyA_site, length of tail, percent of A base in tails.
        """
        cigartuples = read.cigartuples

        # check if clipped and get clip
        if len(cigartuples) <= 1:
            return None, 0, 0

        if read.is_reverse:
            clip = cigartuples[0]
        else:
            clip = cigartuples[-1]

        if clip[0] != 4:
            return None, 0, 0

        # calculate info about polyA
        clip_len = clip[1]

        if read.is_reverse:
            tail_seq = read.seq[:clip_len][::-1]
            polyA_site = read.reference_start
            tail_base = 'T'
        else:
            tail_seq = read.seq[-clip_len:]
            polyA_site = read.reference_end
            tail_base = 'A'

        tail_len, percent_a = TailTesCounter._calculate_tail_seq(
            tail_seq, tail_base)
        return polyA_site, tail_len, percent_a

    @staticmethod
    def _calculate_tail_seq(tail_seq, tail_base):
        score = 0
        best_score = 0

        tail_len = 0
        num_a = 0
        _num_a = 0

        for i, base in enumerate(tail_seq):
            if base == tail_base:
                score += 1
                _num_a += 1
            else:
                score -= 4

            if score >= best_score:
                best_score = score
                tail_len = i + 1
                num_a = _num_a

        percent_a = num_a / tail_len if tail_len > 0 else 0
        return tail_len, percent_a

    def iter_tailed_reads(self):
        """Iterates polyA reads and polyA_site based on polyA filters.
        """
        for read in tqdm(self.bam):

            if read.is_secondary:
                continue

            if read.mapping_quality < self.mapq:
                continue

            polyA_site, tail_len, percent_a = self.detect_polyA_tail(read)

            if (tail_len >= self.min_tail_len) \
               and (percent_a >= self.min_percent_a):
                yield read, polyA_site, tail_len, percent_a

    def save_tailed_reads(self, output_bam):
        tailed_bam = pysam.AlignmentFile(output_bam, "wb", template=self.bam)

        for read, polyA_site, tail_len, percent_a in self.iter_tailed_reads():
            tailed_bam.write(read)

        tailed_bam.close()

    def tail_len_dist(self):
        tail_dist = defaultdict(int)

        for _, polyA_site, tail_len, percent_a in self.iter_tailed_reads():
            tail_dist[tail_len] += 1

        return pd.Series(tail_dist).sort_index()

    def plot_tail_len_dist(self):
        dist = self.tail_len_dist()
        dist.plot()
        dist.cumsum().plot()
        plt.legend([
            'pdf',
            'cdf'
        ])
        plt.xlabel('Length of polyA tail')
        plt.ylabel('Number of reads')

    def count(self):
        tes = Counter()

        for read, polyA_site, tail_len, percent_a in self.iter_tailed_reads():
            strand = '-' if read.is_reverse else '+'
            tes[(read.reference_name, polyA_site, strand)] += 1

        return tes


class EndTesCounter(BaseTesCounter):

    def count(self):
        tes = Counter()

        for read in tqdm(self.bam):
            if read.mapping_quality >= self.mapq:
                if read.is_reverse:
                    strand = '-'
                    polyA_site = read.reference_start
                else:
                    strand = '+'
                    polyA_site = read.reference_end
                tes[(read.reference_name, polyA_site, strand)] += 1

        return tes


def save_tes_count_bw(df, output_dir, chrom_sizes, sample):
    bw_from_pyranges(
        pr.PyRanges(df), 'count',
        chrom_sizes,
        str(output_dir / f'{sample}_tes_counts_pos.bw'),
        str(output_dir / f'{sample}_tes_counts_neg.bw')
    )


def count_tes_bam_samples(df_alignment, method, min_tail_len=10,
                          min_percent_a=0.9, mapq=10):
    assert method in {'tail', 'end'}, \
        '`method` parameter for counting need to either `tail` or `end`'

    df = list()

    for _, row in df_alignment.iterrows():
        if method == 'tail':
            _df = TailTesCounter(row['path'], min_tail_len,
                                 min_percent_a, mapq).to_df()
        elif method == 'end':
            _df = EndTesCounter(row['path'], mapq).to_df()
        else:
            raise ValueError('`method` need to be either `tail` or `end`')

        _df['sample'] = row['sample']
        df.append(_df)

    return pd.concat(df)


def agg_tes_samples(df, chrom_sizes, output_dir):
    columns = ['Chromosome', 'Start', 'End', 'Strand']
    df_all = df.groupby(columns).agg('sum').reset_index()

    tes = dict()
    samples = df['sample'].unique()

    if len(samples) > 1:
        for sample, _df in df.groupby('sample'):
            _df = _df.groupby(columns).agg('sum').reset_index()
            tes[sample] = _df
            save_tes_count_bw(_df, output_dir, chrom_sizes, sample)
    else:
        tes[samples[0]] = df_all

    save_tes_count_bw(df_all, output_dir, chrom_sizes, 'all')
    return df_all, tes
