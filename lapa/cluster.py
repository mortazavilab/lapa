import logging
import pandas as pd
from tqdm import tqdm
from kipoiseq import Interval
from kipoiseq.extractors import FastaStringExtractor
from lapa.utils.common import pad_series, polyA_signal_seqs


def _tqdm_clustering(iterable):
    '''
    Adaptor for tqdm to integrate to logging
    '''
    logger = logging.getLogger('progress')
    file = logger.handlers[0].stream if logger.handlers else None

    return tqdm(iterable, mininterval=5, file=file,
                bar_format='- {n_fmt} Chromosome/Strand clustered...\n')


class Cluster:

    def __init__(self, Chromosome: str, Start: int, End: int,
                 Strand: str, counts=None, fields=None):
        self.Chromosome = str(Chromosome)
        self.Start = Start
        self.End = End
        self.Strand = Strand
        self.counts = counts or list()
        self.fields = fields or dict()

    def extend(self, end, count):
        self.End = end
        self.counts.append((end, count))

    @property
    def total_count(self):
        return sum(i for _, i in self.counts)

    def __len__(self):
        return self.End - self.Start

    @staticmethod
    def _count_arr(counts):
        counts = sorted(counts)
        min_pos = counts[0][0]
        max_pos = counts[-1][0]
        count_arr = [0] * (max_pos - min_pos + 1)
        for pos, c in counts:
            count_arr[pos - min_pos] = c

        return pd.Series(count_arr,
                         index=pd.RangeIndex(min_pos, max_pos + 1))

    def peak(self, window=5, std=1):

        pad_size = window // 2
        counts = self._count_arr(self.counts)

        # counts = pd.Series(dict(self.counts)).sort_index()
        counts = pad_series(counts, pad_size=pad_size)
        moving_sum = counts.rolling(window, center=True,
                                    win_type='gaussian').sum(std=std)
        return moving_sum.idxmax()

    def to_dict(self, fasta):
        total = self.total_count
        return {
            'Chromosome': self.Chromosome,
            'Start': self.Start,
            'End': self.End,
            'peak': self.peak(),
            'count': total,
            'Strand': self.Strand,
        }

    def __str__(self):
        return f'{self.Chromosome}:{self.Start}-{self.End}:{self.Strand}'

    def __repr__(self):
        return f'Cluster({str(self)})'


class PolyACluster(Cluster):

    def polyA_site(self, window=5, std=1):
        return self.peak(window=window, std=std)

    def polyA_signal_sequence(self, fasta, polyA_site):
        # the list of poly(A) signals
        # that are annotated for single 3' end processing sites
        # in the region -60 to +10 nt around them
        # According to Gruber et al., 2016, Genome Research.

        if isinstance(fasta, str):
            fasta = FastaStringExtractor(fasta, use_strand=True)

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
                pos = (start + index) if self.Strand == '+' \
                    else (end - index - 6)
                return pos, signal_seq

        return None, None

    def fraction_A(self, fasta, polyA_site):

        if isinstance(fasta, str):
            fasta = FastaStringExtractor(fasta, use_strand=True)

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
        return f'{row["Chromosome"]}\t{row["Start"]}\t{row["End"]}\t{row["polyA_site"]}\t{row["count"]}\t{row["Strand"]}\t{row["fracA"]}\t{row["signal"]}\n'

    def to_dict(self, fasta):
        cluster = super().to_dict(fasta)
        cluster['polyA_site'] = cluster['peak']
        del cluster['peak']

        cluster['fracA'] = self.fraction_A(fasta, cluster['polyA_site'])
        signal_seq_loc, signal_seq = self.polyA_signal_sequence(
            fasta, cluster['polyA_site'])
        cluster['signal'] = f'{signal_seq_loc}@{signal_seq}'

        return cluster


class TssCluster(Cluster):

    def to_dict(self, fasta):
        cluster = super().to_dict(fasta)
        cluster['tss_site'] = cluster['peak']
        del cluster['peak']
        return cluster


class Clustering:

    Cluster = Cluster

    def __init__(self, fasta, extent_cutoff=3, ratio_cutoff=0.05, window=25,
                 groupby=None, fields=None, progress=True):
        self.fasta = FastaStringExtractor(fasta, use_strand=True)
        self.extent_cutoff = extent_cutoff
        self.ratio_cutoff = ratio_cutoff
        self.window = window
        self.groupby = groupby or ['Chromosome', 'Strand']
        self.fields = fields or list()
        self.progress = progress

    def cluster(self, df_tes):
        _groupby = df_tes.groupby(self.groupby)
        if self.progress:
            _groupby = _tqdm_clustering(_groupby)

        for _, _df in _groupby:
            cluster = None

            for _, row in _df.sort_values('End').iterrows():

                if cluster is not None:
                    if (row['End'] - cluster.End) > self.window:
                        yield cluster
                        cluster = None

                # if enough reads supporting TES, create or extent cluster
                threshold = max(self.extent_cutoff,
                                row['coverage'] * self.ratio_cutoff)
                if (row['count'] >= threshold):
                    if cluster is None:
                        start = max(row['End'] - 1, 0)  # avoid -1
                        cluster = self.Cluster(row['Chromosome'], start,
                                               row['End'], row['Strand'])
                        cluster.fields = {i: list() for i in self.fields}

                    cluster.extend(row['End'], row['count'])
                    for i in self.fields:
                        cluster.fields[i].append(row[i])

            if cluster is not None:
                yield cluster

    def to_bed(self, df_tes, bed_path):
        with open(bed_path, 'w') as f:
            for cluster in self.cluster(df_tes):
                f.write(cluster.bed_line(self.fasta))

    def to_df(self, df_tes):
        return pd.DataFrame([
            cluster.to_dict(self.fasta)
            for cluster in self.cluster(df_tes)
        ])


class PolyAClustering(Clustering):

    Cluster = PolyACluster


class TssClustering(Clustering):

    Cluster = TssCluster
