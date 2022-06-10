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
    '''
    Cluster class representing cluster in genome in a chromosome
      with start, end location and strand. Mainly, each position
      between start and end location has counts as value indicating
      number of reads ending at the location. Based on the counts
      peak calling is performed. Counts are stored sparse format but
      converted to dense format for peak calling.

    Args:
      Chromosome: Chromosome of the cluster
      Start: Chromosome of the cluster
      End: Chromosome of the cluster
      Strand: Chromosome of the cluster
      counts (:obj:`List[Tuple[int, int]]`, optional): Sparse representation of
        counts in the cluster as list of tuple where tuple is (position, count).
      fields (:obj:`Dict[int, List]]`, optional): Fields are dictonary of list
        where information about each read in the cluster can be stored.
        Not used in clustering algorith by default.

    Examples:
      Peak calling based on the counts.

      >>> cluster = Cluster('chr1', 10, 11, '+')
      >>> cluster.extend((12, 5))
      >>> len(cluster)
      3
      >>> cluster.peak()
      12
      >>> cluster.extend((15, 5))
      >>> cluster.extend((16, 3))
      >>> cluster.extend((18, 1))
      >>> len(cluster)
      8
      >>> cluster.peak()
      15
    '''

    def __init__(self, Chromosome: str, Start: int, End: int,
                 Strand: str, counts=None, fields=None):
        self.Chromosome = str(Chromosome)
        self.Start = Start
        self.End = End
        self.Strand = Strand
        self.counts = counts or list()
        self.fields = fields or dict()

    def extend(self, end, count):
        '''
        Extend cluster with end position and count.

        Args:
          end: New end location to extend cluster.
          count: Number of reads supporting position.
        '''
        self.End = end
        self.counts.append((end, count))

    @property
    def total_count(self):
        '''
        Total counts in the cluster.
        '''
        return sum(i for _, i in self.counts)

    def __len__(self):
        '''
        Length of the cluster `End - Start`.
        '''
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
        '''
        Detech peaks position where value are maximum in the cluster.
          Counts are smoothed with moving average before performing
          peak calling.

        Args:
          window: window size for smoothing.
          std: Standard deviation of gaussian kernel used for smoothing.
        '''
        pad_size = window // 2
        counts = self._count_arr(self.counts)

        # counts = pd.Series(dict(self.counts)).sort_index()
        counts = pad_series(counts, pad_size=pad_size)
        moving_sum = counts.rolling(window, center=True,
                                    win_type='gaussian').sum(std=std)
        return moving_sum.idxmax()

    def to_dict(self, fasta: FastaStringExtractor):
        '''
        Convert cluster into dictonary annotate regulatory elements
          in the cluster using fasta file.

        Args:
          fasta: FastaStringExtractor object of kipoiseq.
        '''
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
    '''
    PolyA cluster used in poly(A)-site clustering, performs peak calling to
      obtain exact poly(A)-site, and extract sequence elements in
      the vicinity of cluster.

    Examples:
      Peak calling for poly(A)-site detection.

      >>> cluster = PolyACluster('chr1', 100, 101, '+')
      >>> cluster.extend((112, 5))
      >>> cluster.extend((115, 5))
      >>> cluster.extend((116, 3))
      >>> cluster.extend((118, 1))
      >>> len(cluster)
      8
      >>> cluster.polyA_site()
      115
      >>> cluster.polyA_signal_sequence('hg38.fa', polyA_site=115)
      118, 'AATAAA'
      >>> cluster.fraction_A('hg38.fa', polyA_site=115)
      3
    '''

    def polyA_site(self, window=5, std=1):
        '''
        Detects poly(A)-site with peak calling.
          Counts are smoothed with moving average before performing
          peak calling.

        Args:
          window: window size for smoothing.
          std: Standard deviation of gaussian kernel used for smoothing.
        '''
        return self.peak(window=window, std=std)

    def polyA_signal_sequence(self, fasta: FastaStringExtractor,
                              polyA_site: int):
        '''
        Poly(A) signal sequence in the vicinity of poly(A) site.

        Args:
          fasta: Fasta to extract sequences.
          polyA_site: Poly(A) site based on the peak calling.
        '''
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

    def fraction_A(self, fasta: FastaStringExtractor, polyA_site):
        '''
        Fraction of A following the polyA site.

        Args:
          fasta: Fasta to extract sequences.
          polyA_site: Poly(A) site based on the peak calling.
        '''
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

    def to_dict(self, fasta):
        '''
        Convert cluster into dictonary annotate regulatory elements
          (polyA_signal and fracA) in the cluster using fasta file.

        Args:
          fasta: FastaStringExtractor object of kipoiseq.
        '''
        cluster = super().to_dict(fasta)
        cluster['polyA_site'] = cluster['peak']
        del cluster['peak']

        cluster['fracA'] = self.fraction_A(fasta, cluster['polyA_site'])
        signal_seq_loc, signal_seq = self.polyA_signal_sequence(
            fasta, cluster['polyA_site'])
        cluster['signal'] = f'{signal_seq_loc}@{signal_seq}'

        return cluster


class TssCluster(Cluster):
    '''
    TSS cluster used in tss site clustering, performs peak calling to
      obtain exact tss-site, and extract sequence elements in
      the vicinity of cluster.

    Examples:
      Peak calling for tss site detection.

      >>> cluster = PolyACluster('chr1', 100, 101, '+')
      >>> cluster.extend((112, 5))
      >>> cluster.extend((115, 5))
      >>> cluster.extend((116, 3))
      >>> cluster.extend((118, 1))
      >>> len(cluster)
      8
      >>> cluster.peak()
      115
    '''

    def to_dict(self, fasta):
        '''
        Convert cluster into dictonary.

        Args:
          fasta: FastaStringExtractor object of kipoiseq.
        '''
        cluster = super().to_dict(fasta)
        cluster['tss_site'] = cluster['peak']
        del cluster['peak']
        return cluster


class Clustering:
    '''
    Clustering algorith to obtains regions cluster together based
      on the read end counts.

    Args:
      fasta: path to fasta file which used to extract regulatory elements
        in the vicinity of genome.
      extent_cutoff: Extent cluster if number of read end counts
        above this cutoff.
      ratio_cutoff: Ratio of read end counts to coverage of the region.
      window: Patiance window cluster will be terminated on if read numbers
       below this cutoff for this window size of bps.
      groupby: Groupby reads in the same region and sceen read number
        default of Chromosome and Strand.
      fields: Fields to extract from the counts and store in cluster object.
      progress: Show progress bar for clustering
    '''

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
        '''

        Args:
          df_tes: .
        '''
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

    def to_df(self, df_tes):
        '''
        Perform clustering based on read end counts.

        Args:
          df_tes: Counts per genomics position obtain with counting classes
            in pandas.DataFrame with
            `Chromosome, Start, End, Strand, count, coverage` columns.
        '''
        return pd.DataFrame([
            cluster.to_dict(self.fasta)
            for cluster in self.cluster(df_tes)
        ])


class PolyAClustering(Clustering):
    '''
    Clustering algorith to obtains polyA clusters from read end counts.

    Examples:
      Cluster poly(A)-sites from bam file.

      >>> clustering = PolyAClustering('hg38.fasta')
      >>> counter = ThreePrimeCounter(bam_file)
      >>> df_counts = counter.to_df()
      >>> df_counts.head()
      +--------------+-----------+-----------+--------------+-----------+------------+
      | Chromosome   | Start     | End       | Strand       | count     | coverage   |
      | (category)   | (int32)   | (int32)   | (category)   | (int64)   | (int64)    |
      |--------------+-----------+-----------+--------------+-----------+------------|
      | chr1         | 887771    | 887772    | +            | 5         | 5          |
      | chr1         | 994684    | 994685    | -            | 8         | 10         |
      ...
      >>> df_clusters = clustering.to_df(df_counts)
      >>> df_clusters
      +--------------+-----------+-----------+--------------+-----------+--------------+-----------+---------------+
      | Chromosome   |     Start |       End |   polyA_site |     count | Strand       |     fracA | signal        |
      | (category)   |   (int32) |   (int32) |      (int64) |   (int64) | (category)   |   (int64) | (object)      |
      |--------------+-----------+-----------+--------------+-----------+--------------+-----------+---------------|
      | chr17        |    100099 |    100100 |       100100 |        10 | +            |         6 | 100098@GATAAA |
      | chr17        |    100199 |    100200 |       100200 |         7 | -            |         2 | None@None     |
      | chrM         |      1100 |      1101 |         1101 |        11 | +            |        -1 | None@None     |
      ...
    '''
    Cluster = PolyACluster


class TssClustering(Clustering):
    '''
    Clustering algorith to obtains tss clusters from read start counts.

    Examples:
      Cluster tss sites from bam file.

      >>> clustering = PolyAClustering('hg38.fasta')
      >>> counter = FivePrimeCounter(bam_file)
      >>> df_counts = counter.to_df()
      >>> df_counts.head()
      +--------------+-----------+-----------+--------------+-----------+------------+
      | Chromosome   | Start     | End       | Strand       | count     | coverage   |
      | (category)   | (int32)   | (int32)   | (category)   | (int64)   | (int64)    |
      |--------------+-----------+-----------+--------------+-----------+------------|
      | chr1         | 887771    | 887772    | +            | 5         | 5          |
      | chr1         | 994684    | 994685    | -            | 8         | 10         |
      ...
      >>> df_clusters = clustering.to_df(df_counts)
      >>> df_clusters
      +--------------+-----------+-----------+--------------+-----------+--------------+
      | Chromosome   |     Start |       End |     tss_site |     count | Strand       |
      | (category)   |   (int32) |   (int32) |      (int64) |   (int64) | (category)   |
      |--------------+-----------+-----------+--------------+-----------+--------------|
      | chr17        |    100099 |    100100 |       100100 |        10 | +            |
      | chr17        |    100199 |    100200 |       100200 |         7 | -            |
      | chrM         |      1100 |      1101 |         1101 |        11 | +            |
      ...
    '''
    Cluster = TssCluster
