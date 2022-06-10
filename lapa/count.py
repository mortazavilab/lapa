import logging
from collections import defaultdict, Counter
import pysam
import pandas as pd
import pyranges as pr
from tqdm import tqdm
import matplotlib.pyplot as plt
from lapa.utils.io import bw_from_pyranges, \
    _read_talon_read_annot_five_prime_count, \
    _read_talon_read_annot_three_prime_count


def _tqdm_counting(iterable):
    '''
    Adaptor for tqdm to integrate to logging
    '''
    logger = logging.getLogger('progress')
    file = logger.handlers[0].stream if logger.handlers else None

    return tqdm(iterable, mininterval=5, file=file,
                bar_format='- {n_fmt} reads counted...\n')


class BaseCounter:
    '''
    Base class to count features from alignment file.

    Args:
        bam_file: Path to bam file or `pysam.AlignmentFile` object.
        mapq (:obj:`int`, optional): minimum reads quality
            required to use in counting.
        progress (:obj:`bool`, optional): Show progress in counting.

    Attributes:
        bam (pysam.AlignmentFile): alignment file.

    Examples:
        Count mid location of the read.

        >>> class MidCounter(BaseCounter):
        >>>     def count_read(self, read):
        >>>         return (read.reference_start + read.reference_end) / 2
        >>> counter = MidCounter(bam_file)
        >>> counter.to_bigwig(chrom_sizes, output_dir, 'mid')
        >>> os.listdir(output_dir)
        ['mid_pos.bw', 'mid_neg.bw']
        >>> counter.to_df()
        +--------------+-----------+-----------+--------------+-----------+------------+
        | Chromosome   | Start     | End       | Strand       | count     | coverage   |
        | (category)   | (int32)   | (int32)   | (category)   | (int64)   | (int64)    |
        |--------------+-----------+-----------+--------------+-----------+------------|
        | chr1         | 887771    | 887772    | +            | 5         | 5          |
        | chr1         | 994684    | 994685    | -            | 8         | 10         |
        ...
    '''

    def __init__(self, bam_file, mapq=10, progress=True):
        self.mapq = mapq
        self.progress = progress

        if type(bam_file) is pysam.AlignmentFile:
            self.bam_file = bam_file.filename
        elif type(bam_file) is str:
            self.bam_file = bam_file
        else:
            raise ValueError('Type of bam_file need to be either'
                             ' `str` (file path) or `pysam.AlignmentFile`')

    @property
    def bam(self):
        '''Alignment file used in counting.'''
        return pysam.AlignmentFile(self.bam_file, 'rb')

    def iter_reads(self, chrom=None, strand=None):
        if chrom:
            bam = self.bam.fetch(chrom)
        else:
            bam = self.bam

        if self.progress:
            bam = _tqdm_counting(bam)

        for read in bam:
            if self.filter_read(read):
                if strand:
                    read_strand = '-' if read.is_reverse else '+'
                    if read_strand == strand:
                        yield read
                else:
                    yield read

    def filter_read(self, read):
        '''
        Filter reads from counting if not true
        '''
        return (read.flag in {0, 16}) and (read.mapping_quality >= self.mapq)

    def count(self):
        '''
        Counts reads per positions defined by count_Read_function

        Returns:
            Dict[(chrom, pos, strand), int]: Returns dictionary of
                chromosome, position and strand ad index and counts as values.
        '''
        logging.getLogger('progress').info(f'Counting `{self.bam_file}`:')

        counts = Counter()

        for read in self.iter_reads():
            site = self.count_read(read)
            strand = '-' if read.is_reverse else '+'
            counts[(read.reference_name, site, strand)] += 1

        return dict(counts)

    def count_read(self, read: pysam.AlignedSegment):
        raise NotImplementedError()

    def to_gr(self):
        '''
        Counts as dataframe with columns of
            `['Chromosome', 'Start', 'End', 'Strand', 'count']`
        '''
        df = pd.DataFrame([
            (chrom, site, strand, count)
            for (chrom, site, strand), count in self.count().items()
        ], columns=['Chromosome', 'End', 'Strand', 'count'])
        df['Start'] = df['End'] - 1
        df_bam = pr.read_bam(self.bam_file, mapq=self.mapq, as_df=True)
        df_bam['Start'] -= 1

        gr_bam = pr.PyRanges(df_bam)
        gr_bam = gr_bam[gr_bam.Flag.isin({0, 16})]

        return pr.PyRanges(df).count_overlaps(
            gr_bam,
            overlap_col='coverage',
            strandedness='same')

    def to_df(self):
        return self.to_gr().df.astype({'Chromosome': 'str', 'Strand': 'str'})

    @staticmethod
    def _to_bigwig(gr, chrom_sizes, output_dir, prefix):
        if isinstance(gr, pd.DataFrame):
            gr = pr.PyRanges(gr)

        (output_dir / 'counts').mkdir(exist_ok=True)
        (output_dir / 'coverage').mkdir(exist_ok=True)
        (output_dir / 'ratio').mkdir(exist_ok=True)

        bw_from_pyranges(
            gr, 'count',
            chrom_sizes,
            str(output_dir / 'counts' / f'{prefix}_counts_pos.bw'),
            str(output_dir / 'counts' / f'{prefix}_counts_neg.bw')
        )
        bw_from_pyranges(
            gr, 'coverage',
            chrom_sizes,
            str(output_dir / 'coverage' / f'{prefix}_coverage_pos.bw'),
            str(output_dir / 'coverage' / f'{prefix}_coverage_neg.bw')
        )
        gr = gr.assign('ratio', lambda df: df['count'] / df['coverage'])
        bw_from_pyranges(
            gr, 'ratio',
            chrom_sizes,
            str(output_dir / 'ratio' / f'{prefix}_ratio_pos.bw'),
            str(output_dir / 'ratio' / f'{prefix}_ratio_neg.bw')
        )

    def to_bigwig(self, chrom_sizes, output_dir, prefix='lapa_counts'):
        '''
        Saves counts as bigwig file for each strand

        Args:
            chrom_sizes (str): Chrom sizes files (can be generated with) from
                fasta with `faidx fasta -i chromsizes > chrom_sizes`
            output_dir: Output directory to save bigwig files
            prefix (str): File prefix to used in bigwig the files
        '''
        self._to_bigwig(self.to_gr(), chrom_sizes, output_dir, prefix)


def save_count_bw(df, output_dir, chrom_sizes, prefix):
    '''
    Saves counts as bigwig file for each strand generated with
        instance of BaseCounter object

    Args:
        chrom_sizes (str): Chrom sizes files (can be generated with) from fasta
            with `faidx fasta -i chromsizes > chrom_sizes`
        output_dir: Output directory to save bigwig files
        prefix (str): File prefix to used in bigwig the files
    '''
    BaseCounter._to_bigwig(df, chrom_sizes, output_dir, prefix)


def save_tss_count_bw(df, chrom_sizes, output_dir, prefix):
    '''
    Saves counts of transcript end sites (TSS) as bigwig file
        for each strand generated with TesCounters

    Args:
        chrom_sizes (str): Chrom sizes files (can be generated with) from fasta
            with `faidx fasta -i chromsizes > chrom_sizes`
        output_dir: Output directory to save bigwig files
        prefix (str): File prefix to used in bigwig the files
    '''
    save_count_bw(df, output_dir, chrom_sizes, f'{prefix}_tss')


def save_tes_count_bw(df, chrom_sizes, output_dir, prefix):
    '''
    Saves counts of transcript end sites (TSS) as bigwig file
        for each strand generated with TesCounters

    Args:
        chrom_sizes (str): Chrom sizes files (can be generated with) from fasta
            with `faidx fasta -i chromsizes > chrom_sizes`
        output_dir: Output directory to save bigwig files
        prefix (str): File prefix to used in bigwig the files
    '''
    save_count_bw(df, output_dir, chrom_sizes, f'{prefix}_tes')


class ThreePrimeCounter(BaseCounter):
    '''
    Counts 3' ends of reads (transcript end sites) per position
        from alignment file.

    Args:
        bam_file: Path to bam file or `pysam.AlignmentFile` object.
        mapq (:obj:`int`, optional): minimum reads quality
            required to use in counting.
        progress (:obj:`bool`, optional): Show progress in counting.

    Attributes:
        bam (pysam.AlignmentFile): alignment file.

    Examples:
        Count 3' ends of the read per position.

        >>> counter = ThreePrimeCounter(bam_file)
        >>> counter.to_bigwig(chrom_sizes, output_dir, 'mid')
        >>> os.listdir(output_dir)
        ['lapa_count_pos.bw', 'lapa_count_neg.bw']
        >>> counter.to_df()
        +--------------+-----------+-----------+--------------+-----------+------------+
        | Chromosome   | Start     | End       | Strand       | count     | coverage   |
        | (category)   | (int32)   | (int32)   | (category)   | (int64)   | (int64)    |
        |--------------+-----------+-----------+--------------+-----------+------------|
        | chr1         | 887771    | 887772    | +            | 5         | 5          |
        | chr1         | 994684    | 994685    | -            | 8         | 10         |
        ...
    '''

    def count_read(self, read: pysam.AlignedSegment):
        '''
        Returns 3' end of the read
        '''
        if read.is_reverse:
            polya_site = read.reference_start
            tail_base = 'T'
        else:
            polya_site = read.reference_end
            tail_base = 'A'

        cigartuples = read.cigartuples

        # check if clipped and get clip
        if len(cigartuples) <= 1:
            return polya_site

        if read.is_reverse:
            clip = cigartuples[0]
        else:
            clip = cigartuples[-1]

        if clip[0] != 4:  # soft-clip
            return polya_site

        # calculate info about polyA
        clip_len = clip[1]

        if read.is_reverse:
            read_seq = read.seq[clip_len:]
            tail_base = 'T'
        else:
            read_seq = read.seq[:-clip_len][::-1]
            tail_base = 'A'

        read_tail_len, read_percent_a = PolyaTailCounter._calculate_tail_seq(
            read_seq, tail_base)
        polya_site += read_tail_len if read.is_reverse else -read_tail_len

        return polya_site

    @staticmethod
    def _calculate_tail_seq(tail_seq, tail_base):
        '''Calculate tail seq'''

        assert tail_base in {'A', 'T'}

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
                score -= 5

            if score >= best_score:
                best_score = score
                tail_len = i + 1
                num_a = _num_a

            if score < -20:
                break

        percent_a = num_a / tail_len if tail_len > 0 else 0
        return tail_len, percent_a


class FivePrimeCounter(BaseCounter):
    '''
    Counts 5' ends of reads (transcript end sites) per position
        from alignment file.

    Args:
        bam_file: Path to bam file or `pysam.AlignmentFile` object.
        mapq (:obj:`int`, optional): minimum reads quality
            required to use in counting.
        progress (:obj:`bool`, optional): Show progress in counting.

    Attributes:
        bam (pysam.AlignmentFile): alignment file.

    Examples:
        Count 5' ends of the read per position.

        >>> counter = StartCounter(bam_file)
        >>> counter.to_bigwig(chrom_sizes, output_dir, 'mid')
        >>> os.listdir(output_dir)
        ['lapa_count_pos.bw', 'lapa_count_neg.bw']
        >>> counter.to_df()
        +--------------+-----------+-----------+--------------+-----------+------------+
        | Chromosome   | Start     | End       | Strand       | count     | coverage   |
        | (category)   | (int32)   | (int32)   | (category)   | (int64)   | (int64)    |
        |--------------+-----------+-----------+--------------+-----------+------------|
        | chr1         | 887771    | 887772    | +            | 5         | 5          |
        | chr1         | 994684    | 994685    | -            | 8         | 10         |
        ...
    '''

    def count_read(self, read: pysam.AlignedSegment):
        '''
        Returns 5' end of the read
        '''
        if read.is_reverse:
            return read.reference_end
        else:
            return read.reference_start


class PolyaTailCounter(ThreePrimeCounter):
    '''
    Counts 3' end of reads with polyA-tail (transcript end sites) per position
        from alignment file.

    Args:
        bam_file: Path to bam file or `pysam.AlignmentFile` object.
        mapq (:obj:`int`, optional): minimum reads quality
            required to use in counting.
        progress (:obj:`bool`, optional): Show progress in counting.
        min_tail_len:

    Attributes:
        bam (pysam.AlignmentFile): alignment file.

    Examples:
        Count 3' ends of the read per position.

        >>> counter = PolyaTailCounter(bam_file)
        >>> counter.to_bigwig(chrom_sizes, output_dir, 'mid')
        >>> os.listdir(output_dir)
        ['lapa_count_pos.bw', 'lapa_count_neg.bw']
        >>> counter.to_df()
        +--------------+-----------+-----------+--------------+-----------+------------+
        | Chromosome   | Start     | End       | Strand       | count     | coverage   |
        | (category)   | (int32)   | (int32)   | (category)   | (int64)   | (int64)    |
        |--------------+-----------+-----------+--------------+-----------+------------|
        | chr1         | 887771    | 887772    | +            | 5         | 5          |
        | chr1         | 994684    | 994685    | -            | 8         | 10         |
        ...
    '''

    def __init__(self, bam_file, mapq=10, progress=True,
                 min_tail_len=10, min_percent_a=0.9, count_aligned=False):
        super().__init__(bam_file, mapq, progress)
        self.min_tail_len = min_tail_len
        self.min_percent_a = min_percent_a
        self.count_aligned = count_aligned

    @staticmethod
    def detect_polyA_tail(read: pysam.AlignedSegment, count_aligned=False):
        """Detect polyA tails from a read

        Args:
          read: aligned reads
          count_aligned: Count aligned base pairs (likely internal priming)
            as well in tail length.

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
            read_seq = read.seq[clip_len:]
            tail_seq = read.seq[:clip_len][::-1]
            polyA_site = read.reference_start
            tail_base = 'T'
        else:
            read_seq = read.seq[:-clip_len][::-1]
            tail_seq = read.seq[-clip_len:]
            polyA_site = read.reference_end
            tail_base = 'A'

        tail_len, percent_a = PolyaTailCounter._calculate_tail_seq(
            tail_seq, tail_base)
        read_tail_len, read_percent_a = PolyaTailCounter._calculate_tail_seq(
            read_seq, tail_base)

        polyA_site += read_tail_len if read.is_reverse else -read_tail_len

        if count_aligned:
            total_len = tail_len + read_tail_len
            total_percent = (percent_a * tail_len +
                             read_percent_a * read_tail_len) / total_len

            return polyA_site, total_len, total_percent
        else:
            return polyA_site, tail_len, percent_a

    def _read_is_tailed(self, tail_len, percent_a):
        return (tail_len >= self.min_tail_len) \
            and (percent_a >= self.min_percent_a)

    def iter_tailed_reads(self):
        """Iterates polyA reads and polyA_site based on polyA filters.
        """
        for read in super().iter_reads():
            polyA_site, tail_len, percent_a = self.detect_polyA_tail(
                read, self.count_aligned)

            if self._read_is_tailed(tail_len, percent_a):
                yield read, polyA_site, tail_len, percent_a

    def save_tailed_reads(self, output_bam):
        '''
        Save tailed reads as bam files

        Args:
            output_bam: Path to bam file or `pysam.AlignmentFile` object.
        '''
        tailed_bam = pysam.AlignmentFile(output_bam, "wb", template=self.bam)

        for read, _, _, _ in self.iter_tailed_reads():
            tailed_bam.write(read)

        tailed_bam.close()

    def tail_len_dist(self):
        '''Returns tail length distribution of reads based on the filters.'''
        tail_dist = defaultdict(int)

        for _, polyA_site, tail_len, percent_a in self.iter_tailed_reads():
            tail_dist[tail_len] += 1

        return pd.Series(tail_dist).sort_index()

    def plot_tail_len_dist(self):
        '''Plots pdf and cdf of tail length distribution'''
        dist = self.tail_len_dist()
        df = pd.DataFrame({
            'count': dist,
            'cumulative count': dist.cumsum()
        })
        df = df[df['cumulative count'] / df['count'].sum() < 0.99]
        df['count'].plot()
        df['cumulative count'].plot()
        plt.legend([
            'count',
            'cumulative count'
        ])
        plt.xlabel('Length of polyA tail')
        plt.ylabel('Number of reads')

    def filter_read(self, read):
        '''Filter tailed reads and quality'''
        if super().filter_read(read):
            _, tail_len, percent_a = self.detect_polyA_tail(read)
            return self._read_is_tailed(tail_len, percent_a)
        return False


class BaseMultiCounter:
    '''
    Base class to counts reads from multiple aligment files.

    Args:
      df_alignment: DataFrame with columns of ['sample', 'dataset', 'path']
        where sample is the sample name, dataset is name of the group
        (replicates) of sample belong, path is the path to bam file.
      method: Counting method implemented by child class.
      mapq: minimum mapping quality
      is_read_annot: Talon reads annotate file can be provided to `df_alignment`
        argument in that case this argument need to True.
    '''

    def __init__(self, df_alignment: pd.DataFrame, method: str, mapq=10,
                 is_read_annot=False):
        self.df_alignment = df_alignment
        self.method = method
        self.mapq = mapq

        self.is_read_annot = is_read_annot

    def build_counter(self, bam):
        raise NotImplementedError()

    def _count_read_annot(self):
        raise NotImplementedError()

    @staticmethod
    def _to_bigwig(df_all, tes, chrom_sizes,
                   output_dir, prefix='polyA'):
        save_count_bw(df_all, output_dir, chrom_sizes, f'all_{prefix}')

        if len(tes) > 1:
            for sample, df in tes.items():
                save_count_bw(df, output_dir, chrom_sizes,
                              f'{sample}_{prefix}')

    def to_df(self):
        '''
        Export counts as dataframe.

        Returns:
          (pd.DataFrame, Dict[str, pd.DataFrame]): Counst as tuple
          the first element is dataframe of all the counts and second element
          dictonary where first element is the name of sample and second element
          dataframe of counts.
        '''
        # Counting from alignment files
        if self.is_read_annot:
            # HACKED: to support talon read_annot
            df = self._count_read_annot()
        else:
            df = pd.concat([
                self.build_counter(row['path'])
                .to_df()
                .assign(sample=row['sample'])
                for _, row in self.df_alignment.iterrows()
            ])
        # Aggreate counts per samples and all samples
        cols = ['Chromosome', 'Start', 'End', 'Strand']

        tes = dict()
        samples = df['sample'].unique()
        df_all = df.groupby(cols).agg('sum').reset_index()

        if len(samples) > 1:
            for sample, _df in df.groupby('sample'):
                _df = _df.groupby(cols).agg('sum').reset_index()
                tes[sample] = _df
        else:
            tes[samples[0]] = df_all

        return df_all, tes


class TesMultiCounter(BaseMultiCounter):
    '''
    Counts transcript end sites from multiple aligment files.

    Args:
        df_alignment: DataFrame with columns of ['sample', 'dataset', 'path']
          where sample is the sample name, dataset is name of the group
          (replicates) of sample belong, path is the path to bam file.
        method: either `end` or `tail` see `PolyaTailCounter` \
          and `ThreePrimeCounter` for countering behavior.
        mapq: minimum mapping quality
        is_read_annot: Talon reads annotate file can be provided to `df_alignment`
          argument in that case this argument need to True.

    Examples:
        Counts transcript end files for two samples with two replicates

        >>> df_alignment = pd.DataFrame({
        >>>     'sample': ['s1', 's2', 's3', 's4'],
        >>>     'dataset': ['d1', 'd2', 'd3', 'd4'],
        >>>     'path': ['s1.bam', 's2.bam', 's3.bam', 's4.bam']
        >>> })
        >>> counter = TesMultiCounter(df_alignment)
        >>> counter.to_bigwig(chrom_sizes, output_dir) # export counts as bw
        >>> df_all, samples = counter.to_df() # or export as df
        >>> df_all
        +--------------+-----------+-----------+--------------+-----------+------------+
        | Chromosome   | Start     | End       | Strand       | count     | coverage   |
        | (category)   | (int32)   | (int32)   | (category)   | (int64)   | (int64)    |
        |--------------+-----------+-----------+--------------+-----------+------------|
        | chr1         | 887771    | 887772    | +            | 5         | 5          |
        | chr1         | 994684    | 994685    | -            | 8         | 10         |
        ...
        >>> samples['s1']
        +--------------+-----------+-----------+--------------+-----------+------------+
        | Chromosome   | Start     | End       | Strand       | count     | coverage   |
        | (category)   | (int32)   | (int32)   | (category)   | (int64)   | (int64)    |
        |--------------+-----------+-----------+--------------+-----------+------------|
        | chr1         | 887771    | 887772    | +            | 5         | 5          |
        | chr1         | 994684    | 994685    | -            | 8         | 10         |
        ...
    '''

    def __init__(self, alignment, method='end', mapq=10,
                 min_tail_len=10, min_percent_a=0.9, is_read_annot=False):
        self.min_tail_len = min_tail_len
        self.min_percent_a = min_percent_a
        super().__init__(alignment, method, mapq, is_read_annot)

    def build_counter(self, bam):
        if self.method == 'tail':
            return PolyaTailCounter(bam, self.mapq, self.min_tail_len,
                                    self.min_percent_a)
        elif self.method == 'end':
            return ThreePrimeCounter(bam, self.mapq)
        else:
            raise ValueError(
                f'`method` need to be either `tail` or `end` but {method} given')

    def _count_read_annot(self):
        df = self.df_alignment

        cols = ['Chromosome', 'Start', 'End', 'Strand', 'sample']
        df_count = list()

        for sample, _df in df.groupby('sample'):

            _df_count = _read_talon_read_annot_three_prime_count(_df.copy()) \
                .groupby(cols).agg('sum').reset_index()
            _df['Start'] -= 1

            _df_count = pr.PyRanges(_df_count).count_overlaps(
                pr.PyRanges(_df),
                overlap_col='coverage',
                strandedness='same').df.astype({
                    'Chromosome': 'str', 'Strand': 'str'})

            df_count.append(_df_count)

        return pd.concat(df_count)


class TssMultiCounter(BaseMultiCounter):
    '''
    Counts transcript start sites from multiple aligment files.

    Args:
        df_alignment: DataFrame with columns of ['sample', 'dataset', 'path']
          where sample is the sample name, dataset is name of the group
          (replicates) of sample belong, path is the path to bam file.
        method: either `end` or `tail` see `FiveTailCounter`
        mapq: minimum mapping quality
        is_read_annot: Talon reads annotate file can be provided to
          `df_alignment` argument in that case this argument need to True.

    Examples:
        Counts transcript end files for two samples with two replicates

        >>> df_alignment = pd.DataFrame({
        >>>     'sample': ['s1', 's2', 's3', 's4'],
        >>>     'dataset': ['d1', 'd2', 'd3', 'd4'],
        >>>     'path': ['s1.bam', 's2.bam', 's3.bam', 's4.bam']
        >>> })
        >>> counter = TssMultiCounter(df_alignment)
        >>> counter.to_bigwig(chrom_sizes, output_dir) # export counts as bw
        >>> df_all, samples = counter.to_df() # or export as df
        >>> df_all
        +--------------+-----------+-----------+--------------+-----------+------------+
        | Chromosome   | Start     | End       | Strand       | count     | coverage   |
        | (category)   | (int32)   | (int32)   | (category)   | (int64)   | (int64)    |
        |--------------+-----------+-----------+--------------+-----------+------------|
        | chr1         | 887771    | 887772    | +            | 5         | 5          |
        | chr1         | 994684    | 994685    | -            | 8         | 10         |
        ...
        >>> samples['s1']
        +--------------+-----------+-----------+--------------+-----------+------------+
        | Chromosome   | Start     | End       | Strand       | count     | coverage   |
        | (category)   | (int32)   | (int32)   | (category)   | (int64)   | (int64)    |
        |--------------+-----------+-----------+--------------+-----------+------------|
        | chr1         | 887771    | 887772    | +            | 5         | 5          |
        | chr1         | 994684    | 994685    | -            | 8         | 10         |
        ...
    '''

    def __init__(self, alignment, method='start', mapq=10,
                 is_read_annot=False):
        super().__init__(alignment, method, mapq, is_read_annot)

    def build_counter(self, bam):
        if self.method == 'start':
            return FivePrimeCounter(bam, self.mapq)
        else:
            raise ValueError('`method` need to be either `start`')

    def _count_read_annot(self):
        df = self.df_alignment

        cols = ['Chromosome', 'Start', 'End', 'Strand', 'sample']
        df_count = list()

        for sample, _df in df.groupby('sample'):

            _df_count = _read_talon_read_annot_five_prime_count(_df.copy()) \
                .groupby(cols).agg('sum').reset_index()
            _df['Start'] -= 1

            _df_count = pr.PyRanges(_df_count).count_overlaps(
                pr.PyRanges(_df),
                overlap_col='coverage',
                strandedness='same').df.astype({
                    'Chromosome': 'str', 'Strand': 'str'})

            df_count.append(_df_count)

        return pd.concat(df_count)
