from collections import defaultdict, Counter
import pysam
import pandas as pd
import pyranges as pr
from tqdm import tqdm
import matplotlib.pyplot as plt
from lapa.utils.io import bw_from_pyranges, read_sample_csv, \
    read_talon_read_annot_five_prime_count, \
    read_talon_read_annot_three_prime_count


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
        +--------------+-----------+-----------+--------------+------------+
        | Chromosome   | Start     | End       | Strand       | count      |
        | (category)   | (int32)   | (int32)   | (category)   | (uint16)   |
        |--------------+-----------+-----------+--------------+------------|
        | chr1         | 887771    | 887772    | +            | 5          |
        | chr1         | 994684    | 994685    | -            | 8          |
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

    def iter_reads(self):
        bam = self.bam

        if self.progress:
            bam = tqdm(bam)

        for read in bam:
            if self.filter_read(read):
                yield read

    def filter_read(self, read):
        '''
        Filter reads from counting if not true
        '''
        return (not read.is_secondary) and (read.mapping_quality >= self.mapq)

    def count(self):
        '''
        Counts reads per positions defined by count_Read_function

        Returns:
            Dict[(chrom, pos, strand), int]: Returns dictionary of
                chromosome, position and strand ad index and counts as values.
        '''
        counts = Counter()

        for read in self.iter_reads():
            site = self.count_read(read)
            strand = '-' if read.is_reverse else '+'
            counts[(read.reference_name, site, strand)] += 1

        return dict(counts)

    def count_read(self, read: pysam.AlignedSegment):
        raise NotImplementedError()

    def to_df(self):
        '''
        Counts as dataframe with columns of
            `['Chromosome', 'End', 'Strand', 'count']`
        '''
        df = pd.DataFrame([
            (chrom, site, strand, count)
            for (chrom, site, strand), count in self.count().items()
        ], columns=['Chromosome', 'End', 'Strand', 'count'])

        df['Start'] = df['End'] - 1
        return df

    @staticmethod
    def _to_bigwig(df, chrom_sizes, output_dir, prefix):
        bw_from_pyranges(
            pr.PyRanges(df), 'count',
            chrom_sizes,
            str(output_dir / f'{prefix}_pos.bw'),
            str(output_dir / f'{prefix}_neg.bw')
        )

    def to_bigwig(self, chrom_sizes, output_dir, prefix='lapa_counts'):
        '''
        Saves counts as bigwig file for each strand

        Args:
            chrom_sizes (str): Chrom sizes files (can be generated with) from fasta
                with `faidx fasta -i chromsizes > chrom_sizes`
            output_dir: Output directory to save bigwig files
            prefix (str): File prefix to used in bigwig the files
        '''
        self._to_bigwig(self.to_df(), chrom_sizes, output_dir, prefix)


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

        >>> counter = EndCounter(bam_file)
        >>> counter.to_bigwig(chrom_sizes, output_dir, 'mid')
        >>> os.listdir(output_dir)
        ['lapa_count_pos.bw', 'lapa_count_neg.bw']
        >>> counter.to_df()
        +--------------+-----------+-----------+--------------+------------+
        | Chromosome   | Start     | End       | Strand       | count      |
        | (category)   | (int32)   | (int32)   | (category)   | (uint16)   |
        |--------------+-----------+-----------+--------------+------------|
        | chr1         | 887771    | 887772    | +            | 5          |
        | chr1         | 994684    | 994685    | -            | 8          |
        ...
    '''

    def count_read(self, read: pysam.AlignedSegment):
        '''
        Returns 3' end of the read
        '''
        if read.is_reverse:
            return read.reference_start
        else:
            return read.reference_end


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
        +--------------+-----------+-----------+--------------+------------+
        | Chromosome   | Start     | End       | Strand       | count      |
        | (category)   | (int32)   | (int32)   | (category)   | (uint16)   |
        |--------------+-----------+-----------+--------------+------------|
        | chr1         | 887771    | 887772    | +            | 5          |
        | chr1         | 994684    | 994685    | -            | 8          |
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
        +--------------+-----------+-----------+--------------+------------+
        | Chromosome   | Start     | End       | Strand       | count      |
        | (category)   | (int32)   | (int32)   | (category)   | (uint16)   |
        |--------------+-----------+-----------+--------------+------------|
        | chr1         | 887771    | 887772    | +            | 5          |
        | chr1         | 994684    | 994685    | -            | 8          |
        ...
    '''

    def __init__(self, bam_file, mapq=10, progress=True,
                 min_tail_len=10, min_percent_a=0.9):
        super().__init__(bam_file, mapq, progress)
        self.min_tail_len = min_tail_len
        self.min_percent_a = min_percent_a

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

        tail_len, percent_a = PolyaTailCounter._calculate_tail_seq(
            tail_seq, tail_base)
        return polyA_site, tail_len, percent_a

    @staticmethod
    def _calculate_tail_seq(tail_seq, tail_base):
        '''Calculate tail seq'''
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

    def _read_is_tailed(self, tail_len, percent_a):
        return (tail_len >= self.min_tail_len) \
            and (percent_a >= self.min_percent_a)

    def iter_tailed_reads(self):
        """Iterates polyA reads and polyA_site based on polyA filters.
        """
        for read in super().iter_reads():
            polyA_site, tail_len, percent_a = self.detect_polyA_tail(read)

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
        dist.plot()
        dist.cumsum().plot()
        plt.legend([
            'pdf',
            'cdf'
        ])
        plt.xlabel('Length of polyA tail')
        plt.ylabel('Number of reads')

    def filter_read(self, read):
        '''Filter tailed reads and quality'''
        _, tail_len, percent_a = self.detect_polyA_tail(read)
        return super().filter_read(read) \
            and self._read_is_tailed(tail_len, percent_a)


class BaseMultiCounter:

    def __init__(self, alignment: str, method: str, mapq=10):
        self.alignment = alignment
        self.method = method
        self.mapq = mapq

        self.is_read_annot = False        
        self.df_alignment = self._prepare_alignment(str(alignment))

    def build_counter(self, bam):
        raise NotImplementedError()

    def _count_read_annot(self):
        raise NotADirectoryError()

    def _prepare_alignment(self, alignment):
        '''
        '''
        if alignment.endswith('.csv'):
            df = pd.read_csv(alignment)
            assert all(pd.Series(['sample', 'path']).isin(df.columns)), \
                'provided csv file should be consist of columns `sample` and `path`'
            return df[['sample', 'path']]

        elif alignment.endswith('.bam'):
            alignments = alignment.split(',')
            return pd.DataFrame({
                'sample': ['all'] * len(alignments),
                'path': alignments
            })

        elif alignment.endswith('_read_annot.tsv'):
            # HACKED: to support talon read_annot
            self.is_read_annot = True
            self.method = None
            return None

        else:
            raise ValueError(
                'Unknown file alignment format: supported '
                'file formats are `bam` and `sample.csv`')

    @staticmethod
    def _to_bigwig(df_all, tes, chrom_sizes,
                   output_dir, prefix='lapa_counts'):
        save_count_bw(df_all, output_dir, chrom_sizes, f'all_{prefix}')

        if len(tes) > 1:
            for sample, df in tes.items():
                save_count_bw(df, output_dir, chrom_sizes,
                              f'{sample}_{prefix}')

    def to_df(self):
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

    def __init__(self, alignment, method='end', mapq=10,
                 min_tail_len=10, min_percent_a=0.9):
        self.min_tail_len = min_tail_len
        self.min_percent_a = min_percent_a
        super().__init__(alignment, method, mapq)

    def build_counter(self, bam):
        if self.method == 'tail':
            return PolyaTailCounter(bam, self.mapq, self.min_tail_len,
                                    self.min_percent_a)
        elif self.method == 'end':
            return ThreePrimeCounter(bam, self.mapq)
        else:
            raise ValueError(f'`method` need to be either `tail` or `end` but {method} given')

    def _count_read_annot(self):
        return read_talon_read_annot_three_prime_count(self.alignment)


class TssMultiCounter(BaseMultiCounter):

    def __init__(self, alignment, method='start', mapq=10):
        super().__init__(alignment, method, mapq)

    def build_counter(self, bam):
        if self.method == 'start':
            return FivePrimeCounter(bam, self.mapq)
        else:
            raise ValueError('`method` need to be either `start`')

    def _count_read_annot(self):
        return read_talon_read_annot_five_prime_count(self.alignment)
