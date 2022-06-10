import logging
from pathlib import Path
import pandas as pd
import pyranges as pr
from lapa.cluster import PolyAClustering, TssClustering
from lapa.genomic_regions import PolyAGenomicRegions, TssGenomicRegions
from lapa.count import TesMultiCounter, TssMultiCounter
from lapa.utils.io import cluster_col_order, tss_cluster_col_order, \
    read_talon_read_annot
from lapa.replication import replication_dataset


class _Lapa:

    def __init__(self, fasta, annotation, chrom_sizes, output_dir,
                 method, mapq=10,
                 cluster_extent_cutoff=3, cluster_window=25,
                 cluster_ratio_cutoff=0.05,
                 min_replication_rate=0.95, replication_rolling_size=1000,
                 replication_num_sample=2, replication_min_count=1,
                 non_replicates_read_threhold=10):

        self.fasta = fasta
        self.annotation = annotation
        self.chrom_sizes = chrom_sizes
        self.output_dir = Path(output_dir)

        # parameters
        self.method = method
        self.mapq = mapq

        # clustering parameters
        self.cluster_extent_cutoff = cluster_extent_cutoff
        self.cluster_window = cluster_window
        self.cluster_ratio_cutoff = cluster_ratio_cutoff

        # replication parameters
        self.min_replication_rate = min_replication_rate
        self.replication_rolling_size = replication_rolling_size
        self.replication_num_sample = replication_num_sample
        self.replication_min_count = replication_min_count
        self.non_replicates_read_threhold = non_replicates_read_threhold

        # create file structure
        self.output_dir.mkdir(exist_ok=True)
        self.sample_dir.mkdir()
        self.raw_sample_dir.mkdir()
        self.dataset_dir.mkdir()

        self.is_read_annot = False

        # initilize logs directory and files
        self.make_logs()

    @property
    def sample_dir(self):
        return self.output_dir / 'sample'

    @property
    def raw_sample_dir(self):
        return self.output_dir / 'raw_sample'

    @property
    def dataset_dir(self):
        return self.output_dir / 'dataset'

    def make_logs(self):
        log_dir = self.output_dir / 'logs'
        log_dir.mkdir()

        # logs files
        self.progress_log = logging.getLogger('progress')
        self.progress_log.setLevel(logging.INFO)
        self.progress_log.addHandler(logging.FileHandler(
            log_dir / 'progress.log'))

        self.warn_log = logging.getLogger('warning')
        self.warn_log.setLevel(logging.WARN)
        self.warn_log.addHandler(logging.FileHandler(
            log_dir / 'warnings.log'))

        self.final_log = logging.getLogger('final')
        self.final_log.setLevel(logging.INFO)
        self.final_log.addHandler(logging.FileHandler(
            log_dir / 'final_stats.log'))

    def prepare_alignment(self, alignment):
        '''
        Get sample and dataset mapping and respective bam files
        from the given alignment file.
        '''
        if alignment.endswith('.csv'):
            df = pd.read_csv(alignment)

            if 'dataset' not in df.columns:
                df['dataset'] = 'all'

            assert all(pd.Series(['sample', 'path', 'dataset']).isin(df.columns)), \
                'provided csv file should be consist of columns `sample` and `path`'

            df_alignment = df[['sample', 'path', 'dataset']]

        elif alignment.endswith('.bam'):
            alignments = alignment.split(',')

            df_alignment = pd.DataFrame({
                'sample': [Path(i).stem for i in alignments],
                'dataset': ['all'] * len(alignments),
                'path': alignments
            })

        elif alignment.endswith('_read_annot.tsv'):
            # HACKED: to support talon read_annot
            self.is_read_annot = True
            self.method = None
            self.df_read_annot = read_talon_read_annot(alignment)
            df_alignment = None

        else:
            raise ValueError(
                'Unknown file alignment format: supported '
                'file formats are `bam` and `sample.csv`')

        if not self.is_read_annot:
            sample_dataset_mapping = df_alignment \
                .groupby('dataset')['sample'] \
                .agg(list).to_dict()
        else:
            sample_dataset_mapping = {
                'all': self.df_read_annot['sample'].unique().tolist()
            }

        return df_alignment, sample_dataset_mapping

    def counting(self, alignment):
        if not self.is_read_annot:
            counter = self.create_counter(alignment)
        else:
            counter = self.create_counter(
                self.df_read_annot, is_read_annot=True)

        df_all_count, sample_counts = counter.to_df()
        counter._to_bigwig(df_all_count, sample_counts, self.chrom_sizes,
                           self.output_dir, prefix=self.prefix)

        return df_all_count, sample_counts

    def clustering(self, df_counts):
        return self.create_clustering().to_df(df_counts)

    def annotate_cluster(self, df_cluster):
        gr = pr.PyRanges(df_cluster)

        total = gr.count.sum()
        df = self.create_genomic_regions().annotate(gr)

        df['tpm'] = df['count'] * 1000000 / total

        return df.sort_values(['Chromosome', 'End'])

    def save_cluster(self, df, path):
        return df[self.cluster_col_order].to_csv(
            path, index=False, sep='\t', header=False)

    def save_clusters(self, df_cluster, raw=False):
        filename = f'{self.prefix}_clusters.bed'
        if raw:
            filename = 'raw_' + filename

        self.save_cluster(df_cluster, self.output_dir / filename)

    def save_samples(self, samples, raw):
        if raw:
            sample_dir = self.raw_sample_dir
        else:
            sample_dir = self.sample_dir

        for sample, df in samples.items():
            self.save_cluster(
                df, sample_dir / f'{sample}.bed')

    def save_datasets(self, datasets):
        for dataset, df in datasets.items():
            self.save_cluster(
                df, self.dataset_dir / f'{dataset}.bed')

    def filter_replication(self, sample_clusters, sample_dataset_mapping):

        replicated_samples = dict()

        for _, samples in sample_dataset_mapping.items():

            rep_samples = {
                sample: sample_clusters[sample] for sample in samples
            }

            if len(samples) == 1:
                sample = samples[0]
                self.warn_log.warning(
                    f'Appling non_replicates_read_threhold={self.non_replicates_read_threhold}'
                    ' to filter sample {sample} because sample does not replicates')
                df_cluster = rep_samples[sample]
                rep_samples = {
                    sample: df_cluster[df_cluster['count'] >= self.non_replicates_read_threhold]
                }

            else:
                rep_samples = replication_dataset(
                    rep_samples, 'count',
                    self.replication_rolling_size, self.min_replication_rate,
                    self.replication_num_sample, self.replication_min_count)

            for sample, df in rep_samples.items():
                replicated_samples[sample] = self.calculate_usage(
                    df.drop(['tpm', 'gene_count', 'usage'], axis=1))

        return replicated_samples

    def sample_cluster(self, df_cluster, sample_counts):
        gr = pr.PyRanges(sample_counts)

        columns = ['Chromosome', 'Start', 'End', 'Strand', 'gene_id']

        gr_cluster = pr.PyRanges(df_cluster[columns].drop_duplicates())
        df_join = gr_cluster.join(gr, suffix='_sample').df

        df_join = df_join.groupby(columns, observed=True) \
                         .agg({'count': 'sum', 'coverage': 'mean'}) \
                         .reset_index()

        threshold = df_join['coverage'] * self.cluster_ratio_cutoff
        df_join = df_join[df_join['count'] > threshold]

        df_apa = df_join.reset_index().set_index(columns).join(
            df_cluster.set_index(columns)[self._keep_cols], rsuffix='_cluster'
        ).reset_index()

        return df_apa

    @staticmethod
    def calculate_usage(df_cluster):
        # get number of tes in the gene
        df_cluster = df_cluster.set_index('gene_id').join(
            df_cluster[['gene_id', 'count']]
            .groupby('gene_id').agg('sum').rename(
                columns={'count': 'gene_count'}))

        # calculate usage
        df_cluster['usage'] = df_cluster['count'] / df_cluster['gene_count']

        df_cluster['tpm'] = (df_cluster['count'] * 1000000 /
                             df_cluster['count'].sum()).round(2)

        return df_cluster.reset_index()

    @classmethod
    def aggregate_samples(cls, samples):

        aggregation = {'count': 'sum'}

        for i in cls._keep_cols:
            aggregation[i] = 'first'

        df = pd.concat(samples).groupby([
            'gene_id', 'Chromosome', 'Start', 'End', 'Strand'
        ], observed=True).agg(aggregation).reset_index()

        return cls.calculate_usage(df)

    def create_counter(self, alignment):
        raise NotImplementedError()

    def create_clustering(self):
        raise NotImplementedError()

    def create_genomic_regions(self):
        raise NotADirectoryError()

    def __call__(self, alignment):
        '''
        Main function run all the steps of LAPA
        '''
        # Create alignment dataframe from given input file
        alignment, sample_dataset_mapping = self.prepare_alignment(alignment)

        # Count read numbers
        self.progress_log.info('===== Counting reads (1 / 5) ===== \n')
        df_all_count, sample_counts = self.counting(alignment)

        # Create clusters from read count numbers
        self.progress_log.info(
            '===== Clustering and calling peaks (2 / 5) ===== \n')
        df_cluster = self.clustering(df_all_count)

        # Annotate clusters and save
        self.progress_log.info('===== Annotating cluster (3 / 5) ===== \n')
        df_cluster = self.annotate_cluster(df_cluster)
        df_cluster = self.calculate_usage(df_cluster)

        self.save_clusters(df_cluster, raw=True)

        # Get clusters of samples from all cluster and sample count
        self.progress_log.info(
            '===== Calculating and subseting clusters samples (4 / 5) ===== \n')

        sample_clusters = dict()
        for sample, df_counts in sample_counts.items():
            df_sample_cluster = self.sample_cluster(df_cluster, df_counts)
            df_sample_cluster = self.calculate_usage(df_sample_cluster)
            sample_clusters[sample] = df_sample_cluster

        self.save_samples(sample_clusters, raw=True)

        # Subset samples based on the replication rate
        self.progress_log.info('===== Calculating replication and producing '
                               'replicated clusters (5 / 5) ===== \n')
        replicated_sample_clusters = self.filter_replication(
            sample_clusters, sample_dataset_mapping)
        self.save_samples(replicated_sample_clusters, raw=False)

        # Aggreate to dataset level and save dataset clusters
        datasets = dict()
        for dataset, samples in sample_dataset_mapping.items():
            datasets[dataset] = self.aggregate_samples(
                [replicated_sample_clusters[i] for i in samples])
        self.save_datasets(datasets)

        # Aggreate and save all clusters
        df_all = self.aggregate_samples(
            replicated_sample_clusters.values())
        self.save_clusters(df_all)


class Lapa(_Lapa):

    _keep_cols = [
        'polyA_site', 'fracA', 'signal',
        'Feature', 'annotated_site'
    ]
    
    def __init__(self, fasta, annotation, chrom_sizes, output_dir, method='end',
                 min_tail_len=10, min_percent_a=0.9, mapq=10,
                 cluster_extent_cutoff=3, cluster_window=25, cluster_ratio_cutoff=0.05,
                 min_replication_rate=0.95, replication_rolling_size=1000,
                 filter_internal_priming=True, replication_num_sample=2,
                 replication_min_count=1, non_replicates_read_threhold=10):

        if method not in {'tail', 'end'}:
            raise ValueError(
                f'`method` need to be either `tail` or `end` but {method} given')

        super().__init__(fasta, annotation, chrom_sizes, output_dir,
                         method, mapq,
                         cluster_extent_cutoff, cluster_window,
                         cluster_ratio_cutoff,
                         min_replication_rate, replication_rolling_size,
                         replication_num_sample, replication_min_count)

        self.min_tail_len = min_tail_len
        self.min_percent_a = min_percent_a

        self.filter_internal_priming = filter_internal_priming

        self.prefix = 'polyA'
        self.cluster_col_order = cluster_col_order

    def create_counter(self, alignment, is_read_annot=False):
        return TesMultiCounter(alignment, self.method, self.mapq,
                               self.min_tail_len, self.min_percent_a,
                               is_read_annot=is_read_annot)

    def create_clustering(self):
        return PolyAClustering(self.fasta,
                               extent_cutoff=self.cluster_extent_cutoff,
                               ratio_cutoff=self.cluster_ratio_cutoff,
                               window=self.cluster_window)

    def create_genomic_regions(self):
        return PolyAGenomicRegions(self.annotation)

    def sample_cluster(self, df_cluster, sample_counts):

        if self.filter_internal_priming:
            df_cluster = df_cluster[~((df_cluster['fracA'] > 7)
                                      & (df_cluster['signal'] == 'None@None'))]

        return super().sample_cluster(df_cluster, sample_counts)


class LapaTss(_Lapa):

    _keep_cols = [
            'tss_site', 'Feature', 'annotated_site'
    ]
    
    def __init__(self, fasta, annotation, chrom_sizes, output_dir,
                 method='start', mapq=10,
                 cluster_extent_cutoff=3, cluster_window=25,
                 cluster_ratio_cutoff=0.05,
                 min_replication_rate=0.95, replication_rolling_size=1000,
                 replication_num_sample=2, replication_min_count=1,
                 non_replicates_read_threhold=10):

        if method not in {'start'}:
            raise ValueError(
                f'`method` need to be either `start` but {method} given')

        super().__init__(fasta, annotation, chrom_sizes, output_dir,
                         method, mapq,
                         cluster_extent_cutoff, cluster_window,
                         cluster_ratio_cutoff,
                         min_replication_rate, replication_rolling_size,
                         replication_num_sample, replication_min_count)

        self.prefix = 'tss'
        self.cluster_col_order = tss_cluster_col_order

    def create_counter(self, alignment, is_read_annot=False):
        return TssMultiCounter(alignment, self.method, self.mapq,
                               is_read_annot=is_read_annot)

    def create_clustering(self):
        return TssClustering(self.fasta,
                             extent_cutoff=self.cluster_extent_cutoff,
                             ratio_cutoff=self.cluster_ratio_cutoff,
                             window=self.cluster_window)

    def create_genomic_regions(self):
        return TssGenomicRegions(self.annotation)


def lapa(alignment: str, fasta: str, annotation: str, chrom_sizes :str, output_dir:str,
         method='end', min_tail_len=10, min_percent_a=0.9, mapq=10,
         cluster_extent_cutoff=3, cluster_window=25, cluster_ratio_cutoff=0.05,
         min_replication_rate=0.95, replication_rolling_size=1000,
         replication_num_sample=2, replication_min_count=1,
         non_replicates_read_threhold=10):
    '''
    LAPA high level api for polyA cluster calling.

    Args:
      alignment: Single or multiple bam file paths are separated
        with a comma.Alternatively, CSV file with columns of
        sample, dataset, path where the sample columns contains
        the name of the sample, the dataset is the group of
        samples replicates of each other, and path is the
        path of bam file.
      fasta: Genome reference (GENCODE or ENSEMBL fasta)
      annotation: Standart genome annotation (GENCODE or ENSEMBL gtf).
        GENCODE gtf file do not contains annotation for
        `five_prime_utr` and `three_prime_utr` so need to be
        corrected with `gencode_utr_fix`
      chrom_sizes: Chrom sizes files (can be generated with
      output_dir: See lapa.readthedocs.io/en/latest/output.html)
        for the details of the directory structure and file format.
      method: Counting method either `end` or `tail` where tails
        counting only counts reads with poly(A)-tail with certain
        length defined by `--min_tail_len` parameter. `end`
        counting still detects tails if exists but uses end
        location of all the reads in counting regardless of tail length.
      min_tail_len: Minimum tail length for `tail` counting strategy.
        This parameter will be ignored in `end` counting setting.
      min_percent_a: Minimum percentage of A bp in soft-trimmed segment
        to consider the segment as tails.
        This parameter will be ignored for end counting.
      mapq: Minimum read quality to required for tes calling
      cluster_extent_cutoff: Minimum number of reads to initialized cluster and
        terminated cluster will be terminated if read
        numbers below this cutoff for certain number of base pairs.
      cluster_ratio_cutoff: Percentage of coverage change for initialize cluster.
        At least x% of reads covering the bp need to ended in the
        position to initilized the cluster. This filter implies
        `<x%` of the reads given position could stop by chance
        so filtered as noise.
      cluster_window: Patience threshold to wait for termination cluster.
        If reads counts below the threshold for x bp then cluster
        will be terminated otherwise cluster will be extended.
        if number of reads subceed `the cluster_extent_cutoff.
      min_replication_rate: Minimum replication rate to include cluster in
        replicated clusters. 0.95 is recommended cutoff for
        experimental replication and 75% for biological replication.
      replication_rolling_size: Replication rolling size to calculate replication rate.
      replication_num_sample: Number of samples which region need to be observed for replication
      replication_min_count: Minimum count needed to recognize region as expressed
      non_replicates_read_threhold: Minimum read count need for the samples without replication.
        If there is not replicate samples for the sample, this default cutoff will be applied.
    '''
    _lapa = Lapa(fasta, annotation, chrom_sizes, output_dir, method=method,
                 min_tail_len=min_tail_len, min_percent_a=min_percent_a, mapq=mapq,
                 cluster_extent_cutoff=cluster_extent_cutoff,
                 cluster_window=cluster_window, cluster_ratio_cutoff=cluster_ratio_cutoff,
                 min_replication_rate=min_replication_rate,
                 replication_rolling_size=replication_rolling_size,
                 replication_num_sample=replication_num_sample,
                 replication_min_count=replication_min_count,
                 non_replicates_read_threhold=non_replicates_read_threhold)
    _lapa(alignment)


def lapa_tss(alignment: str, fasta: str, annotation: str, chrom_sizes: str, output_dir: str,
             method='start', mapq=10,
             cluster_extent_cutoff=3, cluster_window=25, cluster_ratio_cutoff=0.05,
             min_replication_rate=0.95, replication_rolling_size=1000,
             replication_num_sample=2, replication_min_count=1,
             non_replicates_read_threhold=10):
    '''
    LAPA TSS high level api for polyA cluster calling.

    Args:
      alignment: Single or multiple bam file paths are separated
        with a comma.Alternatively, CSV file with columns of
        sample, dataset, path where the sample columns contains
        the name of the sample, the dataset is the group of
        samples replicates of each other, and path is the
        path of bam file.
      fasta: Genome reference (GENCODE or ENSEMBL fasta)
      annotation: Standart genome annotation (GENCODE or ENSEMBL gtf).
        GENCODE gtf file do not contains annotation for
        `five_prime_utr` and `three_prime_utr` so need to be
        corrected with `gencode_utr_fix`
      chrom_sizes: Chrom sizes files (can be generated with
      method: Counting method
      output_dir: See lapa.readthedocs.io/en/latest/output.html)
        for the details of the directory structure and file format.
      min_tail_len: Minimum tail length for `tail` counting strategy.
        This parameter will be ignored in `end` counting setting.
      min_percent_a: Minimum percentage of A bp in soft-trimmed segment
        to consider the segment as tails.
        This parameter will be ignored for end counting.
      mapq: Minimum read quality to required for tes calling
      cluster_extent_cutoff: Minimum number of reads to initialized cluster and
        terminated cluster will be terminated if read
        numbers below this cutoff for certain number of base pairs.
      cluster_ratio_cutoff: Percentage of coverage change for initialize cluster.
        At least x% of reads covering the bp need to ended in the
        position to initilized the cluster. This filter implies
        `<x%` of the reads given position could stop by chance
        so filtered as noise.
      cluster_window: Patience threshold to wait for termination cluster.
        If reads counts below the threshold for x bp then cluster
        will be terminated otherwise cluster will be extended.
        if number of reads subceed `the cluster_extent_cutoff.
      min_replication_rate: Minimum replication rate to include cluster in
        replicated clusters. 0.95 is recommended cutoff for
        experimental replication and 75% for biological replication.
      replication_rolling_size: Replication rolling size to calculate replication rate.
      replication_num_sample: Number of samples which region need to be observed for replication
      replication_min_count: Minimum count needed to recognize region as expressed
      non_replicates_read_threhold: Minimum read count need for the samples without replication.
        If there is not replicate samples for the sample, this default cutoff will be applied.
    '''
    _lapa = LapaTss(fasta, annotation, chrom_sizes, output_dir,
                    method=method, mapq=mapq,
                    cluster_extent_cutoff=cluster_extent_cutoff,
                    cluster_window=cluster_window, cluster_ratio_cutoff=cluster_ratio_cutoff,
                    min_replication_rate=min_replication_rate,
                    replication_rolling_size=replication_rolling_size,
                    replication_num_sample=replication_num_sample,
                    replication_min_count=replication_min_count,
                    non_replicates_read_threhold=non_replicates_read_threhold)
    _lapa(alignment)
