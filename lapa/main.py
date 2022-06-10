import click
from lapa.lapa import lapa, lapa_tss
from lapa.link import link_tss_to_tes
from lapa.correction import correct_talon


@click.command()
@click.option('--alignment',
              help='Single or multiple bam file paths are separated '
              'with a comma.Alternatively, CSV file with columns of '
              'sample, dataset, path where the sample columns contains '
              'the name of the sample, the dataset is the group of '
              'samples replicates of each other, and path is the '
              'path of bam file.',
              required=True)
@click.option('--fasta',
              help='Genome reference (GENCODE or ENSEMBL fasta)',
              required=True)
@click.option('--annotation',
              help='Standart genome annotation (GENCODE or ENSEMBL gtf). '
              'GENCODE gtf file do not contains annotation for '
              '`five_prime_utr` and `three_prime_utr` so need to be '
              'corrected with `gencode_utr_fix` '
              '(see https://github.com/MuhammedHasan/gencode_utr_fix.git).',
              required=True)
@click.option('--chrom_sizes',
              help='Chrom sizes files (can be generated with '
              '`faidx fasta -i chromsizes > chrom_sizes`)',
              required=True)
@click.option('--output_dir',
              help='Output directory of LAPA. '
              'See lapa.readthedocs.io/en/latest/output.html) '
              'for the details of the directory structure and file format.',
              required=True)
@click.option('--counting_method',
              help='Counting method either `end` or `tail` where tails '
              'counting only counts reads with poly(A)-tail with certain '
              'length defined by `--min_tail_len` parameter. `end` '
              'counting still detects tails if exists but uses end '
              'location of all the reads in counting regardless of tail length.',
              type=click.Choice(['end', 'tail'], case_sensitive=True),
              default='end')
@click.option('--min_tail_len',
              help='Minimum tail length for `tail` counting strategy.'
              'This parameter will be ignored in `end` counting setting.',
              default=10,
              type=int)
@click.option('--min_percent_a',
              help='Minimum percentage of A bp in soft-trimmed segment'
              ' to consider the segment as tails. '
              'This parameter will be ignored for end counting.',
              default=0.9,
              type=float)
@click.option('--mapq',
              help='Minimum read quality to required for tes calling.',
              default=10,
              type=int)
@click.option('--cluster_extent_cutoff',
              help='Minimum number of reads to initialized cluster and '
              'terminated cluster will be terminated if read'
              ' numbers below this cutoff for certain number of base pairs.',
              default=3,
              type=int)
@click.option('--cluster_ratio_cutoff',
              help='Percentage of coverage change for initialize cluster.'
              'At least x% of reads covering the bp need to ended in the '
              'position to initilized the cluster. This filter implies '
              '`<x%` of the reads given position could stop by chance '
              'so filtered as noise.',
              default=0.05,
              type=click.FloatRange(0, 1))
@click.option('--cluster_window',
              help='Patience threshold to wait for termination cluster.'
              'If reads counts below the threshold for x bp then cluster '
              'will be terminated otherwise cluster will be extended. '
              'if number of reads subceed `the cluster_extent_cutoff`.',
              default=25,
              type=int)
@click.option('--min_replication_rate',
              help='Minimum replication rate to include cluster in '
              'replicated clusters. 0.95 is recommended cutoff for '
              'experimental replication and 75% for biological replication.',
              default=0.95,
              type=click.FloatRange(0, 1))
@click.option('--replication_rolling_size',
              help='Replication rolling size to calculate replication rate.',
              default=1000,
              type=int)
@click.option('--replication_num_sample',
              help='Number of samples which region need to be observed for replication.',
              default=2,
              type=int)
@click.option('--replication_min_count',
              help='Minimum count needed to recognize region as expressed.',
              default=1,
              type=int)
@click.option('--non_replicates_read_threhold',
              help='Minimum read count need for the samples without replication.'
              ' If there is not replicate samples for the sample, '
              ' this default cutoff will be applied.',
              default=10,
              type=int)
def cli_lapa(alignment, fasta, annotation, chrom_sizes, output_dir,
             counting_method, min_tail_len=10, min_percent_a=0.9, mapq=10,
             cluster_extent_cutoff=3, cluster_ratio_cutoff=0.05, cluster_window=25,
             min_replication_rate=0.95, replication_rolling_size=1000,
             replication_num_sample=2, replication_min_count=1,
             non_replicates_read_threhold=10):
    '''
    CLI interface for lapa polyA cluster calling.
    '''
    lapa(alignment, fasta, annotation, chrom_sizes, output_dir,
         counting_method, min_tail_len=min_tail_len,
         min_percent_a=min_percent_a,
         cluster_extent_cutoff=cluster_extent_cutoff,
         cluster_ratio_cutoff=cluster_ratio_cutoff,
         cluster_window=cluster_window, mapq=mapq,
         min_replication_rate=min_replication_rate,
         replication_rolling_size=replication_rolling_size,
         replication_num_sample=replication_num_sample,
         replication_min_count=replication_min_count,
         non_replicates_read_threhold=non_replicates_read_threhold)


@click.command()
@click.option('--alignment',
              help='Single or multiple bam file paths are separated '
              'with a comma.Alternatively, CSV file with columns of '
              'sample, dataset, path where the sample columns contains '
              'the name of the sample, the dataset is the group of '
              'samples replicates of each other, and path is the '
              'path of bam file.',
              required=True)
@click.option('--fasta',
              help='Genome reference (GENCODE or ENSEMBL fasta)',
              required=True)
@click.option('--annotation',
              help='Standart genome annotation (GENCODE or ENSEMBL gtf). '
              'GENCODE gtf file do not contains annotation for '
              '`five_prime_utr` and `three_prime_utr` so need to be '
              'corrected with `gencode_utr_fix` '
              '(see https://github.com/MuhammedHasan/gencode_utr_fix.git)',
              required=True)
@click.option('--chrom_sizes',
              help='Chrom sizes files (can be generated with)'
              '`faidx fasta -i chromsizes > chrom_sizes`)',
              required=True)
@click.option('--output_dir',
              help='Output directory of LAPA. '
              'See lapa.readthedocs.io/en/latest/output.html) '
              'for the details of the directory structure and file format.',
              required=True)
@click.option('--mapq',
              help='Minimum read quality to required for tss calling',
              default=10,
              type=int)
@click.option('--cluster_extent_cutoff',
              help='Minimum number of reads to initialized cluster and '
              'terminated cluster will be terminated if read'
              ' numbers below this cutoff for certain number of base pairs.',
              default=3,
              type=int)
@click.option('--cluster_ratio_cutoff',
              help='Percentage of coverage change for initialize cluster.'
              'At least x% of reads covering the bp need to ended in the '
              'position to initilized the cluster. This filter implies '
              '`<x%` of the reads given position could stop by chance '
              'so ignored as noise ',
              default=0.05,
              type=click.FloatRange(0, 1))
@click.option('--cluster_window',
              help='Patience threshold to wait for termination cluster.'
              'If reads counts below the threshold for x bp then cluster '
              'will be terminated otherwise cluster will be extended. '
              'if number of reads subceed `the cluster_extent_cutoff`',
              default=25,
              type=int)
@click.option('--min_replication_rate',
              help='Minimum replication rate to include cluster in '
              'replicated clusters. 0.95 is recommended cutoff for '
              'experimental replication and 75% for biological replication.',
              default=0.95,
              type=click.FloatRange(0, 1))
@click.option('--replication_rolling_size',
              help='Replication rolling size to calcultate replication rate',
              default=1000,
              type=int)
@click.option('--replication_num_sample',
              help='Number of samples which region need to be observed for replication',
              default=2,
              type=int)
@click.option('--replication_min_count',
              help='Minimum count needed to recognize region as expressed',
              default=1,
              type=int)
@click.option('--non_replicates_read_threhold',
              help='Minimum read count need for the samples without replication.'
              ' If there is not replicate samples for the sample, '
              ' this default cutoff will be applied.',
              default=10,
              type=int)
def cli_lapa_tss(alignment, fasta, annotation, chrom_sizes, output_dir,
                 mapq=10, cluster_extent_cutoff=3, cluster_ratio_cutoff=0.05,
                 cluster_window=25, min_replication_rate=0.95,
                 replication_rolling_size=1000, replication_num_sample=2,
                 replication_min_count=1, non_replicates_read_threhold=10):
    '''
    CLI interface for lapa tss cluster calling.
    '''
    lapa_tss(alignment, fasta, annotation, chrom_sizes, output_dir,
             cluster_extent_cutoff=cluster_extent_cutoff,
             cluster_ratio_cutoff=cluster_ratio_cutoff,
             cluster_window=cluster_window, mapq=mapq,
             min_replication_rate=min_replication_rate,
             replication_rolling_size=replication_rolling_size,
             replication_num_sample=replication_num_sample,
             replication_min_count=replication_min_count,
             non_replicates_read_threhold=non_replicates_read_threhold)


@click.command()
@click.option('--alignment',
              help='Path of the bam file. Start and end position of '
              'each read in the file will be overlaped against the '
              'tss/poly(A) cluster and annotated accordingly.',
              required=True)
@click.option('--lapa_dir',
              help='LAPA output directory of generated before with `lapa` command.',
              required=True)
@click.option('--lapa_tss_dir',
              help='LAPA output directory of generated before with `lapa_tss` command.',
              required=True)
@click.option('--output',
              help='Output path to .csv file which contains linking reads.',
              required=True)
@click.option('--mapq',
              help='Minimum read quality to required for linking.',
              default=10,
              type=int)
@click.option('--min_read_length',
              help='Minimum read quality to required for linking.',
              default=10,
              type=int)
@click.option('--dataset',
              help='Which dataset to use in linking. '
              'Valid options (`all`, `raw_all`, or dataset)',
              default='all')
def cli_lapa_link_tss_to_tes(alignment, lapa_dir, lapa_tss_dir, output,
                             mapq=10, min_read_length=100, dataset='all'):
    '''
    CLI interface for detecting of linking reads. Linking reads are
    the reads start in a tss cluster and in poly(A) cluster. Linking
    reads represents transcript complete isoforms.
    '''
    df = link_tss_to_tes(alignment, lapa_dir, lapa_tss_dir, mapq=mapq,
                         min_read_length=min_read_length, dataset=dataset)
    df.to_csv(output, index=False)


@click.command()
@click.option('--links',
              help='Path to linking read file generated '
              'with `lapa_link_tss_to_tes` command',
              required=True)
@click.option('--read_annot',
              help='read_annot of TALON annotating read,'
              ' transcript assignments.',
              required=True)
@click.option('--gtf_input',
              help='Input gtf file to extract splice chains.',
              required=True)
@click.option('--gtf_output',
              help='Output corrected gtf contains trascripts '
              'with tss/poly(A) end support.',
              required=True)
@click.option('--abundance_input',
              help='Input abundance file of TALON which contains '
              'abundance of each transcript.',
              required=True)
@click.option('--abundance_output',
              help='Update abundance file which calculated '
              'based on abundance of linking reads.',
              required=True)
@click.option('--keep_unsupported',
              help='Keep transcripts without tss and tes support in '
              'the original gtf. If true transcript created with '
              'non-linking reads (partial) in the original files '
              'are kept gtf and abundance.',
              is_flag=True)
def cli_lapa_correct_talon(links, read_annot, gtf_input, gtf_output,
                           abundance_input, abundance_output,
                           keep_unsupported=False):
    '''
    CLI interface for create GTF file with tss/poly(A) cluster support
    based on the linking reads and using splice chain of TALON.
    '''
    correct_talon(links, read_annot, gtf_input, gtf_output,
                  abundance_input, abundance_output,
                  keep_unsupported=keep_unsupported)
