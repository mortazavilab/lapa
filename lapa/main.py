import click
from lapa.lapa import lapa, lapa_tss
from lapa.link import link_tss_to_tes
from lapa.correction import correct_talon


@click.command()
@click.option('--alignment',
              help='Bam or csv file.',
              required=True)
@click.option('--fasta',
              help='Genome reference (Encode or Ensembl fasta)',
              required=True)
@click.option('--annotation',
              help='Standart transcriptome Annotation (Encode or Ensembl gtf)',
              required=True)
@click.option('--chrom_sizes',
              help='Chrom sizes files (can be generated with)'
              '`faidx fasta -i chromsizes > chrom_sizes`)',
              required=True)
@click.option('--output_dir',
              help='Output directory',
              required=True)
@click.option('--counting_method',
              help='Counting method either `end` or `tail`',
              type=click.Choice(['end', 'tail'], case_sensitive=True),
              default='end')
@click.option('--min_tail_len',
              help='Minimum tail length for tail counting strategy.'
              'This parameter will be ignored for end counting',
              default=10,
              type=int)
@click.option('--min_percent_a',
              help='Minimum percentage of A bp while seeking for tails. '
              'This parameter will be ignored for end counting',
              default=0.9,
              type=float)
@click.option('--mapq',
              help='Minimum read quality to required for tes calling',
              default=10,
              type=int)
@click.option('--cluster_extent_cutoff',
              help='Number of reads to initialized and terminated cluster',
              default=3,
              type=int)
@click.option('--cluster_ratio_cutoff',
              help='Percentage of coverage change for initialize cluster',
              default=0.05,
              type=click.FloatRange(0, 1))
@click.option('--cluster_window',
              help='Patience threshold to wait for termination cluster'
              'if number of reads subceed `the cluster_extent_cutoff`',
              default=25,
              type=int)
@click.option('--min_replication_rate',
              help='Minimum replication rate to include cluster in replicated clusters',
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
def cli_lapa(alignment, fasta, annotation, chrom_sizes, output_dir,
             counting_method, min_tail_len=10, min_percent_a=0.9, mapq=10,
             cluster_extent_cutoff=3, cluster_ratio_cutoff=0.05, cluster_window=25,
             min_replication_rate=0.95, replication_rolling_size=1000,
             replication_num_sample=2, replication_min_count=1):
    '''
    Click command line interface for lapa polyA cluster calling.
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
         replication_min_count=replication_min_count)


@click.command()
@click.option('--alignment',
              help='Bam files',
              required=True)
@click.option('--fasta',
              help='Genome reference (Encode or Ensembl fasta)',
              required=True)
@click.option('--annotation',
              help='Standart transcriptome Annotation (Encode or Ensembl gtf)',
              required=True)
@click.option('--chrom_sizes',
              help='Chrom sizes files (can be generated with)'
              '`faidx fasta -i chromsizes > chrom_sizes`)',
              required=True)
@click.option('--output_dir',
              help='Output directory',
              required=True)
@click.option('--mapq',
              help='Minimum read quality to required for tss calling',
              default=10,
              type=int)
@click.option('--cluster_extent_cutoff',
              help='Number of reads to initialized and terminated cluster',
              default=3,
              type=int)
@click.option('--cluster_ratio_cutoff',
              help='Percentage of coverage change for initialize cluster',
              default=0.05,
              type=click.FloatRange(0, 1))
@click.option('--cluster_window',
              help='Patience threshold to wait for termination cluster'
              'if number of reads subceed `the cluster_extent_cutoff`',
              default=25,
              type=int)
@click.option('--min_replication_rate',
              help='Minimum replication rate to include cluster in replicated clusters',
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
def cli_lapa_tss(alignment, fasta, annotation, chrom_sizes, output_dir,
                 mapq=10, cluster_extent_cutoff=3, cluster_ratio_cutoff=0.05,
                 cluster_window=25, min_replication_rate=0.95,
                 replication_rolling_size=1000, replication_num_sample=2,
                 replication_min_count=1):
    '''
    Click command line interface for lapa tss cluster calling.
    '''
    lapa_tss(alignment, fasta, annotation, chrom_sizes, output_dir,
             cluster_extent_cutoff=cluster_extent_cutoff,
             cluster_ratio_cutoff=cluster_ratio_cutoff,
             cluster_window=cluster_window, mapq=mapq,
             min_replication_rate=min_replication_rate,
             replication_rolling_size=replication_rolling_size,
             replication_num_sample=replication_num_sample,
             replication_min_count=replication_min_count)


@click.command()
@click.option('--alignment',
              help='Bam files',
              required=True)
@click.option('--lapa_dir',
              help='LAPA output directory of generated before with `lapa` command',
              required=True)
@click.option('--lapa_tss_dir',
              help='LAPA output directory of generated before with `lapa_tss` command',
              required=True)
@click.option('--output',
              help='Output path',
              required=True)
@click.option('--mapq',
              help='Minimum read quality to required for linking',
              default=10,
              type=int)
@click.option('--min_read_length',
              help='Minimum read quality to required for linking',
              default=10,
              type=int)
@click.option('--dataset',
              help='Which dataset to use in linking. '
              'Valid options (`all`, `raw_all`, or dataset)',
              default='all')
def cli_lapa_link_tss_to_tes(alignment, lapa_dir, lapa_tss_dir, output,
                             mapq=10, min_read_length=100, dataset='all'):
    df = link_tss_to_tes(alignment, lapa_dir, lapa_tss_dir, mapq=mapq,
                         min_read_length=min_read_length, dataset=dataset)
    df.to_csv(output, index=False)


@click.command()
@click.option('--links',
              help='links',
              required=True)
@click.option('--read_annot',
              help='read_annot file output by TALON',
              required=True)
@click.option('--gtf_input',
              help='Input gtf file need to corrected',
              required=True)
@click.option('--gtf_output',
              help='Output corrected gtf file',
              required=True)
@click.option('--abundance_input',
              help='Input abundance file for correction',
              required=True)
@click.option('--abundance_output',
              help='After abundance file after correction',
              required=True)
@click.option('--keep_unsupported',
              help='Keep transcripts without tss and tes support in the original gtf',
              is_flag=True)
def cli_lapa_correct_talon(links, read_annot, gtf_input, gtf_output,
                           abundance_input, abundance_output,
                           keep_unsupported=False):
    correct_talon(links, read_annot, gtf_input, gtf_output,
                  abundance_input, abundance_output,
                  keep_unsupported=keep_unsupported)
