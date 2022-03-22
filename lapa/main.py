import click
from lapa.lapa import lapa
from lapa.read import correct_gtf


@click.command()
@click.option('--alignment',
              help='Bam or Sam files',
              required=True)
@click.option('--fasta',
              help='Genome reference (Encode or Ensembl fasta)',
              required=True)
@click.option('--annotation',
              help='Standart transcriptome Annotation (Encode or Ensembl gtf)',
              required=True)
@click.option('--chrom_sizes',
              help='Chrom sizes files (can be generated with'
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
              help='Minimum percentage of A bp while seeking for tails. ',
              default=10,
              type=int)
@click.option('--cluster_extent_cutoff',
              help='Number of reads to initialized and terminated cluster',
              default=3,
              type=int)
@click.option('--cluster_window',
              help='Patience threshold to wait for termination cluster'
              'if number of reads subceed `the cluster_extent_cutoff`',
              default=25,
              type=int)
def cli(alignment, fasta, annotation, chrom_sizes, output_dir, counting_method,
        min_tail_len=10, min_percent_a=0.9, mapq=10,
        cluster_extent_cutoff=3, cluster_window=25):
    lapa(alignment, fasta, annotation, chrom_sizes, output_dir,
         counting_method, min_tail_len=min_tail_len,
         min_percent_a=min_percent_a, mapq=mapq)


@click.command()
@click.option('--gtf_input',
              help='Input gtf file need to corrected',
              required=True)
@click.option('--gtf_output',
              help='Output corrected gtf file',
              required=True)
@click.option('--lapa_dir',
              help='LAPA output directory of generated before',
              required=True)
@click.option('--read_annot',
              help='read_annot file output by TALON',
              required=True)
@click.option('--fasta',
              help='Genome reference (Encode or Ensembl fasta)',
              required=True)
@click.option('--correct_tss', is_flag=True,
              help="If TSS corrected as well",
              default=True, show_default=True)
def lapa_correct_talon_gtf(gtf_input, gtf_output, lapa_dir,
                           read_annot, fasta, correct_tss):
    correct_gtf(gtf_input, gtf_output, lapa_dir,
                read_annot, fasta, correct_tss)
