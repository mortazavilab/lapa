import click
from lapa.lapa import lapa


@click.command()
@click.option('--alignment', help='Bam or Sam files')
@click.option('--fasta', help='Genome reference (Encode or Ensembl fasta)')
@click.option('--annotation', help='Standart transcriptome Annotation (Encode or Ensembl gtf)')
@click.option('--chrom_sizes', help='Chrom sizes files (can be generated with `faidx fasta -i chromsizes > chrom_sizes`)')
@click.option('--output_dir', help='Output directory')
@click.option('--method', help='Counting method either `end` or `tail`',
              type=click.Choice(['end', 'tail'], case_sensitive=True))
@click.option('--min_tail_len', default=10, type=int)
@click.option('--min_percent_a', default=0.9, type=float)
@click.option('--mapq', default=10, type=int)
def cli(alignment, fasta, annotation, chrom_sizes, output_dir, method,
        min_tail_len=10, min_percent_a=0.9, mapq=10):
    lapa(alignment, fasta, annotation, chrom_sizes, output_dir, method,
         min_tail_len=min_tail_len, min_percent_a=min_percent_a, mapq=mapq)


if __name__ == '__main__':
    cli()
