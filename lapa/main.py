import click
from lapa.lapa import lapa


@click.command()
@click.option('--alignment', help='Bam or Sam files')
@click.option('--annotation', help='Standart transcriptome Annotation (Encode or Ensembl gtf)')
@click.option('--fasta', help='Genome reference (Encode or Ensembl fasta)')
@click.option('--chrom_sizes', help='Chrom sizes files (can be generated with `faidx fasta -i chromsizes > chrom_sizes`)')
@click.option('--output_dir', help='Output directory')
def cli(alignment, fasta, annotation, chrom_sizes, output_dir):
    lapa(alignment, fasta, annotation, chrom_sizes, output_dir)


if __name__ == '__main__':
    cli()
