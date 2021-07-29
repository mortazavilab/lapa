import click
from longread_postprocessing.tes import longread_postprocessing_tes


def longread_postprocessing_tss():
    pass


def longread_postprocessing_splicing():
    pass


@click.group()
def cli():
    pass


@cli.command()
@click.option('--alignment', help='Bam or Sam files')
@click.option('--annotation', help='Standart transcriptome Annotation (Encode or Ensembl gtf)')
@click.option('--fasta', help='Genome reference (Encode or Ensembl fasta)')
@click.option('--chrom_sizes', help='Chrom sizes files (can be generated with `faidx fasta -i chromsizes > chrom_sizes`)')
@click.option('--output_dir', help='Output directory')
def tes(alignment, fasta, annotation, chrom_sizes, output_dir):
    longread_postprocessing_tes(alignment, fasta, annotation,
                                chrom_sizes, output_dir)


@cli.command()
def tss():
    pass


@cli.command()
def splicing():
    pass


@cli.command()
def gtf_fix_tes():
    pass


if __name__ == '__main__':
    cli()
