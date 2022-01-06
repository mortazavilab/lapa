import click
from lapa.lapa import lapa
from lapa.read import correct_gtf


@click.command()
@click.option('--alignment', help='Bam or Sam files')
@click.option('--fasta', help='Genome reference (Encode or Ensembl fasta)')
@click.option('--annotation', help='Standart transcriptome Annotation (Encode or Ensembl gtf)')
@click.option('--chrom_sizes', help='Chrom sizes files (can be generated with `faidx fasta -i chromsizes > chrom_sizes`)')
@click.option('--output_dir', help='Output directory')
@click.option('--counting_method', help='Counting method either `end` or `tail`',
              type=click.Choice(['end', 'tail'], case_sensitive=True), default='end')
@click.option('--min_tail_len', default=10, type=int)
@click.option('--min_percent_a', default=0.9, type=float)
@click.option('--mapq', default=10, type=int)
@click.option('--cluster_extent_cutoff', default=3, type=int)
@click.option('--cluster_window', default=25, type=int)
def cli(alignment, fasta, annotation, chrom_sizes, output_dir, counting_method,
        min_tail_len=10, min_percent_a=0.9, mapq=10,
        cluster_extent_cutoff=3, cluster_window=25):
    lapa(alignment, fasta, annotation, chrom_sizes, output_dir, counting_method,
         min_tail_len=min_tail_len, min_percent_a=min_percent_a, mapq=mapq)


@click.command()
@click.option('--gtf_input', help='Input gtf file need to corrected')
@click.option('--gtf_output', help='Output corrected gtf file')
@click.option('--lapa_dir', help='LAPA output directory of generated before')
@click.option('--read_annot', help='read_annot file output by TALON')
@click.option('--fasta', help='Genome reference (Encode or Ensembl fasta)')
@click.option('--correct_tss', is_flag=True, help="If TSS corrected as well",
              default=True, show_default=True)
def lapa_correct_talon_gtf(gtf_input, gtf_output, lapa_dir,
                           read_annot, fasta, correct_tss):
    correct_gtf(gtf_input, gtf_output, lapa_dir,
                read_annot, fasta, correct_tss)

# if __name__ == '__main__':
#     cli()
