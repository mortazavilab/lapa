import pyranges as pr
from longread_postprocessing import longread_postprocessing_tes, \
    read_polyA_cluster
from conftest import fasta, gtf, read_annot, chrom_sizes


def test_longread_postprocessing_tes(tmp_path):
    output_dir = tmp_path / 'tes'

    longread_postprocessing_tes(
        read_annot, fasta, gtf, chrom_sizes, output_dir)

    df_cluster = read_polyA_cluster(str(output_dir / 'polyA_clusters.bed'))

    # TODO: write better test cases!
    df_cluster.columns == [
        'Chromosome', 'Start', 'End', 'Strand', 'polyA_site', 'count', 'fracA',
        'singal', 'Feature', 'canonical_site', 'canonical', 'tpm']

    assert set(df_cluster['Chromosome']) == {'chr17', 'ERCC-00060'}
