from collections import Counter
import pandas as pd
from lapa import lapa
from lapa.lapa import Lapa
from lapa.utils.io import cluster_col_order, tss_cluster_col_order, \
    read_talon_read_annot, read_tss_cluster, read_polyA_cluster
from conftest import fasta, gtf, chrom_sizes, \
    quantseq_gm12_bam, sample_csv, quantseq_both_gm12_bam, read_annot


def test_lapa_annotate_cluster(tmp_path):

    _lapa = Lapa(fasta, gtf, chrom_sizes, tmp_path / 'lapa')

    df_cluster = pd.DataFrame({
        'Chromosome': ['chr17', 'chrM', 'chr17', 'chr17'],
        'Start': [100099, 1100, 43144800, 43145100],
        'End': [100101, 1101, 43144850, 43145130],
        'polyA_site': [100100, 1101, 43144825, 43145114],
        'count': [15, 11, 5, 10],
        'Strand': ['+', '+', '+', '+'],
        'fracA': [6, 6, 1, 1],
        'signal': ['100157@GATAAA', '100157@GATAAA', 'None@None', 'None@None']
    })

    df_annotate = _lapa.annotate_cluster(df_cluster)

    df_expected = pd.DataFrame({
        'Chromosome': ['chr17', 'chr17', 'chr17', 'chr17', 'chrM'],
        'Start': [100099, 43144800, 43145100, 43145100, 1100],
        'End': [100101, 43144850, 43145130, 43145130, 1101],
        'polyA_site': [100100, 43144825, 43145114, 43145114, 1101],
        'count': [15, 5, 10, 10, 11],
        'Strand': ['+', '+', '+', '+', '+'],
        'fracA': [6, 1, 1, 1, 6],
        'signal': ['100157@GATAAA', 'None@None', 'None@None', 'None@None', '100157@GATAAA'],
        'Feature': ['intergenic', 'exon', 'exon', 'exon', 'intergenic'],
        'gene_id': ['intergenic_0', 'ENSG00000198496.12', 'ENSG00000198496.12', 'ENSG00000267681.1', 'intergenic_1'],
        'gene_name': ['intergenic_0', 'NBR2', 'NBR2', 'CTD-3199J23.6', 'intergenic_1'],
        'annotated_site': [-1, -1, -1, -1, -1],
        'tpm': [
            15 * 1000000 / (15 + 11 + 10 + 5),
            5 * 1000000 / (15 + 11 + 10 + 5),
            10 * 1000000 / (15 + 11 + 10 + 5),
            10 * 1000000 / (15 + 11 + 10 + 5),
            11 * 1000000 / (15 + 11 + 10 + 5)
        ]
    })

    df_annotate = df_annotate.reset_index()
    del df_annotate['index']

    pd.testing.assert_frame_equal(df_annotate, df_expected,
                                  check_dtype=False, check_categorical=False)


def test_lapa_sample_cluster(tmp_path):

    df_sample_count = pd.DataFrame({
        'Chromosome': ['chr17', 'chr17'],
        'Start': [43144811, 43145111],
        'End': [43144812, 43145112],
        'Strand': ['+', '+'],
        'count': [10, 6]
    })

    df_cluster = pd.DataFrame({
        'Chromosome': ['chr17', 'chr17', 'chr17', 'chrM'],
        'Start': [100099, 43144800, 43145100, 1100],
        'End': [100101, 43144850, 43145130, 1101],
        'Strand': ['+', '+', '+', '+'],
        'polyA_site': [100100, 43144825, 43145114, 1101],
        'count': [15, 6, 10, 11],
        'fracA': [6, 1, 10, 1],
        'signal': ['100157@GATAAA', 'None@None', 'None@None', '100157@GATAAA'],
        'Feature': ['intergenic', 'exon', 'exon', 'intergenic'],
        'gene_id': ['', 'ENSG00000198496.12', 'ENSG00000198496.12', ''],
        'gene_name': ['', 'NBR2', 'NBR2', ''],
        'annotated_site': [-1, -1, 43145114, -1],
        'annotated': [False, False, True, False],
        'tpm': [
            15 * 1000000 / (15 + 11 + 10 + 6),
            10 * 1000000 / (15 + 11 + 10 + 6),
            11 * 1000000 / (15 + 11 + 10 + 6),
            6 * 1000000 / (15 + 11 + 10 + 6)
        ]
    })

    _lapa = Lapa(fasta, gtf, chrom_sizes, tmp_path / 'lapa_internal_priming',
                 filter_internal_priming=False)
    df_apa = _lapa.sample_cluster(df_cluster, df_sample_count)
    df_apa = _lapa.calculate_usage(df_apa)

    assert df_apa.shape[0] == 2
    assert all(df_apa['count'] == [10, 6])
    assert all(df_apa['usage'] == [10 / 16, 6 / 16])

    _lapa = Lapa(fasta, gtf, chrom_sizes, tmp_path / 'lapa')

    df_apa = _lapa.sample_cluster(df_cluster, df_sample_count)
    df_apa = _lapa.calculate_usage(df_apa)

    assert df_apa.shape[0] == 1

    row = df_apa.iloc[0]
    assert row['Chromosome'] == 'chr17'
    assert row['Start'] == 43144800
    assert row['End'] == 43144850
    assert row['Strand'] == '+'


def test_lapa_bam_pb(tmp_path):
    output_dir = tmp_path / 'lapa'

    lapa(quantseq_gm12_bam, fasta, gtf, chrom_sizes,
         output_dir, method='tail')

    # ├── raw_polyA_clusters.bed
    # ├── polyA_clusters.bed
    # ├── counts
    # │   ├── all_polyA_counts_neg.bw
    # │   └── all_polyA_counts_pos.bw
    # ├── coverage
    # │   ├── all_polyA_coverage_neg.bw
    # │   └── all_polyA_coverage_pos.bw
    # ├── dataset
    # │   └── quantseq3_gm12878_chr17_rep1.bed
    # ├── ratio
    # │   ├── all_polyA_ratio_neg.bw
    # │   └── all_polyA_ratio_pos.bw
    # ├── raw_sample
    # │   └── quantseq3_gm12878_chr17_rep1.bed
    # ├── logs
    # │   ├── progress.log
    # │   ├── warnings.log
    # │   └── final_stats.log
    # └── sample
    #     └── quantseq3_gm12878_chr17_rep1.bed

    assert (output_dir / 'counts').exists()
    assert (output_dir / 'counts' / 'all_polyA_counts_pos.bw').exists()
    assert (output_dir / 'counts' / 'all_polyA_counts_neg.bw').exists()

    assert (output_dir / 'coverage').exists()
    assert (output_dir / 'coverage' / 'all_polyA_coverage_pos.bw').exists()
    assert (output_dir / 'coverage' / 'all_polyA_coverage_neg.bw').exists()

    assert (output_dir / 'ratio').exists()
    assert (output_dir / 'ratio' / 'all_polyA_ratio_pos.bw').exists()
    assert (output_dir / 'ratio' / 'all_polyA_ratio_neg.bw').exists()

    assert (output_dir / 'dataset').exists()

    assert (output_dir / 'raw_sample').exists()
    assert (output_dir / 'raw_sample' /
            'quantseq3_gm12878_chr17_rep1.bed').exists()

    assert (output_dir / 'sample').exists()
    assert (output_dir / 'sample' /
            'quantseq3_gm12878_chr17_rep1.bed').exists()

    assert (output_dir / 'raw_polyA_clusters.bed').exists()

    assert (output_dir / 'logs' / 'progress.log').exists()
    assert (output_dir / 'logs' / 'warnings.log').exists()
    assert (output_dir / 'logs' / 'final_stats.log').exists()

    df_cluster = read_polyA_cluster(str(output_dir / 'raw_polyA_clusters.bed'))

    assert all(df_cluster.columns == cluster_col_order)

    assert set(df_cluster['Chromosome']) == {'chr17'}
    assert (df_cluster['End'] - df_cluster['Start']).max() < 100

    assert all(df_cluster['fracA'] <= 10)
    assert all(df_cluster['fracA'] >= 0)

    df_apa = read_polyA_cluster(
        str(output_dir / 'sample' / 'quantseq3_gm12878_chr17_rep1.bed'))

    assert all(df_apa.columns == cluster_col_order)

    counts = pd.Series(Counter(df_apa['Feature']))
    assert counts.idxmax() == 'three_prime_utr'

    df_cluster = read_polyA_cluster(str(output_dir / 'polyA_clusters.bed'))
    assert df_cluster.shape == df_apa.shape


def test_lapa_bam_quantseq(tmp_path):

    output_dir = tmp_path / 'lapa'

    lapa(quantseq_both_gm12_bam, fasta, gtf,
         chrom_sizes, output_dir, method='tail',
         cluster_ratio_cutoff=0.01, replication_rolling_size=10,
         min_replication_rate=0.95)

    df_cluster = read_polyA_cluster(str(output_dir / 'raw_polyA_clusters.bed'))

    assert all(df_cluster.columns == cluster_col_order)

    assert set(df_cluster['Chromosome']) == {'chr17'}
    assert (df_cluster['End'] - df_cluster['Start']).max() < 150

    assert all(df_cluster['fracA'] <= 10)
    assert all(df_cluster['fracA'] >= 0)

    df_raw = read_polyA_cluster(
        str(output_dir / 'raw_sample' / 'quantseq3_gm12878_chr17_rep1.bed'))

    assert all(df_raw.columns == cluster_col_order)

    counts = pd.Series(Counter(df_raw['Feature']))
    assert counts.idxmax() == 'three_prime_utr'

    df_replicated = read_polyA_cluster(
        str(output_dir / 'sample' / 'quantseq3_gm12878_chr17_rep1.bed'))

    assert df_raw.shape[0] > df_replicated.shape[0]


def test_lapa_read_csv(tmp_path):
    output_dir = tmp_path / 'lapa'

    lapa(sample_csv, fasta, gtf, chrom_sizes, output_dir,
         method='tail',
         cluster_ratio_cutoff=0.01,
         cluster_extent_cutoff=10)

    df_cluster = read_polyA_cluster(str(output_dir / 'raw_polyA_clusters.bed'))

    assert set(df_cluster['Chromosome']) == {'chr17'}

    assert (df_cluster['End'] - df_cluster['Start']).max() < 150

    assert all(df_cluster['fracA'] <= 10)
    assert all(df_cluster['fracA'] >= 0)

    df_apa = read_polyA_cluster(str(output_dir / 'sample' / 'gm12878_1.bed'))
    counts = pd.Series(Counter(df_apa['Feature']))

    assert all(df_apa.columns == cluster_col_order)

    counts = pd.Series(Counter(df_apa['Feature']))
    assert counts.idxmax() == 'three_prime_utr'


def test_lapa_tss(lapa_tss_read_annot):
    df_cluster = read_tss_cluster(lapa_tss_read_annot / 'raw_tss_clusters.bed')

    assert all(df_cluster.columns == tss_cluster_col_order)

    assert all(df_cluster['tss_site'] <= df_cluster['End'])
    assert all(df_cluster['tss_site'] >= df_cluster['Start'])

    assert (df_cluster['End'] - df_cluster['Start']).max() < 300
    assert set(df_cluster['Chromosome']) == {'chr17', 'ERCC-00060'}

    counts = df_cluster.drop_duplicates([
        'Chromosome', 'Start', 'End', 'Strand'])['count'].sum()

    df_tss = read_tss_cluster(lapa_tss_read_annot /
                              'sample' / 'gm12878.bed')
    assert all(df_tss.columns == tss_cluster_col_order)

    df_tss = read_tss_cluster(lapa_tss_read_annot /
                              'sample' / 'hepg2.bed')
    assert all(df_tss.columns == tss_cluster_col_order)


def test_lapa(lapa_read_annot):
    df_cluster = read_polyA_cluster(
        str(lapa_read_annot / 'raw_polyA_clusters.bed'))

    assert all(df_cluster.columns == cluster_col_order)

    assert all(df_cluster['polyA_site'] <= df_cluster['End'])
    assert all(df_cluster['polyA_site'] >= df_cluster['Start'])

    assert (df_cluster['End'] - df_cluster['Start']).max() < 150

    assert set(df_cluster['Chromosome']) == {'chr17', 'ERCC-00060'}

    counts = df_cluster.drop_duplicates([
        'Chromosome', 'Start', 'End', 'Strand'])['count'].sum()

    df_read_annot = read_talon_read_annot(read_annot)
    # at least 75% of the reads used in clustering
    assert df_read_annot.shape[0] * 0.75 < counts

    df_apa = read_polyA_cluster(
        str(lapa_read_annot / 'sample' / 'gm12878.bed'))
    assert all(df_apa == cluster_col_order)

    counts = pd.Series(Counter(df_apa['Feature']))
    assert counts.idxmax() == 'three_prime_utr'

    counts = df_apa.drop_duplicates([
        'Chromosome', 'Start', 'End', 'Strand'
    ])['count'].sum()

    _df = df_read_annot[df_read_annot['sample'] == 'gm12878']
    assert _df.shape[0] > counts
