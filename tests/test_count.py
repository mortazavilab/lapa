import pytest
import pysam
import pyBigWig
import pandas as pd
import numpy as np
from lapa.utils.io import read_talon_read_annot
from lapa.count import PolyaTailCounter, FivePrimeCounter, \
    ThreePrimeCounter,  FivePrimeCounter, \
    TesMultiCounter, TssMultiCounter
from conftest import quantseq_gm12_bam, short_bam, chr17_chrom_sizes, \
    read_annot, pb_brca1_bam


@pytest.fixture
def three_prime_counter_brca1():
    return ThreePrimeCounter(pb_brca1_bam)


@pytest.fixture
def five_prime_counter_brca1():
    return FivePrimeCounter(pb_brca1_bam)


def test_ThreePrimeCounter_to_df(three_prime_counter_brca1, tmp_path):

    df = three_prime_counter_brca1 \
        .to_df() \
        .set_index(['Chromosome', 'End', 'Strand'])


    three_prime_counter_brca1.to_bigwig(chr17_chrom_sizes, tmp_path)
    
    assert df.loc[('chr17', 43044295, '-'), 'count'] == 25
    assert df.loc[('chr17', 43044295, '-'), 'coverage'] == 31

    assert df.loc[('chr17', 43044301, '-'), 'count'] == 20
    assert df.loc[('chr17', 43044301, '-'), 'coverage'] == 65


def test_FivePrimeCounter_to_df(five_prime_counter_brca1):

    df = five_prime_counter_brca1 \
        .to_df() \
        .set_index(['Chromosome', 'End', 'Strand'])

    assert df.loc[('chr17', 43125360, '-'), 'count'] == 6
    assert df.loc[('chr17', 43125360, '-'), 'coverage'] == 29

    assert df.loc[('chr17', 43125382, '-'), 'count'] == 3
    assert df.loc[('chr17', 43125382, '-'), 'coverage'] == 8


def load_reads(seq, cigar, flag, pos):
    bam = pysam.AlignmentFile(quantseq_gm12_bam)
    read = {
        'name': 'NS500169:839:HMNMLAFX2:1:11309:11029:9457',
        'flag': flag,
        'ref_name': 'chr17',
        'ref_pos': pos,
        'map_quality': '255',
        'cigar': cigar,
        'next_ref_name': '*',
        'next_ref_pos': '0',
        'length': '0',
        'seq': seq,
        'qual': '5' * len(seq),
        'tags': []
    }
    return pysam.AlignedSegment.from_dict(read, bam.header)


def test_ThreePrimeCounter_count_read():
    seq = 'GCGCAGCAGGGGTGGTCGCCATGGAGACGCGTGGCCCTGGCCTGGCGGTCCGCGCTGAGAGTCGCCGATTAGTCGGCATCGGGCCTCGGGCGCCCCCGGGGCGGGTTGGGTTGCAGCCCAGCGGGCGGCTGGACCGCCGCGGTGGGGCGGGGACAATGGGGTACAAGGACAACGACGGCGAGGAGGAGGAGCGGGAGGGCGGCGCCGCGGGCCCGCGGGGGTCTAGACTGCCCCCCATCACAGGCGGCGCCTCCGAGCTGGCCAAACGGAAGGTGAAGAAGAAAAAAAGGAAGAAGAAGACCAAGGGGTCTGGCAAGGGGGACGGGATCTTGCTCTGGGGAAAGCCAGCTGCAGTGTTGTGTGGGCTGCCTGAGGAAAGGTGCCCTCCGGAAAGAAATGATGTCTCCGCCCAACAGCATGGTGGACCTGAGTCCTGCTACCAGCCACATGAGTGTGCTTGGAAGCAGGTGGAGCCCAGTTGAGCCTTGAGGTGCCTCCCTCTGTGAGAGACCCCAAGCCGGAGACATCCGGCCAAGAGGCACCCAGATTCCTGCCCCACAGAAACTGATACAATAACGTTGTTTAAAGCCACTACATTTGGAGGTAATCTTGTTACACAGCAATAACTGACTAATACAGTCGGCCTCTTCCATCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
    cigar = '324M10025N335M100S'
    flag = '0'
    pos = '410323'
    read = load_reads(seq, cigar, flag, pos)

    polyA_site = ThreePrimeCounter.count_read(None, read)
    assert polyA_site == 421001


def test_PolyaTailCounter_detect_polyA():
    seq = 'CGCAATATGACGTTTCCATTTACTTTGGATTATATGTCATTATAAATATTAACAAATAAGACTTAAAAAGGACACCTTCGGGTAGGTCAGACCAAAATACAAAACTTGTCTGTGGGACTGCAGTTTGGA'
    cigar = '127M2S'
    flag = '16'
    pos = '143694'
    read = load_reads(seq, cigar, flag, pos)

    polyA_site, tail_len, percent_A = PolyaTailCounter.detect_polyA_tail(read)

    assert polyA_site is None
    assert tail_len == 0
    assert percent_A == 0

    seq = 'GCGCAGCAGGGGTGGTCGCCATGGAGACGCGTGGCCCTGGCCTGGCGGTCCGCGCTGAGAGTCGCCGATTAGTCGGCATCGGGCCTCGGGCGCCCCCGGGGCGGGTTGGGTTGCAGCCCAGCGGGCGGCTGGACCGCCGCGGTGGGGCGGGGACAATGGGGTACAAGGACAACGACGGCGAGGAGGAGGAGCGGGAGGGCGGCGCCGCGGGCCCGCGGGGGTCTAGACTGCCCCCCATCACAGGCGGCGCCTCCGAGCTGGCCAAACGGAAGGTGAAGAAGAAAAAAAGGAAGAAGAAGACCAAGGGGTCTGGCAAGGGGGACGGGATCTTGCTCTGGGGAAAGCCAGCTGCAGTGTTGTGTGGGCTGCCTGAGGAAAGGTGCCCTCCGGAAAGAAATGATGTCTCCGCCCAACAGCATGGTGGACCTGAGTCCTGCTACCAGCCACATGAGTGTGCTTGGAAGCAGGTGGAGCCCAGTTGAGCCTTGAGGTGCCTCCCTCTGTGAGAGACCCCAAGCCGGAGACATCCGGCCAAGAGGCACCCAGATTCCTGCCCCACAGAAACTGATACAATAACGTTGTTTAAAGCCACTACATTTGGAGGTAATCTTGTTACACAGCAATAACTGACTAATACAGTCGGCCTCTTCCATCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
    cigar = '324M10025N335M100S'
    flag = '0'
    pos = '410323'
    read = load_reads(seq, cigar, flag, pos)

    polyA_site, tail_len, percent_A = PolyaTailCounter.detect_polyA_tail(read)

    assert polyA_site == 421001
    assert tail_len == 100
    assert percent_A == 1

    polyA_site, tail_len, percent_A = PolyaTailCounter.detect_polyA_tail(
        read, count_aligned=True)
    assert tail_len == 105

    seq = 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATGACGTTTCCATTTACTTTGGATTATATGTCATTATAAATATTAACAAATAAGACTTAAAAAGGACACCTTCGGGTAGGTCAGACCAAAATACAAAACTTGTCTGTGGGACTGCAGTTTGAGGACAGTGTCTGCAGCCGTCACATGGCAGCAAAACGGGGTTAAGCAGTGCACGAGAGTCTGCGTCGACGACAGCCAGAGTCCATGCATCGGGAGGTTCACTCGGTTTGCGAAGGAACAACGGGCTCGGCATGCACGGGCTCGGCGGCGGCGGACGGGCCGGGGCGCAGTTCCCCGCGCTCGCCACTAGAGGTCAGGAGGTGACCGCTTCGGGGCTGGAAGACGGGCCCGTCGGGGATTGGCGCAGGCGGCGGGCGGGGCGGCGGGCGGGGCGGCGCTGGAGGCAGCGCCTGGTTACTGACACCTGGAATGACTTTTTTTTTTTGGCATCAGATTTCCTGTCTTTGTGGGGATGATGGACCCGAGTAAAGATGCCCATTCGGGGTCAAAGGCAGAGCCGCTTCTGCAGCTTCTCAAAGCGTTGTTTGTTTGTTTTTTTTTCTGAGACGGAGTCTTGCTCTGTCGCCCAGGCTGGAGTGCAGTGCCGCGATCTTGGCTCACTGCAGCCTTCACCTCCCGGGTTCAAGCGATTCTCGTGCCTCAGCCTCCCTGAGTAGCTGGGACTACAGGCGTGCGCCACCACACCCGGCTAATTTTTGTATTTTTAGTAGAGATGGGGTTTCGCCCTGTTGGCCAGGCGGGTTTCGAACTCCTGACCTCAGGTGATCCGCCAGCCTCGGCCTCCCAAAGTGCTGCAATTACAGGCGTGAGCCATCGTGCCTCAAAACGCTTTAACAGAAAGACAATCTGCACGGGATCTAAAAGGGTGCTGAGATCCTAGGGAAGGAAGGATCCAAACTTCCTGGGGAGTTCCTGCCCGAGTGCCTGTGCTGCCCCTGGGCTGGCTGGCCAGTAAGCCCGCCTCCCAGCCTGACTGTCCCCATCTTTCGGTCCCAGCCCCATTTGCACAGCCTGGGCAGTAGAGGGCCCTGGACTGGGGGGCTGGAGTTCTGGGTTCTTGTCCCAGCTGTGCCCCTTCAGGCCCCTTTCACCCACTGGGCCTCACATTCCCCATGTGCCTAATAAGGAAGTCATGGTAGGTGGGGGGTACAGCCTCTTCTCGCTTTGCCATTCACTTCTGTGTCTTCAGAGAGGCTGTCAAATCTCCTCTCTGAATGAGACCCCCAAAAAAGCAAAGCTAAGAAGATACCCAGAGCTTAACTAGACACCAGGCCTTTTAGAAATAGACCACCTCTTACCTTAGGCCCCCAGAGGGTGCCCATTCTGTGTGGAGAAAAGAGGAGCCCTTGCCTCAGCCCCCAGAGGCTAGGGTGGGGTGGCTGAGTTTTGGGGCCAGGTTAGACACCTCTGGGGAAGCCTGAAGTAGCACGGATGGTTTCAAAGCCAGCGGATGAGGTGGCAGGACAGTGACAAGACCCCAGGTCTCCATCATGCACCTGAGCTTACTGAGCCTTCACCTGGTGTCTCTGCTGAGTCTCCGACAGACCAAGAGGGAAGGGACTGGGGTACCGACCCCCAGAGAGAGAAAGCGGCAGTCAGAGGACCTGGGTTTGAGTCCTGGTTCACCCCTTCCTGGCTGTGTGGCCTTGGCAAACTACTCAGACTCTGGGAACCTGTTTCACCTGCAAGATGGGGATGAGAATCAAACCCACCTTGCAGGGCTGTGAACATGCGTTCAGACCAGTGCATGCAAAAACCTTCACACAAAAACCCTTCCTAAAGGCAGGGAGACAAGTCTCCAGGAGCAGCAGGTGGCCAGGGACTGTGTGGGGGCTGGGGTCCTGTTTTCCCCGCAACCTGGGAAAGGCCTGATGGGCACTTGGTAGGATTCAAATCATAACCCTGGACCTCAGGTGGTGGTGTGTGCTTTGTGTCTGCAGTGGAAGGTCCCACCGTCGTGCCTGTGAGTCCCTCTATCGTGTTGAAGGGGTGGGCAGCAGGCGGGGGAGCCCAGGCTGGCAGCTAGCAGGACCTCTCTGGTGTGAGCTCAGCACGCCAATTCCCCTGAAGCGTGGTGTCCAGCACTCTGGGCTGGGGGCTGTGGATCCTGGTTCCAGCTGTGTGGGGCCTGGAAGGCCCTGGGCAGGTCACCTGACCTCTCTGGGCCTGTTTCCATCACACACTGATGGGCTGAGCACACTGGGACTCTGGCTGTGCAACTCCTTGGCTCTGCATTTGTTCACCCAGCGTTCCTGAGGGGCCCTTGGTAGGCAGAGAAAGTTTGTGGGCTTCAGGCATGGGCCCCAAGATTCAGATACTCTCAAGCCTCCTGGGGAGTCTCACTCAGGGGAGCAGACAGGCCCACCAGCCAGGGTGATTCTTGCTCAGCTCCTGGAGGGTGGAGCTGGCCCCACAGGCCTCCCCAAAGACAAGGCCCTGGACGGCCACTGATTACTCCCAAAAGGGACATCTGTGGCGGTAGTGGGGACCAAGATGCCACGAGCAAGCACCCGGAGGCCCCCAGTGTGAGCCATGGAGTGGAGGGGAGGGGAAAGGGCAGAGTCAGGACGTGTAGGAATGCTTGCTTTTTTTCCAAGCACAAGGGACCCTTTTCTCCACTGCAGCTGACCTGATGCTTATGCCAGGAAGGAGGAGGGGCGGGCTCCGTCCTTGAGGTCCCTCAGGAGTAGAAAGAGATCAGAGTGGGAGACTTGGGTCTGAGGTCTGAACTTGAGCCTACACCAGTTTCTCCATGGTGTTTCCATCACGTCCCCCACCTCCTGCCTCGAGCCTCACACCTTCCTAACAAACCCTCCCCTGGAGAGGAGACCCTGGGTCCACAGCACCCGGCCCCATTGGCTTTCTCTCTCCAGGCGTTTTAGACCACTCCCACTGCCTGGACTGCTTTTAGAAGCCCTCACTCACCCTCCAAGGCCTTACGGAAGCACCGCCTCCTCCAGGAAGCCGTCCCTGACCTCCCGGCAGAGTCAGCAGAGCCCACGGTACTTGCAGACTGGCCAAGCTGGTCTTGGTATTTCTCACCCCAGTTGGAATGGTTTGTTCCCAGCTTCCTGACTAGATGGGAACTCCCTGAGGGCAGGCCCTGTGTCTCATTCACCCCAGGGCTTGTGGAATCATTGAGTGAAGCATGTGCCAATTTACCCATCATCAGGAGCCTCCGGGAACTCTGGCAGACTTTCGGCCGGCAGGCCCGCTTCTTCCATCTGTCCAGCAGCGAAGGAGATGAGGGATGCAGTTAGGCTTTCTTGGGCTGGAGCAGCCAGTCTTCAAGGTCCCATCCTCCACC'
    cigar = '104S257M5D11M2I97M12D186M1I1030M1D357M2I57M1D31M1D474M1I108M1D253M1D272M1D190M'
    flag = '16'
    pos = '143699'
    read = load_reads(seq, cigar, flag, pos)

    polyA_site, tail_len, percent_A = PolyaTailCounter.detect_polyA_tail(
        read, count_aligned=False)

    assert polyA_site == 143699
    assert tail_len == 104
    assert percent_A == 1

    polyA_site, tail_len, percent_A = PolyaTailCounter.detect_polyA_tail(
        read, count_aligned=True)
    assert tail_len == 105

    seq = 'CAAATAGGAGCAGTATGATTTTTTCTTTTATTTTTAATTGACAAATAATAATTATATAAACAATGTAATTTCTAACAACAAAGGAATGCTAAAAAATATACAAAAAAAAAAAAAAAAAAAAGAGCGGAGGCGCCCCCGGGTGCACGCCGG'
    cigar = '1S119M30S'
    flag = '0'
    pos = '716075'
    read = load_reads(seq, cigar, flag, pos)

    polyA_site, tail_len, percent_A = PolyaTailCounter.detect_polyA_tail(read)

    assert polyA_site == 716174
    assert tail_len == 1
    assert percent_A == 1

    polyA_site, tail_len, percent_A = PolyaTailCounter.detect_polyA_tail(
        read, count_aligned=True)
    assert tail_len == 20

    seq = 'GGGCCGGGCGGGCCGCCGGGGGCGCTGCCGCCCTTTTTTTTTTTTTTTTTTAAGAGCAGAGTGAACTCTTTATTGATTATACAAATTACCACTATTTATTTTAAACCCAAAGTGACTTCAAAGTTGTTTTTTGGTTTTTAAAGGGGCTAC'
    flag = '16'
    cigar = '50S100M'
    pos = '515042'
    read = load_reads(seq, cigar, flag, pos)

    polyA_site, tail_len, percent_A = PolyaTailCounter.detect_polyA_tail(read)

    assert polyA_site == 515042
    assert tail_len == 17
    assert percent_A == 1

    polyA_site, tail_len, percent_A = PolyaTailCounter.detect_polyA_tail(
        read, count_aligned=True)
    assert tail_len == 18


@pytest.fixture
def tail_counter_short():
    return PolyaTailCounter(short_bam)


def test_PolyaTailCounter_iter_tailed_reads(tail_counter_short):
    tailed = sum(1 for i in tail_counter_short.iter_tailed_reads())
    assert tailed > 3600


def test_PolyaTailCounter_tail_len_dist(tail_counter_short):
    tail_len = tail_counter_short.tail_len_dist()
    assert tail_len.shape[0] > 20


def test_PolyaTailCounter_save_tailed_reads(tail_counter_short, tmp_path):
    tailed_bam = str(tmp_path / 'tailed.bam')
    tail_counter_short.save_tailed_reads(tailed_bam)

    tailed_bam = pysam.AlignmentFile(tailed_bam, 'rb')
    tailed = sum(1 for i in tailed_bam)
    assert tailed > 3600


def test_PolyaTailCounter_count(tail_counter_short):
    tes = tail_counter_short.count()
    assert len(tes) > 1000


def test_PolyaTailCounter_to_df(tail_counter_short):
    df = tail_counter_short.to_df()
    assert df.shape[0] > 1000


def test_TesMultiCounter_samples(tmp_path):
    # Test multiple samples
    df = pd.DataFrame({
        'sample': ['short_rep1', 'short_rep2', 'short', 'short'],
        'dataset': ['short', 'short', 'short', 'short'],
        'path': [short_bam, short_bam, short_bam, short_bam]
    })

    counter = TesMultiCounter(df, method='tail')
    df_all, tes = counter.to_df()

    pd.testing.assert_frame_equal(
        tes['short_rep1'],
        tes['short_rep2']
    )
    _df = tes['short_rep1'].copy()
    _df['count'] *= 4
    _df['coverage'] *= 4
    pd.testing.assert_frame_equal(df_all, _df)

    _df = tes['short_rep1'].copy()
    _df['count'] *= 2
    _df['coverage'] *= 2
    pd.testing.assert_frame_equal(tes['short'], _df)

    # Test bigwig exporting
    output_dir = tmp_path / 'multi_sample'
    output_dir.mkdir()

    TesMultiCounter._to_bigwig(df_all, tes, chr17_chrom_sizes,
                               output_dir)

    assert {i.name for i in (output_dir / 'counts').iterdir()} == {
        'all_polyA_counts_pos.bw',
        'all_polyA_counts_neg.bw',
        'short_polyA_counts_pos.bw',
        'short_polyA_counts_neg.bw',
        'short_rep1_polyA_counts_pos.bw',
        'short_rep1_polyA_counts_neg.bw',
        'short_rep2_polyA_counts_pos.bw',
        'short_rep2_polyA_counts_neg.bw'
    }

    assert {i.name for i in (output_dir / 'coverage').iterdir()} == {
        'all_polyA_coverage_pos.bw',
        'all_polyA_coverage_neg.bw',
        'short_polyA_coverage_pos.bw',
        'short_polyA_coverage_neg.bw',
        'short_rep1_polyA_coverage_pos.bw',
        'short_rep1_polyA_coverage_neg.bw',
        'short_rep2_polyA_coverage_pos.bw',
        'short_rep2_polyA_coverage_neg.bw'
    }

    assert {i.name for i in (output_dir / 'ratio').iterdir()} == {
        'all_polyA_ratio_pos.bw',
        'all_polyA_ratio_neg.bw',
        'short_polyA_ratio_pos.bw',
        'short_polyA_ratio_neg.bw',
        'short_rep1_polyA_ratio_pos.bw',
        'short_rep1_polyA_ratio_neg.bw',
        'short_rep2_polyA_ratio_pos.bw',
        'short_rep2_polyA_ratio_neg.bw'
    }

    bw_pos = pyBigWig.open(
        str(output_dir / 'counts' / 'all_polyA_counts_pos.bw'))
    bw_neg = pyBigWig.open(
        str(output_dir / 'counts' / 'all_polyA_counts_neg.bw'))
    count = np.nansum(
        bw_pos.values('chr17',
                      int(df_all['End'].min() - 1),
                      int(df_all['End'].max()))
    ) + np.nansum(
        bw_neg.values('chr17',
                      int(df_all['End'].min() - 1),
                      int(df_all['End'].max()))
    )
    assert df_all['count'].sum() == count

    bw_pos = pyBigWig.open(
        str(output_dir / 'counts' / 'short_polyA_counts_pos.bw'))
    bw_neg = pyBigWig.open(
        str(output_dir / 'counts' / 'short_polyA_counts_neg.bw'))
    count = np.nansum(
        bw_pos.values('chr17',
                      int(df_all['End'].min() - 1),
                      int(df_all['End'].max()))
    ) + np.nansum(
        bw_neg.values('chr17',
                      int(df_all['End'].min() - 1),
                      int(df_all['End'].max()))
    )
    assert df_all['count'].sum() / 2 == count


def test_TesMultiCounter_single_sample(tmp_path):
    # Test single sample
    df = pd.DataFrame({
        'sample': ['short', 'short'],
        'path': [short_bam, short_bam],
        'dataset': ['short', 'short']
    })

    counter = TesMultiCounter(df, method='tail')
    df_all, tes = counter.to_df()

    output_dir = tmp_path / 'single_sample'
    output_dir.mkdir()
    TesMultiCounter._to_bigwig(df_all, tes, chr17_chrom_sizes,
                               output_dir)

    assert {i.name for i in (output_dir / 'counts').iterdir()} == {
        'all_polyA_counts_pos.bw',
        'all_polyA_counts_neg.bw'
    }

    assert {i.name for i in (output_dir / 'coverage').iterdir()} == {
        'all_polyA_coverage_pos.bw',
        'all_polyA_coverage_neg.bw'
    }

    assert {i.name for i in (output_dir / 'ratio').iterdir()} == {
        'all_polyA_ratio_pos.bw',
        'all_polyA_ratio_neg.bw'
    }


def test_TesMultiCounter_samples_read_annot(tmp_path):
    df_read_annot = read_talon_read_annot(read_annot)

    counter = TesMultiCounter(df_read_annot, is_read_annot=True)
    df_all, tes = counter.to_df()

    output_dir = tmp_path / 'read_annot_dir'
    output_dir.mkdir()
    TesMultiCounter._to_bigwig(df_all, tes, chr17_chrom_sizes,
                               output_dir, prefix='tes')

    assert {i.name for i in (output_dir / 'counts').iterdir()} == {
        'all_tes_counts_pos.bw',
        'all_tes_counts_neg.bw',
        'gm12878_tes_counts_neg.bw',
        'gm12878_tes_counts_pos.bw',
        'all_tes_counts_pos.bw',
        'all_tes_counts_neg.bw',
        'hepg2_tes_counts_neg.bw',
        'hepg2_tes_counts_pos.bw'
    }

    assert {i.name for i in (output_dir / 'coverage').iterdir()} == {
        'all_tes_coverage_pos.bw',
        'all_tes_coverage_neg.bw',
        'gm12878_tes_coverage_neg.bw',
        'gm12878_tes_coverage_pos.bw',
        'hepg2_tes_coverage_neg.bw',
        'hepg2_tes_coverage_pos.bw'
    }

    assert {i.name for i in (output_dir / 'ratio').iterdir()} == {
        'all_tes_ratio_pos.bw',
        'all_tes_ratio_neg.bw',
        'gm12878_tes_ratio_neg.bw',
        'gm12878_tes_ratio_pos.bw',
        'hepg2_tes_ratio_neg.bw',
        'hepg2_tes_ratio_pos.bw'
    }


def test_TssMultiCounter_samples(tmp_path):
    # Test read_annot
    df_read_annot = read_talon_read_annot(read_annot)

    counter = TssMultiCounter(df_read_annot, is_read_annot=True)
    df_all, tes = counter.to_df()

    output_dir = tmp_path / 'read_annot_dir'
    output_dir.mkdir()
    TssMultiCounter._to_bigwig(df_all, tes, chr17_chrom_sizes,
                               output_dir, prefix='tss')

    assert {i.name for i in (output_dir / 'counts').iterdir()} == {
        'all_tss_counts_pos.bw',
        'all_tss_counts_neg.bw',
        'gm12878_tss_counts_neg.bw',
        'gm12878_tss_counts_pos.bw',
        'hepg2_tss_counts_neg.bw',
        'hepg2_tss_counts_pos.bw'
    }

    assert {i.name for i in (output_dir / 'coverage').iterdir()} == {
        'all_tss_coverage_pos.bw',
        'all_tss_coverage_neg.bw',
        'gm12878_tss_coverage_neg.bw',
        'gm12878_tss_coverage_pos.bw',
        'hepg2_tss_coverage_neg.bw',
        'hepg2_tss_coverage_pos.bw'
    }

    assert {i.name for i in (output_dir / 'ratio').iterdir()} == {
        'all_tss_ratio_pos.bw',
        'all_tss_ratio_neg.bw',
        'gm12878_tss_ratio_neg.bw',
        'gm12878_tss_ratio_pos.bw',
        'hepg2_tss_ratio_neg.bw',
        'hepg2_tss_ratio_pos.bw'

    }
